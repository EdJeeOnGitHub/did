source(here::here("vignettes/setup_sims.R"))
devtools::load_all()
# remotes::install_github("pedrohcgs/DRDID")
# library(did)
library(BMisc)
library(data.table)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(furrr)
library(testthat)
library(data.table)

source("tests/absorbing_binary_did_functions.R")

ncl <- 1
time.periods <- 5
biters <- 200

# Creates simulation params
sim_params = did::reset.sim(
    time.periods = time.periods, 
    n = 2000)

sim_df = did::build_sim_dataset(sp_list = sim_params, panel = TRUE) %>%
    as_tibble()


binary_sim_df = sim_df %>%
    group_by(G) %>%
    mutate(Y_above = Y > quantile(Y, 0.15)) %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(
        first_Y =  min(period[Y_above == TRUE]), 
        first_Y = if_else(!is.finite(first_Y), max(period), first_Y)
    ) %>%
    mutate(Y_binary = period >= first_Y)

df = as.data.table(binary_sim_df)


tictoc::tic()
rc_cs_fit = att_gt(
    data = binary_sim_df,
    yname = "Y_binary",
    tname = "period",
    gname = "G",
    panel = FALSE,
    est_method = "ipw",
    idname = "id",
    print_details = TRUE,
    control_group = "notyettreated" )
tictoc::toc()



tictoc::tic()
cs_fit = att_gt(
    data = binary_sim_df,
    yname = "Y_binary",
    tname = "period",
    gname = "G",
    est_method = "ipw",
    idname = "id",
    print_details = TRUE,
    bstrap = TRUE,
    biter = 10000,
    control_group = "notyettreated" )
tictoc::toc()

tidy_cs_fit = tidy(cs_fit) %>% as_tibble()

inner_join(
rc_cs_fit %>% tidy() %>% as_tibble() %>%
    select(group, time, rc_estimate = estimate, rc_std.error = std.error),
cs_fit %>% tidy() %>% as_tibble() %>%
    select(group, time, cs_estimate = estimate, cs_std.error = std.error), 
    by = c("group", "time")
)


manual_did = estimate_did(
    df,
    y_var = "first_Y",
    group_var = "G",
    t_var = "period",
    id_var = "id"
)$att_df


comp_df = inner_join(
    manual_did %>% rename(manual_estimate = att_g_t),
    tidy_cs_fit %>%
        select(group, time, cs_estimate = estimate), 
    by = c("group","time")
)

# comp_df %>%
#     ggplot(aes(
#         x = cs_estimate,
#         y = manual_estimate
#     )) +
#     geom_point() +
#     geom_abline(linetype = "longdash") +
#     theme_bw()

lookup_max_error = comp_df %>%
    mutate(error = abs(manual_estimate - cs_estimate)) %>%
    summarise(max_error = max(error)) %>%
    pull()

test_that("Manual Estimates OK", {
    expect_lte(lookup_max_error, 1e-8)
}
)


# Profiling Stuff
mb_results = microbenchmark::microbenchmark(
    manual_did = estimate_did(data = df, y_var = "first_Y", group_var = 'G', t_var = "period", id_var = "id" ),
    times = 5
)

mean_time_seconds = mean(mb_results$time / 1e9) 

test_that("Timing okay", {
    expect_lte(mean_time_seconds, 5)
})



#### Influence Function Time ####
N_indiv_dt = create_indiv_per_period_dt(df, "G", "period", c(1:time.periods), unique(df$G))
summ_indiv_dt = create_indiv_first_treat_dt(df, "first_Y", "G", "id")
summ_group_dt = create_group_first_treat_dt(
    summ_indiv_dt, 
    "first_Y", 
    "G", 
    c(1:time.periods), 
    unique(summ_indiv_dt$G))


manual_infs = purrr::map2(
    manual_did$group, 
    manual_did$time,
    ~calculate_influence_function(
        g_val = .x, 
        t_val = .y, 
        summ_indiv_dt,
        prop_score_known = FALSE
    )
)

new_infs = map(manual_infs, "full_inf_func")
new_adjustment = map(manual_infs, "n_adjustment")
new_infs = map2(new_infs, new_adjustment, ~.x*.y)
att_inf_output = cs_fit$inffunc
att_infs = map(1:ncol(att_inf_output), ~as.matrix(att_inf_output[, .x]))


inf_matrix = matrix(unlist(new_infs), nrow = nrow(new_infs[[1]]))
test_that("Influence functions match", {
    map2(new_infs, att_infs, all.equal) %>%
    map(expect_true)
}
)




test_that("Same analytical", {
    ed_V = Matrix::t(inf_matrix) %*% inf_matrix / nrow(inf_matrix)
    diff_df = bind_cols(
        ed_V = ed_V[, 1],
        cs_V = cs_fit$V_analytical[, 1]
    ) %>%
        mutate( 
            diff = ed_V - cs_V
        )
    map(diff_df$diff, expect_lte, 1e-6)
}
)



manual_se = calculate_att_g_t_se(
    inf_matrix,
    biter = 10000
    )

cs_manual_se = calculate_att_g_t_se(
    as.matrix(cs_fit$inffunc),
    biter = 10000
)



comp_se = tibble(
    manual_se = manual_se,
    cs_se = tidy_cs_fit %>%
        pull(std.error), 
    cs_manual_se = cs_manual_se
)


test_that("Homo SEs work", {
    # Results don't tend to match perfectly just due to bootstrapping it seems
    # Even taking cs_fit$inffunc and bootstrapping that tends to give slightly 
    # different results to using inf_matrix.
    se_diff = abs(tidy_cs_fit$std.error - manual_se)/tidy_cs_fit$std.error
    # no worse than 10% deviation from CS package
    map(100*se_diff, expect_lte, 10)
    # On average within 5%
    expect_lte(mean(se_diff), 0.05)

})


#### Aggregations ####

#' @title Compute extra term in influence function due to estimating weights
#'
#' @description A function to compute the extra term that shows up in the
#'  influence function for aggregated treatment effect parameters
#'  due to estimating the weights
#'
#' @param keepers a vector of indices for which group-time average
#'  treatment effects are used to compute a particular aggregated parameter
#' @param pg a vector with same length as total number of group-time average
#'  treatment effects that contains the probability of being in particular group
#' @param weights.ind additional sampling weights (nx1)
#' @param G vector containing which group a unit belongs to (nx1)
#' @param group vector of groups
#'
#' @return nxk influence function matrix
#'
#' @keywords internal
wif <- function(keepers, pg, weights.ind, G, group) {
  # note: weights are all of the form P(G=g|cond)/sum_cond(P(G=g|cond))
  # this is equal to P(G=g)/sum_cond(P(G=g)) which simplifies things here

  # effect of estimating weights in the numerator
  if1 <- sapply(keepers, function(k) {
    (weights.ind * 1*BMisc::TorF(G==group[k]) - pg[k]) /
      sum(pg[keepers])
  })
  # effect of estimating weights in the denominator
  if2 <- base::rowSums( sapply( keepers, function(k) {
    weights.ind*1*BMisc::TorF(G==group[k]) - pg[k]
  })) %*%
    t(pg[keepers]/(sum(pg[keepers])^2))

  # return the influence function for the weights
  if1 - if2
}


calculate_es_estimates = function(att_df, 
                                  inf_matrix,
                                  balance_e = NULL, 
                                  y_var = att_g_t, 
                                  biter = 1000,
                                  n_cores = 8,
                                  prop_score_known = FALSE,
                                  group_vector = NULL
                                  ) {
    # att_df = manual_did
    # biter = 100
    # n_cores = 4
    # prop_score_known = FALSE
    # group_vector = summ_indiv_dt[, G]
    if (!is.null(balance_e)) {
        att_df = att_df %>%
            group_by(group) %>%
            mutate(max_et = mean(event.time)) %>%
            filter(max_et >= balance_e) %>%
            filter(event.time <= balance_e) %>%
            ungroup() %>%
            as.data.table()
    }
    att_df[, wt := pr/sum(pr), .(event.time, treated)]
    agg_pr = att_df[, unique(pr), .(event.time)][, V1]
    es_df = att_df[, .(estimate = sum(get(y_var)*wt)), .(event.time, treated)]

    event_times = unique(es_df$event.time)
    event_time_att_idx = map(
        event_times, 
        ~lst(
            whichones = which(att_df[, event.time == .x]), 
            weights.agg = att_df[event.time == .x, wt],
            pg = att_df[event.time == .x, pr]
            )
        )
    if (prop_score_known == FALSE) {
        weight_if = map(
            event_time_att_idx,
            ~wif(
                keepers = .x$whichones,
                pg = agg_pr,
                weights.ind = rep(1, nrow(inf_matrix)),
                G = group_vector,
                group = att_df[, group]
            )
        )
    } else {
        weight_if = map(event_time_att_idx, NULL)
    }


    et_if = imap(
        event_time_att_idx, 
        ~get_agg_inf_func(
            att = att_df[, get(y_var)], 
            inffunc1 = inf_matrix,
            whichones = .x$whichones,
            weights.agg = .x$weights.agg,
            wif = weight_if[[.y]]
        )
    )
    et_se = map_dbl(
        et_if,
        ~calculate_se(.x, biter = biter, n_cores = n_cores)
    )
    es_df[, std.error := et_se]
    setorder(es_df, event.time)
    return(es_df)
}



manual_es = manual_did %>%
    calculate_es_estimates(
        inf_matrix = inf_matrix, 
        y_var = "att_g_t",
        group_vector = summ_indiv_dt[, G],
        biter = 10000)

tidy_es_fit = cs_fit %>%
    aggte(type = "dynamic", biters = 10000) %>%
    tidy() %>%
    as_tibble()

tidy_es_fit %>%
    select(event.time, estimate, std.error)


test_that("Event Study Matches", {
    es_comp_estimate = bind_cols(
        cs_es = tidy_es_fit$estimate,
        manual_es = manual_es$estimate
    )
    es_comp_se = bind_cols(
        cs_es = tidy_es_fit$std.error,
        manual_es = manual_es$std.error
    )
    map2(
        es_comp_estimate$cs_es,
        es_comp_estimate$manual_es,
        ~expect_equal(.x, .y)
    )
    es_comp_se = es_comp_se %>%
        mutate(diff = abs(manual_es -cs_es)) %>%
        mutate(pct_diff = 100*diff/cs_es)
    map(
        es_comp_se$pct_diff,
        ~expect_lte(.x, 10)
    )

})

compare_timings = TRUE

if (compare_timings) {
comp_mb = microbenchmark::microbenchmark(
    # cs_fit = att_gt(
    #     data = binary_sim_df,
    #     yname = "Y_binary",
    #     tname = "period",
    #     gname = "G",
    #     est_method = "ipw",
    #     bstrap = FALSE,
    #     idname = "id",
    #     control_group = "notyettreated" ),
    just_att = map2(
        manual_did$g,
        manual_did$t,
        ~calculate_att_g_t(
            g_val = .x,
            t_val = .y,
            lookup_table = summ_group_dt,
            N_table = N_indiv_dt 
        )
    ),
    manual_did = estimate_did(
        df,
        y_var = "first_Y",
        group_var = "G",
        t_var = "period",
        id_var = "id"
    ), 
    manual_inf = map2(
        manual_did$g, 
        manual_did$t,
        ~calculate_influence_function(
            g_val = .x, 
            t_val = .y, 
            summ_indiv_dt
        )
    ),
    manual_no_probit_inf = map2(
        manual_did$g, 
        manual_did$t,
        ~calculate_influence_function(
            g_val = .x, 
            t_val = .y, 
            summ_indiv_dt, 
            prop_score_known = TRUE
        )
    ),
    times = 2
)
comp_mb

n_current = manual_did %>%
    nrow()
n_zm = CJ(groups = 1:30, periods = 1:52) %>% nrow()

summ_time_df = comp_mb %>% 
    as_tibble() %>%
    mutate(
        seconds = time / 1e9
    ) %>%
    group_by(
        expr
    ) %>%
    summarise(
        mean_time = mean(seconds)
    ) %>%
    mutate(
        estimated_time_seconds = mean_time*n_zm/n_current, 
        estimated_time_minutes = estimated_time_seconds/60
    )
summ_time_df

}
stop()





str(new_infs)

manual_inf_M = matrix(unlist(new_infs), nrow = nrow(new_infs[[1]]))

# ed = mboot()


profvis::profvis(
    estimate_did(
        df,
        y_var = "first_Y",
        group_var = "G",
        t_var = "period",
        id_var = "id"
    )
)

# profvis::profvis(
#     map2(
#         manual_did$g, 
#         manual_did$t,
#         ~calculate_influence_function(
#             g_val = .x, 
#             t_val = .y, 
#             summ_indiv_dt
#         )
#     )
# )



# cs_fit_se = mboot(att_inf_output, DIDparams = cs_fit$DIDparams)
# cs_fit_se$se[[4]]
# tidy_cs_fit


# att_infs[[1]]

# cs_fit$inffunc[, 1]


# ed = mboot(inf_M, DIDparams = cs_fit$DIDparams)

# ed$se

# inf_matrix = map(manual_infs, "inf_func")
# inf_M = matrix(unlist(inf_matrix), nrow = length(inf_matrix[[1]]))
# inf_M = inf_matrix[[1]]
run_multiplier_bootstrap <- function(inf.func, biters, pl = FALSE, cores = 1) {
  ngroups = ceiling(biters/cores)
  chunks = rep(ngroups, cores)
  # Round down so you end up with the right number of biters
  chunks[1] = chunks[1] + biters - sum(chunks)

  n <- nrow(inf.func)
  parallel.function <- function(biters) {
    BMisc::multiplier_bootstrap(inf.func, biters)
  }
  # From tests, this is about where it becomes worth it to parallelize
  if(n > 2500 & pl == TRUE & cores > 1) {
    results = parallel::mclapply(
      chunks,
      FUN = parallel.function,
      mc.cores = cores
    )
    results = do.call(rbind, results)
  } else {
    results = BMisc::multiplier_bootstrap(inf.func, biters)
  }
  return(results)
}

bres = run_multiplier_bootstrap(matrix(manual_inf_M[, 1]), 1000, pl = TRUE, 8)


V = cov(bres)

bSigma <- apply(bres, 2,
                function(b) (quantile(b, .75, type=1, na.rm = T) -
                                quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))

# critical value for uniform confidence band
bT <- base::suppressWarnings(apply(bres, 1, function(b) max( abs(b/bSigma), na.rm = T)))
bT <- bT[is.finite(bT)]
crit.val <- quantile(bT, 1-alp, type=1, na.rm = T)

se = as.numeric(bSigma) / sqrt(nrow(inf_M))

se

manual_did[, std.error := ed$se]

wide_comp_ci_df = inner_join(
    manual_did %>% rename(
        manual_estimate = att_g_t, 
        manual_std.error = std.error),
    tidy_cs_fit %>%
        select(group, time, cs_estimate = estimate, cs_std.error = std.error), 
    by = c("group","time")
)

wide_comp_ci_df %>%
    ggplot(aes(
        x = cs_std.error, 
        y = manual_std.error
    )) +
    geom_point() +
    geom_abline()

kAse
#se <- rep(0, length(ndg.dim))
se <- rep(NA, length(ndg.dim))
se[ndg.dim] <- as.numeric(bSigma) / sqrt(n_clusters)

stop()


cs_inf_func = cs_fit$inffunc
unique_cs_inf_func = apply(cs_inf_func, 2, unique)

one_we_want = 361

sort(unique(manual_infs[[one_we_want]][[3]]))

sort(unique_cs_inf_func[[one_we_want]])

inf_df

dim(cs_inf_func)

stop()





# plan(multicore, workers = 16)
# manual_boot_att_gt = bootstrap_did(
#     df = df, 
#     y_var = "first_Y", 
#     group_var = 'G', 
#     t_var = "period", 
#     id_var = "id" ,
#     cluster_var = "G",
#     B_draws = 300)

# manual_did_ci = manual_boot_att_gt %>% 
#     calculate_bootstrap_ci(
#         boot_y_var = boot_att_g_t, 
#         y_var = att_g_t, 
#         group, time)


# wide_comp_ci_df = inner_join(
#     manual_did_ci %>% rename(
#         manual_estimate = att_g_t, 
#         manual_conf.low = conf.low, 
#         manual_conf.high = conf.high),
#     tidy_cs_fit %>%
#         select(group, time, cs_estimate = estimate), 
#     by = c("group","time")
# )


# comp_ci_df = bind_rows(
#     manual_did_ci %>% mutate(type = "manual") %>% 
#         rename(estimate = att_g_t), 
#     tidy_cs_fit %>%
#         select(group, time, estimate, conf.low, conf.high) %>%
#         mutate(type = "cs")
# ) 



# long_comp_ci_df = comp_ci_df %>%
#     gather(variable, value, estimate, conf.low, conf.high)  %>%
#     spread(type, value)

# long_comp_ci_df %>%
#     filter(time >= group) %>%
#     ggplot(aes(
#         x = cs, 
#         y = manual
#     )) +
#     geom_point() +
#     facet_wrap(~variable, ncol = 1) +
#     geom_abline() +
#     theme_bw()


# comp_ci_df$group %>% unique()

# comp_ci_df %>%
#     mutate(event.time = time - group) %>%
#     filter(event.time >= 0) %>%
#     filter(group == 15) %>%
#     ggplot(aes(
#         x = event.time, 
#         y = estimate, 
#         ymin = conf.low, 
#         ymax = conf.high, 
#         colour = type
#     )) +
#     geom_pointrange(
#         position = position_dodge(0.5)
#     ) +
#     facet_wrap(~group)


# manual_boot_es = manual_boot_att_gt %>%
#     group_by(B_draw) %>%
#     group_split() %>%
#     map_dfr(calculate_event_study, .id = "B_draw")

# manual_boot_es

# resid_boot_att_gt = merge(
#     manual_boot_att_gt,
#     manual_did[, .(group, time, orig_att_g_t = att_g_t)],
#     all.x = TRUE,
#     by = c("group", "time")
# )
# n_total = df[, uniqueN(id)]
# resid_boot_att_gt[, resid := sqrt(n_total)*(orig_att_g_t - att_g_t)]
# manual_att_se = resid_boot_att_gt[
#     , 
#     .(
#         std.error = (quantile(resid, 0.75) - quantile(resid, 0.25))/(qnorm(0.75) - qnorm(0.25)), 
#         estimate = unique(orig_att_g_t)
#         ), 
#     .(group, time)
#     ]

# bind_rows(
#     manual_att_se %>% mutate(type = "manual"), 
#     tidy_cs_fit %>% mutate(type = "cs")
# ) %>%
#     select( 
#         group, 
#         time, 
#         estimate, 
#         std.error, 
#         type
#     ) %>%
#     select(-estimate) %>%
#     spread(
#         type, std.error
#     )  %>%
#     as_tibble() %>%
#     ggplot(aes(
#         x = cs, 
#         y = manual
#     )) +
#     geom_point() +
#     geom_abline()


# manual_es_se = manual_boot_es %>%
#     group_by(
#         event.time
#     ) %>%
#     summarise(
#         conf.low = quantile(estimate, 0.025, na.rm = TRUE), 
#         conf.high = quantile(estimate, 0.975, na.rm = TRUE)
#     )  %>%
#     left_join(
#         manual_es, 
#         by = "event.time"
#     ) 

# comp_es_se = bind_rows(
#     manual_es_se %>% mutate(type = "manual"),
#     tidy_cs_es %>% 
#         select(estimate, conf.low = point.conf.low, conf.high = point.conf.high, event.time) %>%
#         mutate(type = "cs")
# )


# comp_es_se %>%
#     ggplot(aes(
#         x = event.time, 
#         y = estimate, 
#         ymin = conf.low, 
#         ymax = conf.high, 
#         colour = type
#     )) +
#     geom_pointrange(
#         position = position_dodge(0.5)
#     )
# manual_es_se

# tidy_cs_es

# manual_boot_es %>%
#     group_by(
#         event.time
#     ) %>%
#     summarise(
#         conf.low = quantile(estimate, 0.025), 
#         conf.high = quantile(estimate, 0.975)
#     )  %>%
#     left_join(
#         manual_es
#     ) %>%
#     ggplot(aes(
#         x = event.time, 
#         y = estimate,
#         ymin = conf.low,
#         ymax = conf.high
#     )) +
#     geom_pointrange()


# bind_rows(manual_boot_es) %>%
#     ggplot(aes(
#         x = event.time,
#         y = estimate
#     )) +
#     geom_point


    






# manual_did = estimate_did(data = df, y_var = "first_Y", group_var = 'G', t_var = "period", id_var = "id" )




# ```


# ```{r}



# profvis::profvis(
#     {
#     manual_did = estimate_did(data = df, y_var = "first_Y", group_var = 'G', t_var = "period", id_var = "id" )
#     }
# )


# ```
# cs_es = aggte(cs_fit, type = "dynamic")
# tidy_cs_es = tidy(cs_es) %>% as_tibble()

# comp_es = bind_rows(
#     manual_es %>% mutate(type = "manual"),
#     tidy_cs_es %>% mutate(type = "cs")
# )

# n_es_no_match = comp_es %>%
#     select(
#         event.time, 
#         estimate, 
#         type
#     ) %>%
#     spread(
#         type, estimate
#     ) %>%
#     mutate(diff = abs(cs - manual)) %>%
#     dplyr::filter(diff > 1e-6 ) %>%
#     nrow()

# if (n_es_no_match > 0) {
#     stop("Event study estimates don't match")
# }