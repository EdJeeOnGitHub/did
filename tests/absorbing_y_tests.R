


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

source("tests/absorbing_binary_did_functions.R")

ncl <- 1
time.periods <- 4
biters <- 200

# Creates simulation params
sim_params = did::reset.sim(time.periods = 5, n = 20 )

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




cs_fit = att_gt(
    data = binary_sim_df,
    yname = "Y_binary",
    tname = "period",
    gname = "G",
    est_method = "ipw",
    idname = "id",
    control_group = "notyettreated" )

tidy_cs_fit = tidy(cs_fit) %>% as_tibble()

df = as.data.table(binary_sim_df)

manual_did = estimate_did(
    df,
    y_var = "first_Y",
    group_var = "G",
    t_var = "period",
    id_var = "id"
)


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
N_indiv_dt = create_indiv_per_period_dt(df, "G", "period", c(1:5), unique(df$G))
summ_indiv_dt = create_indiv_first_treat_dt(df, "first_Y", "G", "id")
summ_group_dt = create_group_first_treat_dt(
    summ_indiv_dt, 
    "first_Y", 
    "G", 
    c(1:5), 
    unique(summ_indiv_dt$G))


calculate_influence_function = function(g_val, 
                             t_val, 
                             lookup_table, 
                             N_table, 
                             lookup_indiv_table,
                             verbose = FALSE) {

    # g_val = 3
    # t_val = 3
    # lookup_table = summ_group_dt
    # N_table = N_indiv_dt
    # lookup_indiv_table = summ_indiv_dt
    if (t_val >= g_val) {
        lag_t_val = g_val - 1
    } else {
        lag_t_val = t_val - 1
    }


    

    lookup_indiv_table[, treated := factor(G == t_val, levels = c(TRUE, FALSE))]
    pr_treat = lookup_indiv_table[, mean(treated == TRUE)]
    lookup_indiv_table[, switcher := factor(t_val == first_Y, levels = c(TRUE, FALSE))]
    deltaY = lookup_indiv_table[, as.logical(switcher)]

    D = lookup_indiv_table[, as.logical(treated)]

    n = nrow(lookup_indiv_table)

    PS = stats::glm(D ~ 1, family = "binomial")
    ps.fit = as.vector(PS$fitted.value)

    w.treat = D
    w.cont = ps.fit * (1 - D) / (1 - ps.fit)    

    score.ps = (D - ps.fit) 
    Hessian.ps = stats::vcov(PS)*n

    asy.lin.rep.ps = score.ps %*% Hessian.ps


    eta_df = lookup_indiv_table[, table(treated, switcher)]

    eta_t =  eta_df["TRUE", "TRUE"] / sum(eta_df["TRUE", ])
    eta_c = eta_df["FALSE", "TRUE"] / sum(eta_df["FALSE", ])

    att.treat = w.treat*deltaY
    att.cont = w.cont*deltaY

    inf.treat = (att.treat - w.treat * eta_t) / mean(w.treat)
    inf.cont.1 = (att.cont - w.cont*eta_c)


    M2 = mean(w.cont * (deltaY - eta_c))
    inf.cont.2 = asy.lin.rep.ps %*% M2

    inf.control = (inf.cont.1 + inf.cont.2) / mean(w.cont)
    att.inf.func = inf.treat - inf.control

    inf_treat = D*(deltaY - eta_t) / pr_treat

    inf_control =  (1 - D)*(deltaY - eta_c) / (1 - pr_treat) 



    
    inf_func = as.matrix(inf_treat - inf_control)


    y1 = binary_sim_df %>%
        dplyr::filter(period == t_val) %>%
        pull(Y_binary)
    y0 = binary_sim_df %>%
        dplyr::filter(period == lag_t_val) %>%
        pull(Y_binary)

    cs_infunc = DRDID::std_ipw_did_panel(
        y1,
        y0,
        D, 
        covariates = binary_sim_df %>%
            dplyr::filter(period == t_val) %>%
            mutate(const = 1) %>%
            pull(const), 
        inffunc = TRUE
    )$att.inf.func


    print("got here")

    if (verbose == TRUE) {
        return(lst(
            n_g_t_treated, 
            N_g_t,
            n_g_gm1_treated,
            N_g_gm1,
            n_g_t_nyt,
            N_g_t_nyt,
            n_g_gm1_nyt, 
            N_g_gm1_nyt,
            att_g_t
            ))
    }
    return(lst(g = g_val, t = t_val, inf_func, cs_infunc, att.inf.func))
}

binary_sim_df %>% dplyr::filter(G <= 3 | G == 0) %>% dplyr::filter(period == 3) %>% pull(Y_binary)
binary_sim_df %>% dplyr::filter(G <= 3 | G == 0) %>% dplyr::filter(period == 2) %>% pull(Y_binary)

all.equal(binary_sim_df %>% dplyr::filter(period == 3) %>% mutate(G_1 = G == 3) %>% pull(G_1), as.logical(G))
str(attgt)
# this match for some but not others. Investigate ed. Also length doesn't match
manual_did
bind_cols(attgt$att.inf.func, calculate_influence_function(g_val = 3, t_val = 4, lookup_table = summ_group_dt, N_table = N_indiv_dt, summ_indiv_dt)$inf_func)


binary_sim_df
Ypost
Ypre

binary_sim_df

manual_did

devtools::load_all()
cs_fit = att_gt(
    data = binary_sim_df,
    yname = "Y_binary",
    tname = "period",
    gname = "G",
    est_method = "ipw",
    idname = "id",
    control_group = "notyettreated" )

calculate_influence_function(
    g_val = manual_did[16, group], 
    t_val = manual_did[16, time], 
    lookup_table = summ_group_dt,
    N_table = N_indiv_dt,
    summ_indiv_dt
)

manual_infs = purrr::map2(
    manual_did$g, 
    manual_did$t,
    ~calculate_influence_function(
        g_val = .x, 
        t_val = .y, 
        lookup_table = summ_group_dt,
        N_table = N_indiv_dt,
        summ_indiv_dt
    )
)


ed_infs = map(manual_infs, "inf_func")
cs_infs = map(manual_infs, "cs_infunc")
new_infs = map(manual_infs, "att.inf.func")
att_inf_output = cs_fit$inffunc
att_infs = map(1:ncol(att_inf_output), ~as.matrix(att_inf_output[, .x]))

map2(ed_infs, cs_infs, all.equal)
map2(ed_infs, att_infs, all.equal)

manual_did





manual_did


ed_infs[[2]]
att_infs[[2]]

new_infs

bind_cols(a = ed_infs[[1]][, 1], b = att_infs[[1]][, 1]) %>% 
    ggplot(aes(
        x = a, 
        y = b 
    )) +
    geom_abline() +
    geom_point()


cs_fit_se = mboot(att_inf_output, DIDparams = cs_fit$DIDparams)
cs_fit_se$se[[4]]
tidy_cs_fit


att_infs[[1]]

cs_fit$inffunc[, 1]


ed = mboot(inf_M, DIDparams = cs_fit$DIDparams)

ed$se

inf_matrix = map(manual_infs, "inf_func")
inf_M = matrix(unlist(inf_matrix), nrow = length(inf_matrix[[1]]))
inf_M = inf_matrix[[1]]
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

bres = run_multiplier_bootstrap(inf_M, 1000, FALSE, 1)

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