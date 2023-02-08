


source(here::here("vignettes/setup_sims.R"))
# devtools::load_all()
# remotes::install_github("pedrohcgs/DRDID")
library(did)
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
sim_params = did::reset.sim(time.periods = 20, n = 1000)

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
    binary_sim_df,
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
N_indiv_dt = create_indiv_per_period_dt(df, "G", "period", c(1:20), unique(df$G))
summ_indiv_dt = create_indiv_first_treat_dt(df, "first_Y", "G", "id")
summ_group_dt = create_group_first_treat_dt(
    summ_indiv_dt, 
    "first_Y", 
    "G", 
    c(1:20), 
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


    # n_g_t_treated = lookup_table[t == t_val & G == g_val, n*w]
    # N_g_t = N_table[t == t_val & G == g_val, N*w]
    # y_g_t_treated  = n_g_t_treated /  N_g_t

    # n_g_gm1_treated = lookup_table[t == lag_t_val & G == g_val, n*w] 
    # N_g_gm1 = N_table[t == lag_t_val & G == g_val, N*w]
    # y_g_gm1_treated = n_g_gm1_treated / N_g_gm1

    # n_g_t_nyt = lookup_table[t == t_val & G != g_val & (t_val < G | G == 0)][, sum(n*w)] 
    # N_g_t_nyt = N_table[t == t_val & G != g_val & (t_val < G | G == 0)][, sum(N*w)]
    # y_g_t_nyt = n_g_t_nyt / N_g_t_nyt


    # n_g_gm1_nyt = lookup_table[t == lag_t_val & G != g_val & (t_val < G | G == 0)][, sum(n*w)] 
    # N_g_gm1_nyt = N_table[t == lag_t_val & G != g_val & (t_val < G | G == 0)][, sum(N*w)]
    # y_g_gm1_nyt = n_g_gm1_nyt / N_g_gm1_nyt

    # att_g_t = (y_g_t_treated - y_g_gm1_treated) - (y_g_t_nyt - y_g_gm1_nyt)
    
    

    lookup_indiv_table[, treated := factor(G == t_val, levels = c(TRUE, FALSE))]
    pr_treat = lookup_indiv_table[, mean(treated == TRUE)]
    lookup_indiv_table[, switcher := factor(t_val == first_Y, levels = c(TRUE, FALSE))]
    deltaY = lookup_indiv_table[, as.logical(switcher)]

    D = lookup_indiv_table[, as.logical(treated)]
    
    eta_df = lookup_indiv_table[, table(treated, switcher)]

    eta_t =  eta_df["TRUE", "TRUE"] / sum(eta_df["TRUE", ])
    eta_c = eta_df["FALSE", "TRUE"] / sum(eta_df["FALSE", ])

    inf_treat = D*(deltaY - eta_t) / pr_treat

    inf_control =  (1 - D)*(deltaY - eta_c) / (1 - pr_treat) 



    
    inf_func = inf_treat - inf_control


    # inf_control
    # unique(inf.control)


    # inf_treat
    # unique(inf.treat)

    # eta_control
    # eta.cont

    # # assume treat_i = 1
    # eta.treat
    # 1 - eta_treat / mean(w.treat)


    # eta_control

    # eta.cont


    # lookup_table[G == g_val & t == t_val]

    # (n_g_t_treated - n_g_gm1_treated) / nrow(test_indiv_dt)

    # sum(deltaY[D == TRUE])


    # sum(deltaY[D == FALSE])

    ##

#     test_indiv_dt = copy(summ_indiv_dt)
#     test_indiv_dt[, treated := G == t_val ]
#     test_indiv_dt[, .N, treated]




#     PS = stats::glm(treated ~ 1, data = test_indiv_dt, family = "binomial")

#     ps.fit = as.vector(PS$fitted.values)
#     ps.fit = pmin(ps.fit, 1 - 1e-16)
#     D = test_indiv_dt[, treated]
    
#     mean(D)

#     w.treat = D
#     n = nrow(test_indiv_dt)
#     w.cont = ( 1 - D ) 

#     deltaY = test_indiv_dt[, switcher := t_val == first_Y][, switcher]

#     sum(deltaY)

#     test_indiv_dt[, switcher := t_val == first_Y]
#     test_indiv_dt[, .N, .(switcher, treated)]

#     mean(att.treat)

#     att.treat = w.treat * deltaY
#     att.cont = w.cont * deltaY

#     eta.treat = mean(att.treat) / mean(w.treat)
#     eta.cont = mean(att.cont) / mean(w.cont)

#     ipw.att = eta.treat - eta.cont


#     score.ps = (D - ps.fit)
#     Hessian.ps = stats::vcov(PS)*n
#     asy.lin.rep.ps = score.ps %*% Hessian.ps

#     inf.treat = (att.treat - w.treat * eta.treat) / mean(w.treat)

#     inf.cont.1 = (att.cont - w.cont * eta.cont)
#     M2 = mean(w.cont * (deltaY - eta.cont))
#     inf.cont.2 = asy.lin.rep.ps %*% M2


#     inf.control = (inf.cont.1 + 0*inf.cont.2) / mean(w.cont)


#     att.inf.func = inf.treat - inf.control

#     enframe(att.inf.func[, 1], value = "ed")  %>%
#         group_by(ed) %>%
#         count()

# y1 = binary_sim_df %>%
#     filter(period == t_val) %>%
#     pull(Y_binary)

# y0 = binary_sim_df %>%
#     filter(period == t_val - 1) %>%
#     pull(Y_binary)

# D_cs = binary_sim_df %>%
#     filter(period == t_val) %>%
#     mutate(treated = t_val == G) %>%
#     pull(treated)



# cs_infunc = DRDID::std_ipw_did_panel(
#     y1,
#     y0,
#     D_cs, 
#     covariates = binary_sim_df %>%
#         filter(period == t_val) %>%
#         mutate(const = 1) %>%
#         pull(const), 
#     inffunc = TRUE
# )$att.inf.func


# all.equal(cs_infunc, att.inf.func)

#    unique(att.inf.func)





#     pr_treated = N_g_t + N_g_t_nyt


#     eta_treat = (y_g_t_treated - y_g_gm1_treated) / pr_treated
#     eta_control = (y_g_t_nyt - y_g_gm1_nyt) / (1 - pr_treated)



#     lookup_table
# eta_treat
    

#     inf_treat = (test_indiv_dt[, switcher*treated] - test_indiv_dt[, treated * eta_treat]) / pr_treated

#     inf_control = (test_indiv_dt[, switcher*(1 - treated)] - test_indiv_dt[, (1 - treated)*eta_control]) / (1 - pr_treated)


    inf_func = inf_treat - inf_control



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
    return(lst(g = g_val, t = t_val, inf_func))
}

manual_did
calculate_influence_function(3, 3, summ_group_dt, N_indiv_dt, summ_indiv_dt)

manual_did[16]

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