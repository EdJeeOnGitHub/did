


source(here::here("vignettes/setup_sims.R"))
devtools::load_all()
# remotes::install_github("pedrohcgs/DRDID")
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
        first_Y =  min(period[Y_above == TRUE])
    ) %>%
    mutate(Y_binary = period >= first_Y)

cs_fit = att_gt(
    binary_sim_df,
    yname = "Y_binary",
    tname = "period",
    gname = "G",
    idname = "id",
    xformla = ~1, 
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
    expect_lte(mean_time_seconds, 4)
})


#### Influence Function Time ####
stop()

N_indiv_dt = create_indiv_per_period_dt(df, "G", "period", c(1:20, Inf), unique(df$G))
summ_indiv_dt = create_indiv_first_treat_dt(df, "first_Y", "G", "id")
summ_group_dt = create_group_first_treat_dt(
    summ_indiv_dt, 
    "first_Y", 
    "G", 
    c(1:20, Inf), 
    unique(summ_indiv_dt$G))



calculate_influence_function = function(g_val, 
                             t_val, 
                             lookup_table, 
                             N_table, 
                             verbose = FALSE) {

    g_val = 3
    t_val = 3
    lookup_table = summ_group_dt
    N_table = N_indiv_dt
    if (t_val >= g_val) {
        lag_t_val = g_val - 1
    } else {
        lag_t_val = t_val - 1
    }


    n_g_t_treated = lookup_table[t == t_val & G == g_val, n*w]
    N_g_t = N_table[t == t_val & G == g_val, N*w]
    y_g_t_treated  = n_g_t_treated /  N_g_t

    n_g_gm1_treated = lookup_table[t == lag_t_val & G == g_val, n*w] 
    N_g_gm1 = N_table[t == lag_t_val & G == g_val, N*w]
    y_g_gm1_treated = n_g_gm1_treated / N_g_gm1

    n_g_t_nyt = lookup_table[t == t_val & G != g_val & (t_val < G | G == 0)][, sum(n*w)] 
    N_g_t_nyt = N_table[t == t_val & G != g_val & (t_val < G | G == 0)][, sum(N*w)]
    y_g_t_nyt = n_g_t_nyt / N_g_t_nyt


    n_g_gm1_nyt = lookup_table[t == lag_t_val & G != g_val & (t_val < G | G == 0)][, sum(n*w)] 
    N_g_gm1_nyt = N_table[t == lag_t_val & G != g_val & (t_val < G | G == 0)][, sum(N*w)]
    y_g_gm1_nyt = n_g_gm1_nyt / N_g_gm1_nyt

    att_g_t = (y_g_t_treated - y_g_gm1_treated) - (y_g_t_nyt - y_g_gm1_nyt)
    

    n_g_t_treated - n_g_gm1_treated

    pr_treated = N_g_t + N_g_t_nyt



    eta_treat = (y_g_t_treated - y_g_gm1_treated) / pr_treated
    eta_control = (y_g_t_nyt - y_g_gm1_nyt) / (1 - pr_treated)

    inf_treat = 


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
    return(lst(g = g_val, t = t_val, att_g_t))
}

calculate_influence_function

G_val = 3
t_val = 3

eta_treat = summ_group_dt[G == G_val & t == t_val, n - n_lag]
eta_control = summ_group_dt[G == G_val & t == ]



plan(multicore, workers = 16)
manual_boot_att_gt = bootstrap_did(
    df = df, 
    y_var = "first_Y", 
    group_var = 'G', 
    t_var = "period", 
    id_var = "id" ,
    cluster_var = "G",
    B_draws = 300)

manual_did_ci = manual_boot_att_gt %>% 
    calculate_bootstrap_ci(
        boot_y_var = boot_att_g_t, 
        y_var = att_g_t, 
        group, time)


wide_comp_ci_df = inner_join(
    manual_did_ci %>% rename(
        manual_estimate = att_g_t, 
        manual_conf.low = conf.low, 
        manual_conf.high = conf.high),
    tidy_cs_fit %>%
        select(group, time, cs_estimate = estimate), 
    by = c("group","time")
)


comp_ci_df = bind_rows(
    manual_did_ci %>% mutate(type = "manual") %>% 
        rename(estimate = att_g_t), 
    tidy_cs_fit %>%
        select(group, time, estimate, conf.low, conf.high) %>%
        mutate(type = "cs")
) 



long_comp_ci_df = comp_ci_df %>%
    gather(variable, value, estimate, conf.low, conf.high)  %>%
    spread(type, value)

long_comp_ci_df %>%
    filter(time >= group) %>%
    ggplot(aes(
        x = cs, 
        y = manual
    )) +
    geom_point() +
    facet_wrap(~variable, ncol = 1) +
    geom_abline() +
    theme_bw()


comp_ci_df$group %>% unique()

comp_ci_df %>%
    mutate(event.time = time - group) %>%
    filter(event.time >= 0) %>%
    filter(group == 15) %>%
    ggplot(aes(
        x = event.time, 
        y = estimate, 
        ymin = conf.low, 
        ymax = conf.high, 
        colour = type
    )) +
    geom_pointrange(
        position = position_dodge(0.5)
    ) +
    facet_wrap(~group)


manual_boot_es = manual_boot_att_gt %>%
    group_by(B_draw) %>%
    group_split() %>%
    map_dfr(calculate_event_study, .id = "B_draw")

manual_boot_es

resid_boot_att_gt = merge(
    manual_boot_att_gt,
    manual_did[, .(group, time, orig_att_g_t = att_g_t)],
    all.x = TRUE,
    by = c("group", "time")
)
n_total = df[, uniqueN(id)]
resid_boot_att_gt[, resid := sqrt(n_total)*(orig_att_g_t - att_g_t)]
manual_att_se = resid_boot_att_gt[
    , 
    .(
        std.error = (quantile(resid, 0.75) - quantile(resid, 0.25))/(qnorm(0.75) - qnorm(0.25)), 
        estimate = unique(orig_att_g_t)
        ), 
    .(group, time)
    ]

bind_rows(
    manual_att_se %>% mutate(type = "manual"), 
    tidy_cs_fit %>% mutate(type = "cs")
) %>%
    select( 
        group, 
        time, 
        estimate, 
        std.error, 
        type
    ) %>%
    select(-estimate) %>%
    spread(
        type, std.error
    )  %>%
    as_tibble() %>%
    ggplot(aes(
        x = cs, 
        y = manual
    )) +
    geom_point() +
    geom_abline()


manual_es_se = manual_boot_es %>%
    group_by(
        event.time
    ) %>%
    summarise(
        conf.low = quantile(estimate, 0.025, na.rm = TRUE), 
        conf.high = quantile(estimate, 0.975, na.rm = TRUE)
    )  %>%
    left_join(
        manual_es, 
        by = "event.time"
    ) 

comp_es_se = bind_rows(
    manual_es_se %>% mutate(type = "manual"),
    tidy_cs_es %>% 
        select(estimate, conf.low = point.conf.low, conf.high = point.conf.high, event.time) %>%
        mutate(type = "cs")
)


comp_es_se %>%
    ggplot(aes(
        x = event.time, 
        y = estimate, 
        ymin = conf.low, 
        ymax = conf.high, 
        colour = type
    )) +
    geom_pointrange(
        position = position_dodge(0.5)
    )
manual_es_se

tidy_cs_es

manual_boot_es %>%
    group_by(
        event.time
    ) %>%
    summarise(
        conf.low = quantile(estimate, 0.025), 
        conf.high = quantile(estimate, 0.975)
    )  %>%
    left_join(
        manual_es
    ) %>%
    ggplot(aes(
        x = event.time, 
        y = estimate,
        ymin = conf.low,
        ymax = conf.high
    )) +
    geom_pointrange()


bind_rows(manual_boot_es) %>%
    ggplot(aes(
        x = event.time,
        y = estimate
    )) +
    geom_point


    






manual_did = estimate_did(data = df, y_var = "first_Y", group_var = 'G', t_var = "period", id_var = "id" )




```


```{r}



profvis::profvis(
    {
    manual_did = estimate_did(data = df, y_var = "first_Y", group_var = 'G', t_var = "period", id_var = "id" )
    }
)


```
cs_es = aggte(cs_fit, type = "dynamic")
tidy_cs_es = tidy(cs_es) %>% as_tibble()

comp_es = bind_rows(
    manual_es %>% mutate(type = "manual"),
    tidy_cs_es %>% mutate(type = "cs")
)

n_es_no_match = comp_es %>%
    select(
        event.time, 
        estimate, 
        type
    ) %>%
    spread(
        type, estimate
    ) %>%
    mutate(diff = abs(cs - manual)) %>%
    dplyr::filter(diff > 1e-6 ) %>%
    nrow()

if (n_es_no_match > 0) {
    stop("Event study estimates don't match")
}