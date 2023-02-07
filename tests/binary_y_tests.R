
source(here::here("vignettes/setup_sims.R"))
devtools::load_all()
library(BMisc)
library(data.table)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(furrr)
library(testthat)

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


binary_sim_df %>%
    group_by(period, G) %>%
    summarise(n_outcome = sum(Y_binary)) %>%
    ggplot(aes(
        x = period, 
        y = n_outcome, 
        colour = factor(G)
    )) +
    geom_line()

cs_fit = att_gt(
    binary_sim_df,
    yname = "Y_binary",
    tname = "period",
    gname = "G",
    idname = "id",
    xformla = ~1, 
    print_details = TRUE,
    control_group = "notyettreated" )

df = as.data.table(binary_sim_df)


calc_expectation = function(df){
    # N.B. will impute first time period - 1 so don't include very first time period (CS drop this automatically) 
    if (nrow(df[Y_binary == TRUE]) == 0){
        expectation = 0
    } else {
        expectation = df[Y_binary == TRUE, n]/sum(df[, n])
    }
    return(expectation)
}

att_g_t = function(df, t, g, verbose = FALSE){
    if (t >= g) {
        y_t_treated = df[G == g & period == t, .(n = .N, g = g, t = g), by = Y_binary]
        y_gm1_treated = df[G == g & period == g - 1, .(n = .N, g = g, t = g), by = Y_binary]

        # are you treated by period t
        y_t_nyt = df[(t < G | G == 0) & period == t, .(n = .N, g = g, t = g), by = Y_binary]
        y_gm1_nyt = df[(t < G | G == 0) & period == g-1, .(n = .N, g = g, t = g), by = Y_binary]
    }
    
    if (t < g) {
        y_t_treated = df[G == g & period == t, .(n = .N, g = g, t = g), by = Y_binary]
        y_gm1_treated = df[G == g & period == t - 1, .(n = .N, g = g, t = g), by = Y_binary]
        
        y_t_nyt = df[period == t & G != g & (G > t | G == 0), .(n = .N, g = g, t = g), by = Y_binary]
        y_gm1_nyt = df[period == t -1 & G != g & (G > t | G == 0), .(n = .N, g = g, t = g), by = Y_binary]
    }


    e_y_t_treated = calc_expectation(y_t_treated)
    e_y_gm1_treated = calc_expectation(y_gm1_treated)

    e_y_t_nyt = calc_expectation(y_t_nyt)
    e_y_gm1_nyt = calc_expectation(y_gm1_nyt)

    treat_diff = e_y_t_treated - e_y_gm1_treated
    control_diff = e_y_t_nyt - e_y_gm1_nyt

    att_g_t = treat_diff - control_diff
    if (verbose == TRUE) {
        return(lst(
            att_g_t, 
            g, 
            t, 
            y_t_treated, 
            y_gm1_treated, 
            y_t_nyt,
            y_gm1_nyt))
    }
    return(lst(att_g_t, g, t))
}

gs_and_ts = tidy(cs_fit) %>%
    select(G = group, period = time) %>%
    unique() %>%
    as_tibble()

    

manual_fit = map2_dfr(
    gs_and_ts$G,
    gs_and_ts$period,
    ~att_g_t(df, g = .x, t = .y)
) %>%
    rename(estimate = att_g_t, group = g, time = t)

tidy_cs_fit = cs_fit %>%
    tidy() %>%
    as_tibble()




wide_comp_df = inner_join(
    tidy_cs_fit %>%
        rename(cs_estimate = estimate), 
    manual_fit  %>% rename(manual_estimate = estimate),
        by = c("group", "time")
)


wide_comp_df %>%
    ggplot(aes(
        x = cs_estimate,
        y = manual_estimate
    )) +
    geom_point() +
    geom_abline(linetype = "longdash") +
    theme_bw()

max_error = wide_comp_df %>%
    mutate(error = abs(manual_estimate - cs_estimate)) %>%
    summarise(max_error = max(error)) %>%
    pull()


if (max_error > 1e-8) {
    stop("Manual Estimates Off")
}

test_that("Binary Manual Estimates OK", {
    expect_lte(max_error, 1e-8)
})



