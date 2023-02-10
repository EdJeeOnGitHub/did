

#' Estimate CS DiD with Binary, Absorbing Treatment
#'
#' Estimates CS DiD when Y data is binary and 'absorbing' (once you're vaccinated
#' you cannot become unvaccinated).
#'
#' @param data Dataset to use for estimation
#' @param y_var Integer Y variable that denotes the *first* time period an individual is vaccinated
#' @param group_var Treatment group, in the style of CS package, must be first time period group is treated
#' @param t_var Time period variable.
#' @param id_var Unique ID for individuals.
estimate_did = function(data, 
                        y_var, 
                        group_var, 
                        t_var, 
                        id_var,
                        weight_df = NULL){
    data = as.data.table(data)
    t_levels = sort(unique(data[, get(t_var)]))
    group_levels = data[,.(G = unique(get(group_var)))][order(G), G]
    setkeyv(data, c(id_var, group_var))
    N_indiv_data = create_indiv_per_period_dt(
        df = data,
        group_var = group_var,
        t_var = t_var,
        t_levels = t_levels, 
        group_levels = group_levels,
        weight_df = weight_df
    )
    
    summ_indiv_data = create_indiv_first_treat_dt(
        dt = data,
        y_var = y_var,
        group_var = group_var,
        id_var = id_var
    )
    setkeyv(summ_indiv_data, c(group_var, y_var))
    summ_group_data = create_group_first_treat_dt(
        dt = summ_indiv_data,
        y_var = y_var,
        group_var = group_var,
        t_levels = t_levels,
        group_levels = group_levels,
        weight_df = weight_df
    )
    setkeyv(summ_group_data, c("G", "t"))
    gs_and_ts_we_want = CJ(group = group_levels[group_levels != 0], time = t_levels[t_levels != min(t_levels)])
    # Just put
    att_estimates = map2(
        gs_and_ts_we_want$group,
        gs_and_ts_we_want$time,
        ~calculate_att_g_t(
            g_val = .x,
            t_val = .y,
            lookup_table = summ_group_data,
            N_table = N_indiv_data
        ), 
        .progress = TRUE
    )
    gs_and_ts_we_want[, att_g_t := map_dbl(att_estimates, "att_g_t")]
    gs_and_ts_we_want[, event.time := time - group]
    
    pr_df = N_indiv_data[, .(pr = unique(pr)), G]

    gs_and_ts_we_want = merge(
        gs_and_ts_we_want,
        pr_df,
        by.x = "group",
        by.y = "G",
        all.x = TRUE
    )
    return(gs_and_ts_we_want)
}

create_indiv_first_treat_dt = function(dt, 
                                       y_var, 
                                       group_var, 
                                       id_var,
                                       birth_var = NULL) {
    if (!is.null(birth_var)) {
        summ_dt = dt[
            , 
            .(
                first_Y = unique(get(y_var)),
                G = unique(get(group_var)), 
                born_period = unique(get(birth_var))
            ), 
            by = id_var
            ]
    } else {
        summ_dt = dt[
            , 
            .(
                first_Y = unique(get(y_var)),
                G = unique(get(group_var))
            ), 
            by = id_var
            ]

    }
    summ_dt[, rowid := 1:.N]
    return(summ_dt)
}


create_group_first_treat_dt = function(dt, y_var, group_var, t_levels, group_levels, weight_df = NULL) {
    
    if (is.null(weight_df)) {
        weight_df = data.frame(
            w = 1,
            G = group_levels
        )
    }


    summ_group_dt = dt[, .N, by = c(y_var, group_var)][, t := get(y_var)]
     
    
    zero_t_group_dt = CJ(t = t_levels, group_var = group_levels, n_zero = 0)
    full_group_dt = merge(
        summ_group_dt,
        zero_t_group_dt,
        by.x = c("t", group_var),
        by.y = c("t", "group_var"), 
        all.y = TRUE
    )
    full_group_dt[is.na(N), N := n_zero]
    full_group_dt[, n_zero := NULL]
    full_cumsum_group_dt = full_group_dt[
        order(get(group_var), t),
        .(n = cumsum(N), t = unique(t)),
        group_var
    ]
    setcolorder(full_cumsum_group_dt, c(group_var, "t", "n"))
    full_cumsum_group_dt[, n_lag := shift(n), by = G]
    full_cumsum_group_dt = merge(
        full_cumsum_group_dt,
        weight_df,
        by = "G", 
        all.x = TRUE
    )
    return(full_cumsum_group_dt)
}

create_indiv_per_period_dt = function(df, group_var, t_var, t_levels, group_levels, weight_df = NULL) {
    if (is.null(weight_df)) {
        weight_df = data.frame(
            w = 1,
            G = group_levels
        )
    }
    empty_dt = CJ(G = group_levels, t = t_levels, N_0 = 0 )
    n_indiv = df[, .N, by = c(group_var, t_var)]

    full_dt = merge(
        empty_dt,
        n_indiv,
        all.x = TRUE, 
        by.y = c(group_var, t_var), 
        by.x = c("G", "t")
    )

    full_dt[is.na(N), N := N_0]
    full_dt[, N_0 := NULL]
    full_dt[, N_lag := shift(N), by = group_var]

    full_dt = merge(
        full_dt,
        weight_df,
        by = "G", 
        all.x = TRUE
    )
    pr_df = df[, .(n = .N), by = group_var][, .(pr = n/sum(n), G = get(group_var))]
    full_dt = merge(
        full_dt,
        pr_df,
        by = "G"
    )
    return(full_dt)
}

calculate_att_g_t = function(g_val, 
                             t_val, 
                             lookup_table, 
                             N_table, 
                             verbose = FALSE) {
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


#### Inference ####

rdirichlet = function(n){
    gam_rvar = rgamma(n, 1)
    dir_rvar = gam_rvar / sum(gam_rvar)
    return(dir_rvar)
}

rrademacher = function(n){
    binom_rvar = rbinom(n, 1, 0.5)
    r_rvar = (-1)^binom_rvar
    return(r_rvar)
}




bootstrap_did = function(df, 
                         y_var,
                         group_var,
                         t_var,
                         id_var,
                         cluster_var,
                         B_draws){
    estim_did = estimate_did(
        data = df, 
        y_var = y_var,
        group_var = group_var,
        t_var = t_var,
        id_var = id_var
    )
    boot_draw_df = future_map_dfr(
        1:B_draws,
        ~{
            cluster_df = df[, .(cluster_id = unique(get(cluster_var)))][, .(w = rdirichlet(.N), cluster_id)]
            weight_df = df[
                , 
                .(
                    G = unique(get(group_var)), 
                    cluster_id = unique(get(cluster_var))
                ), 
                by = cluster_var][cluster_df, on = "cluster_id"]
            estimate_did(
                data = df,
                y_var = y_var,
                group_var = group_var,
                t_var = t_var,
                id_var = id_var,
                weight_df = weight_df
            ) %>% mutate(B_draw = .x)

        },
        .progress = TRUE,
        .options = furrr_options(seed = TRUE)
    )

    boot_draw_df = boot_draw_df %>%
        rename(
            boot_att_g_t = att_g_t
        ) %>%
        select(group, time, boot_att_g_t) %>%
        left_join(
            estim_did, 
            by = c("group", "time")
        )

    return(boot_draw_df)

}

calculate_bootstrap_ci = function(boot_df, boot_y_var, y_var, ...){
    
    boot_df %>%
        group_by(...) %>%
        summarise(
            {{ y_var }} := unique({{ y_var }}), 
            conf.low = quantile({{ boot_y_var }}, 0.975, na.rm = TRUE), 
            conf.high = quantile({{ boot_y_var }}, 0.025, na.rm = TRUE)
        ) %>%
        ungroup()


    ## Studentized Bootstrap
    # boot_df %>%
    #     group_by(...) %>%
    #     mutate(
    #         se_hat = sd({{ y_var }}), # can't calculate this cheaply :(
    #         q = ({{ boot_y_var }} - {{ y_var }}) / se_hat
    #         ) %>%
    #     group_by(...) %>%
    #     summarise(
    #         {{ y_var }} := unique({{ y_var }}), 
    #         se_hat = unique(se_hat),
    #         q_975 = quantile(q, 0.975, na.rm = TRUE), 
    #         q_025 = quantile(q, 0.025, na.rm = TRUE)
    #     ) %>%
    #     mutate(
    #         conf.low = att_g_t - q_975*se_hat, 
    #         conf.high = att_g_t - q_025*se_hat
    #     ) %>%
    #     select(
    #         -se_hat, 
    #         -q_975,
    #         -q_025
    #         ) %>%
    #     ungroup()
}


#### Aggregation Stuff ####

calculate_event_study = function(att_pr_df){
    manual_es = att_pr_df %>%
        group_by(
            event.time
        ) %>%
        mutate(
            wt = pr/sum(pr)
        ) %>%
        summarise(
            estimate = sum(wt*att_g_t)
        )
    return(manual_es)
}
