

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
                        weight_df = NULL,
                        prop_score_known = FALSE, 
                        biter = 1000, 
                        n_cores = 8){
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

    inf_func_output = map2(
        gs_and_ts_we_want$g, 
        gs_and_ts_we_want$t,
        ~calculate_influence_function(
            g_val = .x, 
            t_val = .y, 
            summ_indiv_data,
            prop_score_known = prop_score_known
        )
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
    gs_and_ts_we_want[, treated := event.time >= 0]
    return(lst(att_df = gs_and_ts_we_want, inf_func_output))
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

calculate_rc_influence_function = function(g_val, 
                                           t_val, 
                                           lookup_indiv_table,
                                           row_id_var,
                                           verbose = FALSE,
                                           check = FALSE,
                                           prop_score_known = FALSE) {

    # g_val = 5
    # t_val = 5
    # lookup_table = summ_group_dt
    # N_table = N_indiv_dt
    # lookup_indiv_table = summ_indiv_dt
    # id_var = "id"
    if (t_val >= g_val) {
        lag_t_val = g_val - 1
    } else {
        lag_t_val = t_val - 1
    }



    people_we_want = lookup_indiv_table[, (G == g_val | (t_val < G | G == 0)) & born_period <= t_val]
    subset_lookup_indiv_table = lookup_indiv_table[people_we_want]
    subset_lookup_indiv_table[, treated := factor(G == g_val, levels = c(TRUE, FALSE))]
    subset_lookup_indiv_table[, Y_post := first_Y <= t_val]
    subset_lookup_indiv_table[, Y_pre := first_Y <= lag_t_val]
    deltaY = subset_lookup_indiv_table[, Y_post - Y_pre]


    Y_post = subset_lookup_indiv_table[, Y_post]
    Y_pre = subset_lookup_indiv_table[born_period < g_val, Y_pre]

    rc_ids = c(subset_lookup_indiv_table[, get(row_id_var)], subset_lookup_indiv_table[born_period < g_val, get(row_id_var)])
    y = c(
        Y_post,
        Y_pre
    )

    post = c(rep(1, length(Y_post)), rep(0, length(Y_pre)))

    D = c(
        subset_lookup_indiv_table[, as.logical(treated)],
        subset_lookup_indiv_table[born_period < g_val, as.logical(treated)]
        )
    n = length(D)
    i.weights = NULL
    int.cov = matrix(1, n)


  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")

  if (!prop_score_known) {
    #Pscore estimation (logit) and also its fitted values
    PS <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights))
    ps.fit <- as.vector(PS$fitted.values)
    # Do not divide by zero
    ps.fit <- pmin(ps.fit, 1 - 1e-16)
    w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
    w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  }
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  #Compute IPW estimator
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  if (prop_score_known) {
    w.cont.pre <- i.weights * (1 - D) * (1 - post)
    w.cont.post <- i.weights * (1 - D) * post
  }

  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * y / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * y / mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * y / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * y / mean(w.cont.post)

  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)

  # ATT estimator
  ipw.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)
  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  if (!prop_score_known) {
    # Asymptotic linear representation of logit's beta's
    score.ps <- i.weights * (D - ps.fit) * int.cov
    Hessian.ps <- stats::vcov(PS) * n
    asy.lin.rep.ps <-  score.ps %*% Hessian.ps
    # Estimation effect from gamma hat (pscore)
    # Derivative matrix (k x 1 vector)
    M2.pre <- base::colMeans(w.cont.pre *(y - att.cont.pre) * int.cov)/mean(w.cont.pre)
    M2.post <- base::colMeans(w.cont.post *(y - att.cont.post) * int.cov)/mean(w.cont.post)
    # Now the influence function related to estimation effect of pscores
    inf.cont.ps <- asy.lin.rep.ps %*% (M2.post - M2.pre)

  } else {
    inf.cont.ps = 0
  }
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function: no estimation effect
  inf.treat.pre <- eta.treat.pre - w.treat.pre * att.treat.pre/mean(w.treat.pre)
  inf.treat.post <- eta.treat.post - w.treat.post * att.treat.post/mean(w.treat.post)
  inf.treat <- inf.treat.post - inf.treat.pre
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect
  inf.cont.pre <- eta.cont.pre - w.cont.pre * att.cont.pre/mean(w.cont.pre)
  inf.cont.post <- eta.cont.post - w.cont.post * att.cont.post/mean(w.cont.post)
  inf.cont <- inf.cont.post - inf.cont.pre


  # Influence function for the control component
  inf.cont <- inf.cont + inf.cont.ps

  #get the influence function of the DR estimator (put all pieces together)
  att.inf.func <- inf.treat - inf.cont


names(att.inf.func) = rc_ids
    
    if (check == TRUE) {

    saved_stuff = read_rds(stringr::str_glue("temp-data/rc-attgt-{g_val - 1}-{t_val - 1}.rds"))

    attgt = saved_stuff$attgt

    saved_stuff$post
    ss_G = c(saved_stuff$G[saved_stuff$post == 1], saved_stuff$G[saved_stuff$post == 0])
    ss_Y = c(saved_stuff$Y[saved_stuff$post == 1], saved_stuff$Y[saved_stuff$post == 0])
    saved_stuff$post
    saved_stuff$w

    mean(saved_stuff$post)

    length(Y_pre)
    length(Y_post)


    mean(post)

    ss_inf_func = c(
        saved_stuff$attgt$att.inf.func[saved_stuff$post == 1],
        saved_stuff$attgt$att.inf.func[saved_stuff$post == 0]
        )
    ss_post = post


    saved_stuff$post
        all.equal(
            ss_G,
            as.numeric(D)) %>%
            print()

        all.equal(
            ss_Y,
           y 
        ) %>%
        print()

        
    length(attgt$att.inf.func)


    bind_cols(
        a = sort(ss_inf_func), b = sort(att.inf.func[, 1])
    )  %>%
    head(20)

        all.equal(
            ss_inf_func, att.inf.func[, 1]
        ) 


    }

    full_inf_func = matrix(0, nrow(lookup_indiv_table))
    full_inf_func[rc_ids] = att.inf.func
    n_all = nrow(lookup_indiv_table)
    n_subset = nrow(subset_lookup_indiv_table)

    return(lst(g = g_val, t = t_val, full_inf_func, n_adjustment = n_all/n_subset))
}

calculate_influence_function = function(g_val, 
                                        t_val, 
                                        lookup_indiv_table,
                                        verbose = FALSE,
                                        check = FALSE,
                                        prop_score_known = FALSE) {

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




    people_we_want = lookup_indiv_table[, G == g_val | (t_val < G | G == 0)]
    subset_lookup_indiv_table = lookup_indiv_table[G == g_val | (t_val < G | G == 0)]
    subset_lookup_indiv_table[, treated := factor(G == g_val, levels = c(TRUE, FALSE))]
    pr_treat = subset_lookup_indiv_table[, mean(treated == TRUE)]
    subset_lookup_indiv_table[, Y_post := first_Y <= t_val]
    subset_lookup_indiv_table[, Y_pre := first_Y <= lag_t_val]
    deltaY = subset_lookup_indiv_table[, Y_post - Y_pre]




    n_all =  nrow(lookup_indiv_table)
    n_subset = nrow(subset_lookup_indiv_table)

    D = subset_lookup_indiv_table[, as.logical(treated)]



    n = nrow(subset_lookup_indiv_table)
    if (prop_score_known == FALSE){
        PS = stats::glm(D ~ 1, family = "binomial")
        ps.fit = as.vector(PS$fitted.value)
        w.cont = ps.fit * (1 - D) / (1 - ps.fit)    
    } else {
        w.cont = (1 - D)
    }

    w.treat = D



    att.treat = w.treat*deltaY
    att.cont = w.cont*deltaY


    eta.treat = mean(att.treat) / mean(w.treat)
    eta.cont = mean(att.cont) / mean(w.cont)

    inf.treat = (att.treat - w.treat * eta.treat) / mean(w.treat)
    inf.cont.1 = (att.cont - w.cont*eta.cont)

    if (prop_score_known == FALSE) {
        score.ps = (D - ps.fit)
        Hessian.ps = stats::vcov(PS) * n
        asy.lin.rep.ps = score.ps %*% Hessian.ps
        M2 = mean(w.cont * (deltaY - eta.cont))
        inf.cont.2 = asy.lin.rep.ps %*% M2
    } else {
        inf.cont.2 = matrix(0, n_subset)
    }

    inf.control = (inf.cont.1 + inf.cont.2) / mean(w.cont)
    att.inf.func = inf.treat - inf.control
    att.inf.func = att.inf.func[, 1]

    rowids = which(people_we_want == TRUE, arr.ind = FALSE, useNames = TRUE)
    names(att.inf.func) = rowids
    
    if (check == TRUE) {

    saved_stuff = read_rds(stringr::str_glue("temp-data/attgt-{g_val - 1}-{t_val - 1}.rds"))

    attgt = saved_stuff$attgt


        all.equal(
            saved_stuff$G,
            subset_lookup_indiv_table[, as.numeric(G == g_val)]) %>%
            print()

        all.equal(
            saved_stuff$Ypost,
            subset_lookup_indiv_table[, Y_post]
        ) %>%
        print()

        all.equal(
            saved_stuff$Ypre,
            subset_lookup_indiv_table[, Y_pre]
        ) %>%
        print()
        

        all.equal(
            saved_stuff$Ypost - saved_stuff$Ypre, 
            as.numeric(deltaY)
        ) %>%
        print()

        all.equal(
            attgt$att.inf.func[, 1], att.inf.func
        ) %>%
        print()


    }



    full_inf_func = matrix(0, nrow(lookup_indiv_table))
    full_inf_func[rowids] = att.inf.func

    return(lst(g = g_val, t = t_val, full_inf_func, n_adjustment = n_all/n_subset))
}

calculate_se = function(inf_matrix, biter = 2000, pl = TRUE, n_cores = 8, alp = 0.05){
    n = nrow(inf_matrix)
    bres = sqrt(n) * run_multiplier_bootstrap(
        inf_matrix, 
        biter, 
        pl = pl, 
        n_cores)
    V = cov(bres)
    bSigma <- apply(bres, 2,
                    function(b) (quantile(b, .75, type=1, na.rm = T) -
                                    quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))
    # critical value for uniform confidence band
    bT <- base::suppressWarnings(apply(bres, 1, function(b) max( abs(b/bSigma), na.rm = T)))
    bT <- bT[is.finite(bT)]
    crit.val <- quantile(bT, 1-alp, type=1, na.rm = T)
    se = as.numeric(bSigma) / sqrt(nrow(inf_matrix))
    return(se)
}
