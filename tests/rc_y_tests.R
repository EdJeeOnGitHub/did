
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
    n = 20
)

sim_df = did::build_sim_dataset(
    sp_list = sim_params, 
    panel = TRUE) %>%
    as_tibble()


binary_sim_df = sim_df %>%
    group_by(G) %>%
    mutate(Y_above = Y > quantile(Y, 0.333)) %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(
        first_Y =  min(period[Y_above == TRUE]), 
        first_Y = if_else(!is.finite(first_Y), max(period), first_Y)
    ) %>%
    mutate(Y_binary = period >= first_Y)

# introduce births!
binary_sim_df = binary_sim_df %>%
    group_by(id) %>%
    mutate(
        birth_period = round(runif(1, 0, first_Y))
    )  %>%
    ungroup()


rc_sim_df = binary_sim_df %>%
    dplyr::filter(period >= birth_period)


nrow(rc_sim_df)

df = as.data.table(rc_sim_df)

df[, .(N = uniqueN(id)), birth_period][order(birth_period)]

tictoc::tic()
cs_fit = att_gt(
    data = rc_sim_df,
    yname = "Y_binary",
    tname = "period",
    gname = "G",
    est_method = "ipw",
    idname = "id",
    print_details = TRUE,
    panel = FALSE,
    control_group = "notyettreated" )
tictoc::toc()

tidy_cs_fit = tidy(cs_fit) %>% as_tibble()


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

comp_df

comp_df %>%
    ggplot(aes(
        x = cs_estimate,
        y = manual_estimate
    )) +
    geom_point() +
    geom_abline(linetype = "longdash") +
    theme_bw()

lookup_max_error = comp_df %>%
    dplyr::filter(!(is.na(manual_estimate) & is.na(cs_estimate))) %>%
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

summ_indiv_dt = create_indiv_first_treat_dt(
    dt = df, 
    y_var = "first_Y", 
    group_var = "G", 
    id_var = "id", 
    birth_var = "birth_period"
)

summ_group_dt = create_group_first_treat_dt(
    summ_indiv_dt, 
    "first_Y", 
    "G", 
    c(1:time.periods), 
    unique(summ_indiv_dt$G))
stop()
ipw_did_rc <-function(y, post, D, covariates, i.weights = NULL,
                      boot = FALSE, boot.type = "weighted", nboot = NULL,
                      inffunc = FALSE){
#   D = saved_stuff$G
#   y = saved_stuff$Y 
#   post = saved_stuff$post
    D = ss_G
    y = ss_Y
    post = ss_post
  covariates = NULL
  inffunc = TRUE
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add constant to covariate vector
  int.cov <- as.matrix(rep(1,n))
  if (!is.null(covariates)){
    if(all(as.matrix(covariates)[,1]==rep(1,n))){
      int.cov <- as.matrix(covariates)
    } else {
      int.cov <- as.matrix(cbind(1, covariates))
    }
  }
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Pscore estimation (logit) and also its fitted values
  PS <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights))
  if(PS$converged == FALSE){
    warning(" glm algorithm did not converge")
  }
  if(anyNA(PS$coefficients)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  ps.fit <- as.vector(PS$fitted.values)
  # Do not divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #-----------------------------------------------------------------------------
  #Compute IPW estimator
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post) / (1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/ (1 - ps.fit)

  Pi.hat <- mean(i.weights * D)
  lambda.hat <- mean(i.weights * post)
  one.minus.lambda.hat <- mean(i.weights * (1 - post))

  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * y / (Pi.hat * one.minus.lambda.hat)
  eta.treat.post <- w.treat.post * y / (Pi.hat * lambda.hat)
  eta.cont.pre <- w.cont.pre * y / (Pi.hat * one.minus.lambda.hat)
  eta.cont.post <- w.cont.post * y / (Pi.hat * lambda.hat)



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
  # Asymptotic linear representation of logit's beta's
  score.ps <- i.weights * (D - ps.fit) * int.cov
  Hessian.ps <- stats::vcov(PS) * n
  asy.lin.rep.ps <-  score.ps %*% Hessian.ps
  #-----------------------------------------------------------------------------
  # Influence function of the treated components
  inf.treat.post1 <- eta.treat.post - att.treat.post
  inf.treat.post2 <- -(i.weights * D - Pi.hat) * att.treat.post/Pi.hat
  inf.treat.post3 <- -(i.weights * post - lambda.hat) * att.treat.post/lambda.hat
  inf.treat.post <- inf.treat.post1 + inf.treat.post2 + inf.treat.post3

  inf.treat.pre1 <- eta.treat.pre - att.treat.pre
  inf.treat.pre2 <- -(i.weights * D - Pi.hat) * att.treat.pre/Pi.hat
  inf.treat.pre3 <- -(i.weights * (1 - post) - one.minus.lambda.hat) *
    att.treat.pre/one.minus.lambda.hat
  inf.treat.pret <- inf.treat.pre1 + inf.treat.pre2 + inf.treat.pre3
  # Now, get the influence function of control component
  # First, terms of the inf. funct as if pscore was known
  inf.cont.post1 <- eta.cont.post - att.cont.post
  inf.cont.post2 <- -(i.weights * D - Pi.hat) * att.cont.post/Pi.hat
  inf.cont.post3 <- -(i.weights * post - lambda.hat) * att.cont.post/lambda.hat
  inf.cont.post <- inf.cont.post1 + inf.cont.post2 + inf.cont.post3

  inf.cont.pre1 <- eta.cont.pre - att.cont.pre
  inf.cont.pre2 <- -(i.weights * D - Pi.hat) * att.cont.pre/Pi.hat
  inf.cont.pre3 <- -(i.weights * (1 - post) - one.minus.lambda.hat) *
    att.cont.pre/one.minus.lambda.hat
  inf.cont.pret <- inf.cont.pre1 + inf.cont.pre2 + inf.cont.pre3


    test.inf.treat.post = inf.treat.post
    test.inf.treat.pret = inf.treat.pret
    test.inf.cont.post = inf.cont.post
    test.inf.cont.pret = inf.cont.pret
  # Estimation effect from the propensity score parametes
  # Derivative matrix (k x 1 vector)
  mom.logit.pre <- -eta.cont.pre * int.cov
  mom.logit.pre <- base::colMeans(mom.logit.pre)

  mom.logit.post <- -eta.cont.post * int.cov
  mom.logit.post <- base::colMeans(mom.logit.post)

  # Now the influence function related to estimation effect of pscores
  inf.logit <- asy.lin.rep.ps %*% (mom.logit.post - mom.logit.pre)
  #get the influence function of the DR estimator (put all pieces together)
  att.inf.func <- (inf.treat.post - inf.treat.pret) - (inf.cont.post - inf.cont.pret) +
    inf.logit
  #-----------------------------------------------------------------------------
  if (boot == FALSE) {
    # Estimate of standard error
    se.att <- stats::sd(att.inf.func)/sqrt(n)
    # Estimate of upper boudary of 95% CI
    uci <- ipw.att + 1.96 * se.att
    # Estimate of lower doundary of 95% CI
    lci <- ipw.att - 1.96 * se.att

    #Create this null vector so we can export the bootstrap draws too.
    ipw.boot <- NULL
  }
  if (boot == TRUE) {
    if (is.null(nboot) == TRUE) nboot = 999
    if(boot.type == "multiplier"){
      # do multiplier bootstrap
      ipw.boot <- mboot.did(att.inf.func, nboot)
      # get bootstrap std errors based on IQR
      se.att <- stats::IQR(ipw.boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs(ipw.boot/se.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- ipw.att + cv * se.att
      # Estimate of lower doundary of 95% CI
      lci <- ipw.att - cv * se.att
    } else {
      # do weighted bootstrap
      ipw.boot <- unlist(lapply(1:nboot, wboot_ipw_rc,
                                n = n, y = y, post = post, D = D, int.cov = int.cov, i.weights = i.weights))
      # get bootstrap std errors based on IQR
      se.att <- stats::IQR((ipw.boot - ipw.att)) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs((ipw.boot - ipw.att)/se.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- ipw.att + cv * se.att
      # Estimate of lower doundary of 95% CI
      lci <- ipw.att - cv * se.att

    }

  }
  if(inffunc == FALSE) att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot.type <- ifelse(argu$boot.type=="multiplier", "multiplier", "weighted")
  boot <- ifelse(argu$boot == TRUE, TRUE, FALSE)
  argu <- list(
    panel = FALSE,
    normalized = FALSE,
    boot = boot,
    boot.type = boot.type,
    nboot = nboot,
    type = "ipw"
  )
  ret <- (list(ATT = ipw.att,
               se = se.att,
               uci = uci,
               lci = lci,
               boots = ipw.boot,
               att.inf.func = att.inf.func,
               call.param = call.param,
               argu = argu))
  # Define a new class
  class(ret) <- "drdid"

  # return the list
  return(ret)
}

calculate_rc_influence_function = function(g_val, 
                                           t_val, 
                                           lookup_indiv_table,
                                           verbose = FALSE,
                                           check = FALSE,
                                           prop_score_known = FALSE) {

    g_val = 3
    t_val = 3
    lookup_table = summ_group_dt
    N_table = N_indiv_dt
    lookup_indiv_table = summ_indiv_dt
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
  #-----------------------------------------------------------------------------
  #Pscore estimation (logit) and also its fitted values
  PS <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights))
  if(PS$converged == FALSE){
    warning(" glm algorithm did not converge")
  }
  if(anyNA(PS$coefficients)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  ps.fit <- as.vector(PS$fitted.values)
  # Do not divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #-----------------------------------------------------------------------------
  #Compute IPW estimator
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post) / (1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/ (1 - ps.fit)

  Pi.hat <- mean(i.weights * D)
  lambda.hat <- mean(i.weights * post)
  one.minus.lambda.hat <- mean(i.weights * (1 - post))

  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * y / (Pi.hat * one.minus.lambda.hat)
  eta.treat.post <- w.treat.post * y / (Pi.hat * lambda.hat)
  eta.cont.pre <- w.cont.pre * y / (Pi.hat * one.minus.lambda.hat)
  eta.cont.post <- w.cont.post * y / (Pi.hat * lambda.hat)

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
  # Asymptotic linear representation of logit's beta's
  score.ps <- i.weights * (D - ps.fit) * int.cov
  Hessian.ps <- stats::vcov(PS) * n
  asy.lin.rep.ps <-  score.ps %*% Hessian.ps
  #-----------------------------------------------------------------------------
  # Influence function of the treated components
  inf.treat.post1 <- eta.treat.post - att.treat.post
  inf.treat.post2 <- -(i.weights * D - Pi.hat) * att.treat.post/Pi.hat
  inf.treat.post3 <- -(i.weights * post - lambda.hat) * att.treat.post/lambda.hat
  inf.treat.post <- inf.treat.post1 + inf.treat.post2 + inf.treat.post3

  inf.treat.pre1 <- eta.treat.pre - att.treat.pre
  inf.treat.pre2 <- -(i.weights * D - Pi.hat) * att.treat.pre/Pi.hat
  inf.treat.pre3 <- -(i.weights * (1 - post) - one.minus.lambda.hat) *
    att.treat.pre/one.minus.lambda.hat
  inf.treat.pret <- inf.treat.pre1 + inf.treat.pre2 + inf.treat.pre3

  # Now, get the influence function of control component
  # First, terms of the inf. funct as if pscore was known
  inf.cont.post1 <- eta.cont.post - att.cont.post
  inf.cont.post2 <- -(i.weights * D - Pi.hat) * att.cont.post/Pi.hat
  inf.cont.post3 <- -(i.weights * post - lambda.hat) * att.cont.post/lambda.hat
  inf.cont.post <- inf.cont.post1 + inf.cont.post2 + inf.cont.post3

  inf.cont.pre1 <- eta.cont.pre - att.cont.pre
  inf.cont.pre2 <- -(i.weights * D - Pi.hat) * att.cont.pre/Pi.hat
  inf.cont.pre3 <- -(i.weights * (1 - post) - one.minus.lambda.hat) *
    att.cont.pre/one.minus.lambda.hat
  inf.cont.pret <- inf.cont.pre1 + inf.cont.pre2 + inf.cont.pre3


    # all.equal(inf.treat.post, test.inf.treat.post)
    # all.equal(inf.treat.pret, test.inf.treat.pret)
    # all.equal(inf.cont.pret, test.inf.cont.pret)
    # all.equal(inf.cont.post, test.inf.cont.post)


  # Estimation effect from the propensity score parametes
  # Derivative matrix (k x 1 vector)
  mom.logit.pre <- -eta.cont.pre * int.cov
  mom.logit.pre <- base::colMeans(mom.logit.pre)

  mom.logit.post <- -eta.cont.post * int.cov
  mom.logit.post <- base::colMeans(mom.logit.post)

  # Now the influence function related to estimation effect of pscores
  inf.logit <- asy.lin.rep.ps %*% (mom.logit.post - mom.logit.pre)
  #get the influence function of the DR estimator (put all pieces together)
  att.inf.func <- (inf.treat.post - inf.treat.pret) - (inf.cont.post - inf.cont.pret) +
    inf.logit












    # n = nrow(subset_lookup_indiv_table)
    # if (prop_score_known == FALSE){
    #     PS = stats::glm(D ~ 1, family = "binomial")
    #     ps.fit = as.vector(PS$fitted.value)
    #     w.cont.pre = ps.fit * (1 - D) * (1 - post) / (1 - ps.fit)     
    #     w.cont.post = ps.fit * (1 - D) * post / (1 - ps.fit)

    #     w.treat.pre = D * (1 - post)
    #     w.treat.post = D * post

    # } else {
    #     w.cont.pre = (1 - D) * (1 - post) 
    #     w.cont.post = (1 - D) * post
    # }

    
    # w.treat.pre = D * (1 - post)
    # w.treat.post = D * post

    # Pi.hat = mean(D)
    # lambda.hat = mean(post)
    # one.minus.lambda.hat = mean(1 - post)

    # eta.treat.pre = w.treat.pre * Y / (Pi.hat * one.minus.lambda.hat)
    # eta.treat.post = w.treat.post * Y / (Pi.hat * lambda.hat)

    # eta.cont.pre = w.cont.pre * Y / (Pi.hat * one.minus.lambda.hat)
    # eta.cont.post = w.cont.post * Y / (Pi.hat * lambda.hat)

    # att.treat.pre = mean(eta.treat.pre)
    # att.treat.post = mean(eta.treat.post)
    # att.cont.pre = mean(eta.cont.pre)
    # att.cont.post = mean(eta.cont.post)

    # ipw.att = (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)


    # # Influence function of the treated components
    # inf.treat.post1 <- eta.treat.post - att.treat.post
    # inf.treat.post2 <- -( D - Pi.hat) * att.treat.post/Pi.hat
    # inf.treat.post3 <- -( post - lambda.hat) * att.treat.post/lambda.hat
    # inf.treat.post <- inf.treat.post1 + inf.treat.post2 + inf.treat.post3

    # inf.treat.pre1 <- eta.treat.pre - att.treat.pre
    # inf.treat.pre2 <- -(D - Pi.hat) * att.treat.pre/Pi.hat
    # inf.treat.pre3 <- -((1 - post) - one.minus.lambda.hat) * att.treat.pre/one.minus.lambda.hat
    # inf.treat.pret <- inf.treat.pre1 + inf.treat.pre2 + inf.treat.pre3

    # # Now, get the influence function of control component
    # # First, terms of the inf. funct as if pscore was known
    # inf.cont.post1 <- eta.cont.post - att.cont.post
    # inf.cont.post2 <- -( D - Pi.hat) * att.cont.post/Pi.hat
    # inf.cont.post3 <- -( post - lambda.hat) * att.cont.post/lambda.hat
    # inf.cont.post <- inf.cont.post1 + inf.cont.post2 + inf.cont.post3

    # inf.cont.pre1 <- eta.cont.pre - att.cont.pre
    # inf.cont.pre2 <- -( D - Pi.hat) * att.cont.pre/Pi.hat
    # inf.cont.pre3 <- -( (1 - post) - one.minus.lambda.hat) * att.cont.pre/one.minus.lambda.hat
    # inf.cont.pret <- inf.cont.pre1 + inf.cont.pre2 + inf.cont.pre3


    # mom.logit.pre <- -eta.cont.pre 
    # mom.logit.pre <- mean(mom.logit.pre)

    # mom.logit.post <- -eta.cont.post 
    # mom.logit.post <- mean(mom.logit.post)

    # score.ps <- (D - ps.fit) 
    # Hessian.ps <- stats::vcov(PS) * n
    # asy.lin.rep.ps <-  score.ps %*% Hessian.ps

    # # Now the influence function related to estimation effect of pscores
    # inf.logit <- asy.lin.rep.ps %*% (mom.logit.post - mom.logit.pre)
    # #get the influence function of the DR estimator (put all pieces together)
    # att.inf.func <- (inf.treat.post - inf.treat.pret) - (inf.cont.post - inf.cont.pret) +
    #     inf.logit



    # rowids = which(people_we_want == TRUE, arr.ind = FALSE, useNames = TRUE)
    # names(att.inf.func) = rowids
    
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



    # full_inf_func = matrix(0, nrow(lookup_indiv_table))
    # full_inf_func[rowids] = att.inf.func

    return(lst(g = g_val, t = t_val, att.inf.func))
}

manual_infs = purrr::map2(
    manual_did$g, 
    manual_did$t,
    ~calculate_influence_function(
        g_val = .x, 
        t_val = .y, 
        summ_indiv_dt,
        prop_score_known = TRUE
    )
)

new_infs = map(manual_infs, "full_inf_func")
new_adjustment = map(manual_infs, "n_adjustment")
new_infs = map2(new_infs, new_adjustment, ~.x*.y)
att_inf_output = cs_fit$inffunc
att_infs = map(1:ncol(att_inf_output), ~as.matrix(att_inf_output[, .x]))


test_that("Influence functions match", {
    map2(new_infs, att_infs, all.equal) %>%
    map(expect_true)
}
)


compare_timings = TRUE

if (compare_timings) {
comp_mb = microbenchmark::microbenchmark(
    # cs_fit = att_gt(
    #     data = binary_sim_df,
    #     yname = "Y_binary",
    #     tname = "period",
    #     gname = "G",
    #     est_method = "ipw",
    #     idname = "id",
    #     control_group = "notyettreated" ),
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