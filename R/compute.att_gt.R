
custom_est_known_ps = function (y1, 
                          y0, 
                          D, 
                          covariates, 
                          i.weights = NULL, 
                          boot = FALSE, 
                          boot.type = "weighted", 
                          nboot = NULL, 
                          inffunc = FALSE, 
                          g_val,
                          t_val,
                          prop_score_df,
                          panel_idx
                          ) {
    D <- as.vector(D)
    n <- length(D)
    deltaY <- as.vector(y1 - y0)
    int.cov <- as.matrix(rep(1, n))
    if (!is.null(covariates)) {
        if (all(as.matrix(covariates)[, 1] == rep(1, n))) {
            int.cov <- as.matrix(covariates)
        }
        else {
            int.cov <- as.matrix(cbind(1, covariates))
        }
    }
    if (is.null(i.weights)) {
        i.weights <- as.vector(rep(1, n))
    }
    else if (min(i.weights) < 0) 
        stop("i.weights must be non-negative")
    # PS <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", 
    #     weights = i.weights))
    # ps.fit <- as.vector(PS$fitted.values)
    # ps.fit <- pmin(ps.fit, 1 - 1e-16)

    # calculating propensity score by hand
    subset_condition = g_val == prop_score_df$g_period & t_val == prop_score_df$t_period
    subset_prop_score_df = prop_score_df[subset_condition, ]


    id_idx_pscore = match(panel_idx, subset_prop_score_df$id_var)
    matched_prop_score_df = subset_prop_score_df[id_idx_pscore,]
    ps.fit = matched_prop_score_df$prop_score


    w.treat <- i.weights * D
    w.cont <- i.weights * ps.fit * (1 - D)/(1 - ps.fit)
    att.treat <- w.treat * deltaY
    att.cont <- w.cont * deltaY
    eta.treat <- mean(att.treat)/mean(w.treat)
    eta.cont <- mean(att.cont)/mean(w.cont)
    ipw.att <- eta.treat - eta.cont
    score.ps <- i.weights * (D - ps.fit) * int.cov

    Hessian.ps <-  matrix(0, nrow = ncol(score.ps), ncol = ncol(score.ps))

    asy.lin.rep.ps <- score.ps %*% Hessian.ps
    inf.treat <- (att.treat - w.treat * eta.treat)/mean(w.treat)
    inf.cont.1 <- (att.cont - w.cont * eta.cont)
    M2 <- base::colMeans(w.cont * (deltaY - eta.cont) * int.cov)
    inf.cont.2 <- asy.lin.rep.ps %*% M2
    inf.control <- (inf.cont.1 + inf.cont.2)/mean(w.cont)
    att.inf.func <- inf.treat - inf.control
    if (boot == FALSE) {
        se.att <- stats::sd(att.inf.func)/sqrt(n)
        uci <- ipw.att + 1.96 * se.att
        lci <- ipw.att - 1.96 * se.att
        ipw.boot <- NULL
    }
    if (boot == TRUE) {
        if (is.null(nboot) == TRUE) 
            nboot = 999
        if (boot.type == "multiplier") {
            ipw.boot <- mboot.did(att.inf.func, nboot)
            se.att <- stats::IQR(ipw.boot)/(stats::qnorm(0.75) - 
                stats::qnorm(0.25))
            cv <- stats::quantile(abs(ipw.boot/se.att), probs = 0.95)
            uci <- ipw.att + cv * se.att
            lci <- ipw.att - cv * se.att
        }
        else {
            ipw.boot <- unlist(lapply(1:nboot, wboot.std.ipw.panel, 
                n = n, deltaY = deltaY, D = D, int.cov = int.cov, 
                i.weights = i.weights))
            se.att <- stats::IQR(ipw.boot - ipw.att)/(stats::qnorm(0.75) - 
                stats::qnorm(0.25))
            cv <- stats::quantile(abs((ipw.boot - ipw.att)/se.att), 
                probs = 0.95)
            uci <- ipw.att + cv * se.att
            lci <- ipw.att - cv * se.att
        }
    }
    if (inffunc == FALSE) 
        att.inf.func <- NULL
    call.param <- match.call()
    argu <- mget(names(formals()), sys.frame(sys.nframe()))
    boot.type <- ifelse(argu$boot.type == "multiplier", "multiplier", 
        "weighted")
    boot <- ifelse(argu$boot == TRUE, TRUE, FALSE)
    argu <- list(panel = TRUE, normalized = TRUE, boot = boot, 
        boot.type = boot.type, nboot = nboot, type = "ipw")
    ret <- (list(ATT = ipw.att, se = se.att, uci = uci, lci = lci, 
        boots = ipw.boot, att.inf.func = att.inf.func, call.param = call.param, 
        argu = argu))
    class(ret) <- "drdid"
    return(ret)
}



#' @title Compute Group-Time Average Treatment Effects
#'
#' @description `compute.att_gt` does the main work for computing
#'  multiperiod group-time average treatment effects
#'
#' @param dp A DIDparams object
#'
#' @return a list with length equal to the number of groups times the
#'  number of time periods; each element of the list contains an
#'  object that contains group-time average treatment effect as well
#'  as which group it is for and which time period it is for. It also exports
#'  the influence function which is used externally to compute
#'  standard errors.
#'
#' @keywords internal
#'
#' @export
compute.att_gt <- function(dp) {
  #-----------------------------------------------------------------------------
  # unpack DIDparams
  #-----------------------------------------------------------------------------
  data <- as.data.frame(dp$data)
  yname <- dp$yname
  tname <- dp$tname
  idname <- dp$idname
  xformla <- dp$xformla
  weightsname <- dp$weightsname
  est_method <- dp$est_method
  base_period <- dp$base_period
  panel <- dp$panel
  true_repeated_cross_sections <- dp$true_repeated_cross_sections
  print_details <- dp$print_details
  control_group <- dp$control_group
  anticipation <- dp$anticipation
  gname <- dp$gname
  n  <- dp$n
  nT <- dp$nT
  nG <- dp$nG
  tlist <- dp$tlist
  glist <- dp$glist
  # prop_score_df = expand.grid(
  #   g_period = glist,
  #   t_period = tlist,
  #   id_var = unique(data[,idname])
  # ) %>% data.table() 
  # prop_score_df$prop_score = runif(nrow(prop_score_df), 0, 1)

  prop_score_df = dp$prop_score_df

  #-----------------------------------------------------------------------------
  # main computations
  #-----------------------------------------------------------------------------

  # will populate with all att(g,t)
  attgt.list <- list()

  # place holder in lists
  counter <- 1

  # number of time periods
  tlist.length <- length(tlist)
  tfac <- 0

  if (base_period != "universal") {
    tlist.length <- tlist.length - 1
    tfac <- 1
  }

  # influence function
  inffunc <- Matrix::Matrix(data=0,nrow=n, ncol=nG*(nT-tfac), sparse=TRUE)

  # never treated option
  nevertreated <- (control_group[1] == "nevertreated")

  if(nevertreated) {
    data$.C <- 1*(data[,gname] == 0)
  }

  # rename yname to .y
  data$.y <- data[,yname]

  # loop over groups
  for (g in 1:nG) {

    # Set up .G once
    data$.G <- 1*(data[,gname] == glist[g])

    # loop over time periods
    for (t in 1:tlist.length) {

      #-----------------------------------------------------------------------------
      # Set pret

      # varying base period
      pret <- t

      # universal base period
      if (base_period == "universal") {
        # use same base period as for post-treatment periods
        pret <- tail(which( (tlist+anticipation) < glist[g]),1)
      }

      # use "not yet treated as control"
      # that is, never treated + units that are eventually treated,
      # but not treated by the current period (+ anticipation)
      if(!nevertreated) {
        data$.C <- 1 * ((data[,gname] == 0) |
                            ((data[,gname] > (tlist[max(t,pret)+tfac]+anticipation)) &
                               (data[,gname] != glist[g])))
      }


      # check if in post-treatment period
      if ((glist[g]<=tlist[(t+tfac)])) {

        # update pre-period if in post-treatment period to
        # be  period (g-delta-1)
        pret <- tail(which( (tlist+anticipation) < glist[g]),1)

        # print a warning message if there are no pre-treatment period
        if (length(pret) == 0) {
          warning(paste0("There are no pre-treatment periods for the group first treated at ", glist[g], "\nUnits from this group are dropped"))

          # if there are not pre-treatment periods, code will
          # jump out of this loop
          break
        }
      }


      #-----------------------------------------------------------------------------
      # if we are in period (g-1), normalize results to be equal to 0
      # and break without computing anything
      if (base_period == "universal") {
        if (tlist[pret] == tlist[(t+tfac)]) {
          attgt.list[[counter]] <- list(att=0, group=glist[g], year=tlist[(t+tfac)], post=0)
          inffunc[,counter] <- rep(0,n)
          counter <- counter+1
          next
        }
      }

      # print the details of which iteration we are on
      if (print_details) {
        cat(paste("current period:", tlist[(t+tfac)]), "\n")
        cat(paste("current group:", glist[g]), "\n")
        cat(paste("set pretreatment period to be", tlist[pret]), "\n")
      }

      #-----------------------------------------------------------------------------
      # results for the case with panel data
      #-----------------------------------------------------------------------------

      # post treatment dummy variable
      post.treat <- 1*(glist[g] <= tlist[t+tfac])

      # total number of units (not just included in G or C)
      disdat <- data[data[,tname] == tlist[t+tfac] | data[,tname] == tlist[pret],]


      if (panel) {
        # transform  disdat it into "cross-sectional" data where one of the columns
        # contains the change in the outcome over time.
        disdat <- panel2cs2(disdat, yname, idname, tname, balance_panel=FALSE)

        # still total number of units (not just included in G or C)
        n <- nrow(disdat)

        # pick up the indices for units that will be used to compute ATT(g,t)
        disidx <- disdat$.G==1 | disdat$.C==1

        # pick up the data that will be used to compute ATT(g,t)
        disdat <- disdat[disidx,]

        n1 <- nrow(disdat) # num obs. for computing ATT(g,t)

        # drop missing factors
        disdat <- droplevels(disdat)

        # give short names for data in this iteration
        G <- disdat$.G
        C <- disdat$.C

        # handle pre-treatment universal base period differently
        # we need to do this because panel2cs2 just puts later period
        # in .y1, but if we are in a pre-treatment period with a universal
        # base period, then the "base period" is actually the later period
        Ypre <- if(tlist[(t+tfac)] > tlist[pret]) disdat$.y0 else disdat$.y1
        Ypost <- if(tlist[(t+tfac)] > tlist[pret]) disdat$.y1 else disdat$.y0
        w <- disdat$.w

        # matrix of covariates
        covariates <- model.matrix(xformla, data=disdat)

        #-----------------------------------------------------------------------------
        # more checks for enough observations in each group

        # if using custom estimation method, skip this part
        custom_est_method <- class(est_method) == "function"

        if (!custom_est_method) {
          pscore_problems_likely <- FALSE
          reg_problems_likely <- FALSE

          # checks for pscore based methods
          if (est_method %in% c("dr", "ipw")) {
            preliminary_logit <- glm(G ~ -1 + covariates, family=binomial(link=logit))
            preliminary_pscores <- predict(preliminary_logit, type="response")

          if (print_details) {
            print(paste0("g: ", glist[g], " t: ", tlist[t+tfac], " p_scores: ", unique(preliminary_pscores)))
          }
            

            if (max(preliminary_pscores) >= 0.999) {
              pscore_problems_likely <- TRUE
              warning(paste0("overlap condition violated for ", glist[g], " in time period ", tlist[t+tfac]))
            }
          }

          # check if can run regression using control units
          if (est_method %in% c("dr", "reg")) {
            control_covs <- covariates[G==0,,drop=FALSE]
            #if (determinant(t(control_covs)%*%control_covs, logarithm=FALSE)$modulus < .Machine$double.eps) {
            if ( rcond(t(control_covs)%*%control_covs) < .Machine$double.eps) {
              reg_problems_likely <- TRUE
              warning(paste0("Not enough control units for group ", glist[g], " in time period ", tlist[t+tfac], " to run specified regression"))
            }
          }

          if (reg_problems_likely | pscore_problems_likely) {
            attgt.list[[counter]] <- list(att=NA, group=glist[g], year=tlist[(t+tfac)], post=post.treat)
            inffunc[,counter] <- NA
            counter <- counter+1
            next
          }
        }

        #-----------------------------------------------------------------------------
        # code for actually computing att(g,t)
        #-----------------------------------------------------------------------------

        if (inherits(est_method,"function")) {
          # user-specified function
          attgt <- est_method(y1=Ypost, y0=Ypre,
                              D=G,
                              covariates=covariates,
                              i.weights=w,
                              inffunc=TRUE,
                              g_val = glist[g],
                              t_val = tlist[t],
                              prop_score_df = prop_score_df,
                              panel_idx = disdat[, idname] 
                              )
        } else if (est_method == "ipw") {
          # inverse-probability weights
          attgt <- DRDID::std_ipw_did_panel(Ypost, Ypre, G,
                                            covariates=covariates,
                                            i.weights=w,
                                            boot=FALSE, inffunc=TRUE)
        } else if (est_method == "reg") {
          # regression
          attgt <- DRDID::reg_did_panel(Ypost, Ypre, G,
                                        covariates=covariates,
                                        i.weights=w,
                                        boot=FALSE, inffunc=TRUE)
        } else {
          # doubly robust, this is default
          attgt <- DRDID::drdid_panel(Ypost, Ypre, G,
                                      covariates=covariates,
                                      i.weights=w,
                                      boot=FALSE, inffunc=TRUE)
        }

        # adjust influence function to account for only using
        # subgroup to estimate att(g,t)
        attgt$att.inf.func <- (n/n1)*attgt$att.inf.func

      } else { # repeated cross sections / unbalanced panel

        # pick up the indices for units that will be used to compute ATT(g,t)
        # these conditions are (1) you are observed in the right period and
        # (2) you are in the right group (it is possible to be observed in
        # the right period but still not be part of the treated or control
        # group in that period here
        rightids <- disdat$.rowid[ disdat$.G==1 | disdat$.C==1]

        # this is the fix for unbalanced panels; 2nd criteria shouldn't do anything
        # with true repeated cross sections, but should pick up the right time periods
        # only with unbalanced panel
        disidx <- (data$.rowid %in% rightids) & ( (data[,tname] == tlist[t+tfac]) | (data[,tname]==tlist[pret]))

        # pick up the data that will be used to compute ATT(g,t)
        disdat <- data[disidx,]

        # drop missing factors
        disdat <- droplevels(disdat)

        # give short names for data in this iteration
        G <- disdat$.G
        C <- disdat$.C
        Y <- disdat[,yname]
        post <- 1*(disdat[,tname] == tlist[t+tfac])
        # num obs. for computing ATT(g,t), have to be careful here
        n1 <- sum(G+C)
        w <- disdat$.w

        #-----------------------------------------------------------------------------
        # checks to make sure that we have enough observations
        skip_this_att_gt <- FALSE
        if ( sum(G*post) == 0 ) {
          warning(paste0("No units in group ", glist[g], " in time period ", tlist[t+tfac]))
          skip_this_att_gt <- TRUE
        }
        if ( sum(G*(1-post)) == 0) {
          warning(paste0("No units in group ", glist[g], " in time period ", tlist[t]))
          skip_this_att_gt <- TRUE
        }
        if (sum(C*post) == 0) {
          warning(paste0("No available control units for group ", glist[g], " in time period ", tlist[t+tfac]))
          skip_this_att_gt <- TRUE
        }
        if (sum(C*(1-post)) == 0) {
          warning(paste0("No availabe control units for group ", glist[g], " in time period ", tlist[t]))
          skip_this_att_gt <- TRUE
        }

        if (skip_this_att_gt) {
          attgt.list[[counter]] <- list(att=NA, group=glist[g], year=tlist[(t+tfac)], post=post.treat)
          inffunc[,counter] <- NA
          counter <- counter+1
          next
        }


        # matrix of covariates
        covariates <- model.matrix(xformla, data=disdat)

        #-----------------------------------------------------------------------------
        # code for actually computing att(g,t)
        #-----------------------------------------------------------------------------

        if (inherits(est_method, "function")) {
          # user-specified function
          attgt <- est_method(y=Y,
                              post=post,
                              D=G,
                              covariates=covariates,
                              i.weights=w,
                              inffunc=TRUE, glist[g], tlist[t])
        } else if (est_method == "ipw") {
          # inverse-probability weights
          attgt <- DRDID::std_ipw_did_rc(y=Y,
                                         post=post,
                                         D=G,
                                         covariates=covariates,
                                         i.weights=w,
                                         boot=FALSE, inffunc=TRUE)
        } else if (est_method == "reg") {
          # regression
          attgt <- DRDID::reg_did_rc(y=Y,
                                     post=post,
                                     D=G,
                                     covariates=covariates,
                                     i.weights=w,
                                     boot=FALSE, inffunc=TRUE)
        } else {
          # doubly robust, this is default
          attgt <- DRDID::drdid_rc(y=Y,
                                   post=post,
                                   D=G,
                                   covariates=covariates,
                                   i.weights=w,
                                   boot=FALSE, inffunc=TRUE)
        }

        # n/n1 adjusts for estimating the
        # att_gt only using observations from groups
        # G and C
        attgt$att.inf.func <- (n/n1)*attgt$att.inf.func

        # If ATT is NaN, replace it with NA, and make Influence functions equal to zero
        if(is.nan(attgt$ATT)){
          attgt$ATT <- NA
          attgt$att.inf.func <- 0 * attgt$att.inf.func
        }

      } #end panel if

      # save results for this att(g,t)
      attgt.list[[counter]] <- list(
        att = attgt$ATT, group = glist[g], year = tlist[(t+tfac)], post = post.treat
      )

      # recover the influence function
      # start with vector of 0s because influence function
      # for units that are not in G or C will be equal to 0
      inf.func <- rep(0, n)

      # populate the influence function in the right places
      if(panel) {
        inf.func[disidx] <- attgt$att.inf.func
      } else {
        # aggregate inf functions by id (order by id)
        aggte_inffunc = suppressWarnings(stats::aggregate(attgt$att.inf.func, list(rightids), sum))
        disidx <- (unique(data$.rowid) %in% aggte_inffunc[,1])
        inf.func[disidx] <- aggte_inffunc[,2]
      }


      # save it in influence function matrix
      # inffunc[g,t,] <- inf.func
      inffunc[,counter] <- inf.func

      # update counter
      counter <- counter+1
    } # end looping over t
  } # end looping over g

  return(list(attgt.list=attgt.list, inffunc=inffunc))
}
