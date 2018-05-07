# ("dplyr")
# ("qcc")
# ("randtests")
# ("changepoint")

# Calculate Moving Range for SPC
Moving.Range <- function(inset, bygrp="lvl", xvar="dpmo", dtvar="dt.ym"){
  t.set <- data.frame(inset)
  do.call(c, by(t.set, t.set[, bygrp], function(x){
    hold <- x[[xvar]][order(x[[dtvar]])]
    if(length(hold) <= 12){
      c(rep(NA, length(hold)))
    } else{ # Must have > 12 months history to calculate
      c(rep(NA, 12), sapply(seq(13, length(hold)), function(y) {
        qcc::sd.xbar.one(hold[(y - 12):(y - 1)])
      }))
    }
  }))
}

# Calculate Rolling Average for SPC
Rolling.Average <- function(inset, bygrp="lvl", xvar="dpmo", dtvar="dt.ym"){
  t.set <- data.frame(inset)
  t.set$ra <- do.call(c, by(t.set, t.set[, bygrp], function(x){
    hold <- x[[xvar]][order(x[[dtvar]])]
    if(length(hold) <= 12){
      c(rep(NA, length(hold)))
    } else{ # Must have > 12 months history to calculate
      rm <- as.integer(zoo::rollmean(hold, 12, fill=NA, align="right"))
      c(NA, rm[1:(length(rm) - 1)])
    }
  }))
}

# Calculate SPC Special
SPC.Sp <- function(inset, bylvl="code", family="ff"){
  # Set the months at which UCL evaluation occurs
  t.mons <- seq(min(inset$dt.ym), max(inset$dt.ym), by="months")
  inset$updateUCL <- inset$dt.ym %in%
    t.mons[which(months(t.mons) %in% c("June", "December"))]
  # Set the family and code/cluster levels to trend by
  inset$fam <- inset[[family]]
  inset$lvl <- inset[[bylvl]]
  t.clus <- unique(select(inset, fam, lvl)) %>% filter(!is.na(lvl))
  # Compute trending for every value of the unique level
  t.out <- data.frame()
  for(j in 1:nrow(t.clus)){
    t.code <- filter(inset, fam == t.clus$fam[j] & lvl == t.clus$lvl[j] &
                       !is.na(dpmo))
    # Add in moving ranges and rolling averages
    t.code$mr <- Moving.Range(t.code)
    t.code$ra <- Rolling.Average(t.code)
    # Reduce current month range
    t.mons <- seq(min(t.code$dt.ym), max(t.code$dt.ym), by="months")
    # Assess when evaluation is possible based on 12 months from 1st non-zero DPMO
    t.code$canEval <- as.vector(sapply(unique(t.code$lvl), function(x){
      dpmo <- t.code[t.code$lvl == x, "dpmo"]
      if (all(dpmo == 0)){
        rep(F, length(t.mons))
      } else if (is.na(first(which(dpmo > 0)) + 11 > length(t.mons)) |
          first(which(dpmo > 0)) + 11 > length(t.mons)) {
        rep(F, length(t.mons))
      } else {index(dpmo) >= first(which(dpmo > 0)) + 11}
    }))
    # Based on when evaulation occurs and when it is possible, determine when UCL can be set
    t.code$canSet <- with(t.code, updateUCL & canEval)
    # Calculate UCL as 3-Sigma
    t.code$ucl <- with(t.code, ifelse(canSet, ra + 3 * mr, NA))
    # Carry forward the every-6-month UCLs, then calculate the minimum all-time UCL
    # Also include number of sigmas over the mean for ranking
    t.code <- t.code %>% group_by(lvl) %>%
      mutate(ucl.cf=na.locf(ucl, na.rm=F),
             ucl.min=c(rep(NA, sum(is.na(ucl.cf))),
                       cummin(ucl.cf[!is.na(ucl.cf)])),
             n.sigma=as.numeric(c(rep(NA, length(dpmo) - 1),
                                  ifelse(length(which.min(ucl.min)) == 0, NA,
                                         (dpmo[length(dpmo)] -
                                            ra[which.min(ucl.min)]) /
                                           ifelse(mr[which.min(ucl.min)] == 0, 1,
                                                  mr[which.min(ucl.min)]))))) %>%
      ungroup()
    # Finally, assess for breach
    t.code$breach <- with(t.code, ifelse(dpmo > ucl.min, "Breached UCL", "OK"))
    t.out <- rbind(t.out, select(t.code, -fam, -lvl))
  }
  return(t.out)
}

# Calculate additional trending algorithms
More.Trending <- function(inset, bylvl="code", family="ff"){
  # Set the family and code/cluster levels to trend by
  inset$fam <- inset[[family]]
  inset$lvl <- inset[[bylvl]]
  t.clus <- unique(select(inset, fam, lvl))
  # Compute trending for every value of the unique level
  t.out <- data.frame()
  for(j in 1:nrow(t.clus)){
    t.code <- filter(inset, fam == t.clus$fam[j] & lvl == t.clus$lvl[j])
    t.row <-  cbind.data.frame(t.clus[j, ], dt.ym=dataMonth,
                               tidy(t(Trend.Signals.OC(t.code, dataMonth))))
    t.out <- rbind(t.out, t.row)
  }
  return(merge(inset, t.out, by=c("fam", "lvl", "dt.ym"), all.x=T) %>%
           select(-fam, -lvl))
}


####################################
# STANDARDIZED TRENDING FUNCTIONS
####################################

#################
#' Poisson for Rare Events
#'
#' Returns:
#' Logical of significant trend found
#'
#' Parameter List:
#'
#' inset
#' -----
#' Input data frame with the following columns:
#' dt.ym = Date format with each row being a unique month
#' count = Absolute number of complaints
#'
#' eval.date
#' ---------
#' Which Date-format date in inset$dt.ym to run algorithm on
#'
#' zerorate
#' --------
#' Minimum proportion of months containing zeroes for this algorithm to run
#'
#' prate
#' -----
#' Hypothesized Poisson monthly rate at which the test is perform (null vs. greater)
#'
#' pcrit
#' -----
#' Critical p-value for test
#################

Poisson.Rare <- function(inset, eval.date, zerorate=(2/3), prate=0.2, pcrit=0.05){
  this.mo <- which(inset$dt.ym == eval.date)
  if(length(this.mo) == 0) stop("Date not found")
  # Control period starts at first non-zero count month
  ctrl.period <- inset$count[which.max(inset$count > 0):this.mo]

  if(length(ctrl.period) < 4){
    warning("Poisson did not run - 3 months history required")
    return(NA)
  } else if(sum(ctrl.period != 0) < 2){
    warning("Poisson did not run - At least 2 non-zero months required")
    return(NA)
  } else if(sum(ctrl.period == 0) / length(ctrl.period) >= zerorate){
    pt <- stats::poisson.test(round(sum(ctrl.period)), length(ctrl.period), r=prate,
                       alternative=("greater"))
    if(pt$p.value <= pcrit) {return(TRUE)
      } else{return(FALSE)}
  } else{
    warning("Poisson did not run - Too few zero count months")
    return(NA)

  }
}


# Poisson - Continuous Event Rate Return
Poisson.Rare.C <- function(inset, eval.date, zerorate=(2/3)){
  this.mo <- which(inset$dt.ym == eval.date)
  if(length(this.mo) == 0) stop("Date not found")
  # Control period starts at first non-zero count month
  ctrl.period <- inset$count[which.max(inset$count > 0):this.mo]

  if(length(ctrl.period) < 4){
    warning("Poisson did not run - 3 months history required")
    return(NA)
  } else if(sum(ctrl.period != 0) < 2){
    warning("Poisson did not run - At least 2 non-zero months required")
    return(NA)
  } else if(sum(ctrl.period == 0) / length(ctrl.period) >= zerorate){
    pt <- stats::poisson.test(round(sum(ctrl.period)), length(ctrl.period),
                       alternative=("greater"))
    return(pt$estimate)
  } else{
    warning("Poisson did not run - Too few zero count months")
    return(NA)
  }
}

#################
#' Shewart Signals - 4 Western Electric Rules
#'
#' Returns:
#' Vector of 4 with a logical of whether the Zone Rules 1-4 were met
#'
#' Parameter List:
#'
#' inset
#' -----
#' Input data frame with the following columns:
#' dt.ym = Date format with each row being a unique month
#' dpmo = Calculated complaint rate
#'
#' eval.date
#' ---------
#' Which Date-format date in inset$dt.ym to run algorithm on
#'
#' n.mos
#' -----
#' Number of months to look back when calculating algorithm
#'
#' zerorate
#' --------
#' Maximum proportion of months containing zeroes for this algorithm to run
#################

Shewart <- function(inset, eval.date, n.mos=12, zerorate=(1/3)){
  d2 <- 2 / sqrt(pi)
  if(n.mos < 12) stop("At least 12 months history is required")

  this.mo <- which(inset$dt.ym == eval.date)
  if(length(this.mo) == 0) stop("Date not found")

  if(length(inset$dt.ym[1:this.mo]) >= n.mos){
    ctrl.period <- inset$dpmo[(this.mo - n.mos + 1):this.mo]
    if(sum(ctrl.period == 0) / length(ctrl.period) >= zerorate){
      warning("Shewart did not run - too many zero months")
      return(rep(NA, 4))
    } else{
      mu <- mean(ctrl.period[1:(length(ctrl.period) - 1)])
      sigma <- mean(abs(ctrl.period - lag(ctrl.period))[1:(length(ctrl.period) - 1)],
                    na.rm=T)
      tr1 <- ctrl.period[n.mos] > (mu + 3 * sigma / d2)
      tr2 <- (ctrl.period[n.mos] > (mu + 2 * sigma / d2)) &
        any(ctrl.period[(n.mos - 2):(n.mos - 1)] > (mu + 2 * sigma / d2))
      tr3 <- (ctrl.period[n.mos] > (mu + sigma / d2)) &
        (sum(ctrl.period[(n.mos - 4):(n.mos - 1)] > (mu + sigma / d2)) >= 3)
      tr4 <- sum(ctrl.period[(n.mos - 8):n.mos] > mu) == 9
      return(c(tr1, tr2, tr3, tr4))
    }
  } else{
    warning("Shewart did not run - insufficient months")
    return(rep(NA, 4))
  }
}

# Shewhart - Continuous N-Sigma Return
Shewart.C <- function(inset, eval.date, n.mos=12, zerorate=(1/2)){
  d2 <- 2 / sqrt(pi)

  this.mo <- which(inset$dt.ym == eval.date)
  if(length(this.mo) == 0) stop("Date not found")

  if(length(inset$dt.ym[1:this.mo]) >= n.mos){
    ctrl.period <- inset$dpmo[(this.mo - n.mos + 1):this.mo]
    if(sum(ctrl.period == 0) / length(ctrl.period) >= zerorate){
      warning("Shewart did not run - too many zero months")
      return(NA)
    } else{
      mu <- mean(ctrl.period[1:(length(ctrl.period) - 1)])
      sigma <- mean(abs(ctrl.period - lag(ctrl.period))[1:(length(ctrl.period) - 1)],
                    na.rm=T)
      return((ctrl.period[n.mos] - mu) / (sigma / d2))
    }
  } else{
    warning("Shewart did not run - insufficient months")
    return(NA)
  }
}


#################
#' CUSUM Signals - 4 Western Electric Rules
#'
#' Returns:
#' Vector of 4 with a logical of whether the Zone Rules 1-4 were met
#'
#' Parameter List:
#'
#' inset
#' -----
#' Input data frame with the following columns:
#' dt.ym = Date format with each row being a unique month
#' dpmo = Calculated complaint rate
#'
#' eval.date
#' ---------
#' Which Date-format date in inset$dt.ym to run algorithm on
#'
#' nsigma
#' -----
#' Number of sigmas for which to detect a mean shift in CUSUM
#'
#' zerorate
#' --------
#' Maximum proportion of months containing zeroes for this algorithm to run
#################

CUSUM <- function(inset, eval.date, nsigma=1.5, zerorate=(2/3)){
  d2 <- 2 / sqrt(pi)

  this.mo <- which(inset$dt.ym == eval.date)
  if(length(this.mo) == 0) stop("Date not found")
  # Control period starts at first non-zero DPMO month
  ctrl.period <- inset$dpmo[which.max(inset$dpmo > 0):this.mo]

  if(length(ctrl.period) < 4){
    warning("CUSUM did not run - 3 months history required")
    return(rep(NA, 4))
  } else if(sum(ctrl.period == 0) / length(ctrl.period) >= zerorate){
    warning("CUSUM did not run - too many zero months")
    return(rep(NA, 4))
  } else{
    mu <- mean(ctrl.period[1:(length(ctrl.period) - 1)])
    sigma <- mean(abs(ctrl.period - lag(ctrl.period))[1:(length(ctrl.period) - 1)],
                  na.rm=T)
    ctrl.ctr <- ctrl.period - mu

    # Calculate CUSUM Values
    cusum <- vector()
    for(i in 1:length(ctrl.ctr)){
      if(i==1){
        cusum[1] <- max(0, ctrl.ctr[1] - nsigma * sigma)
      } else{
        cusum[i] <- max(0, ctrl.ctr[i] - nsigma * sigma + cusum[i - 1])
      }
    }

    # Calculate WE Trend Rules
    tr1 <- last(cusum) > (3 * sigma / d2)
    tr2 <- (last(cusum) > (2 * sigma / d2)) &
      any(cusum[(length(cusum) - 2):(length(cusum) - 1)] > (2 * sigma / d2))
    if(length(cusum) >= 5){
      tr3 <- (last(cusum) > (sigma / d2)) &
        (sum(cusum[(length(cusum) - 4):(length(cusum) - 1)] > (sigma / d2)) >= 3)
    } else{tr3 <- NA}
    if(length(cusum) >= 9){
      tr4 <- sum(cusum[(length(cusum) - 8):length(cusum)] > 0) == 9
    } else{tr4 <- NA}
    return(c(tr1, tr2, tr3, tr4))
  }
}


#################
#' EWMA Signals - 4 Western Electric Rules
#'
#' Returns:
#' Vector of 4 with a logical of whether the Zone Rules 1-4 were met
#'
#' Parameter List:
#'
#' inset
#' -----
#' Input data frame with the following columns:
#' dt.ym = Date format with each row being a unique month
#' dpmo = Calculated complaint rate
#'
#' eval.date
#' ---------
#' Which Date-format date in inset$dt.ym to run algorithm on
#'
#' lambda
#' -----
#' Weighting variable
#'
#' zerorate
#' --------
#' Maximum proportion of months containing zeroes for this algorithm to run
#################

EWMA <- function(inset, eval.date, lambda=0.5, zerorate=(3/4)){
  d2 <- 2 / sqrt(pi)

  this.mo <- which(inset$dt.ym == eval.date)
  if(length(this.mo) == 0) stop("Date not found")
  # Control period starts at first non-zero DPMO month
  ctrl.period <- inset$dpmo[which.max(inset$dpmo > 0):this.mo]

  if(length(ctrl.period) < 4){
    warning("EWMA did not run - 3 months history required")
    return(rep(NA, 4))
  } else if(sum(ctrl.period == 0) / length(ctrl.period) >= zerorate){
    warning("EWMA did not run - too many zero months")
    return(rep(NA, 4))
  } else{
    mu <- mean(ctrl.period[1:(length(ctrl.period) - 1)])
    sigma <- mean(abs(ctrl.period - lag(ctrl.period))[1:(length(ctrl.period) - 1)],
                  na.rm=T)
    stderr <- sigma * sqrt(lambda / (2 - lambda))

    # Calculate EWMA Values
    ewma <- vector()
    for(i in 1:length(ctrl.period)){
      if(i==1){
        ewma[1] <- lambda * ctrl.period[1]
      } else{
        ewma[i] <- lambda * ctrl.period[i] + (1 - lambda) * ewma[i -1]
      }
    }

    # Calculate WE Trend Rules
    tr1 <- last(ewma) > mu + (3 * stderr / d2)
    tr2 <- (last(ewma) > (mu + 2 * stderr / d2)) &
      any(ewma[(length(ewma) - 2):(length(ewma) - 1)] > (mu + 2 * stderr / d2))
    if(length(ewma) >= 5){
      tr3 <- (last(ewma) > (mu + stderr / d2)) &
        (sum(ewma[(length(ewma) - 4):(length(ewma) - 1)] > (mu + stderr / d2)) >= 3)
    } else{tr3 <- NA}
    if(length(ewma) >= 9){
      tr4 <- sum(ewma[(length(ewma) - 8):length(ewma)] > mu) == 9
    } else{tr4 <- NA}
    return(c(tr1, tr2, tr3, tr4))
  }
}


# EWMA - Continuous Output
EWMA.C <- function(inset, eval.date, lambda=0.85, zerorate=(1/2)){
  d2 <- 2 / sqrt(pi)

  this.mo <- which(inset$dt.ym == eval.date)
  if(length(this.mo) == 0) stop("Date not found")
  # Control period starts at first non-zero DPMO month
  ctrl.period <- inset$dpmo[which.max(inset$dpmo > 0):this.mo]

  if(length(ctrl.period) < 4){
    warning("EWMA did not run - 3 months history required")
    return(NA)
  } else if(sum(ctrl.period == 0) / length(ctrl.period) >= zerorate){
    warning("EWMA did not run - too many zero months")
    return(NA)
  } else{
    mu <- mean(ctrl.period[1:(length(ctrl.period) - 1)])
    sigma <- mean(abs(ctrl.period - lag(ctrl.period))[1:(length(ctrl.period) - 1)],
                  na.rm=T)
    stderr <- sigma * sqrt(lambda / (2 - lambda))

    # Calculate EWMA Values
    ewma <- vector()
    for(i in 1:length(ctrl.period)){
      if(i==1){
        ewma[1] <- lambda * ctrl.period[1]
      } else{
        ewma[i] <- lambda * ctrl.period[i] + (1 - lambda) * ewma[i -1]
      }
    }

    return((last(ewma) - mu) / (stderr / d2))
  }
}


#################
#' Cox-Stuart Signal
#'
#' Returns:
#' Logical of significant trend found
#'
#' Parameter List:
#'
#' inset
#' -----
#' Input data frame with the following columns:
#' dt.ym = Date format with each row being a unique month
#' dpmo = Calculated complaint rate
#'
#' eval.date
#' ---------
#' Which Date-format date in inset$dt.ym to run algorithm on
#'
#' pcrit
#' -----
#' Critical p-value to evaluate test
#'
#' zerorate
#' --------
#' Maximum proportion of months containing zeroes for this algorithm to run
#################

Cox.Stuart <- function(inset, eval.date, pcrit=0.05, zerorate=(3/4)){
  this.mo <- which(inset$dt.ym == eval.date)
  if(length(this.mo) == 0) stop("Date not found")
  # Control period starts at first non-zero DPMO month
  ctrl.period <- inset$dpmo[which.max(inset$dpmo > 0):this.mo]

  if(length(ctrl.period) < 12){
    warning("Cox-Stuart did not run - 11 months history required")
    return(NA)
  } else if(sum(ctrl.period == 0) / length(ctrl.period) >= zerorate){
    warning("Cox-Stuart did not run - too many zero months")
    return(NA)
  } else{
    test <- randtests::cox.stuart.test(ctrl.period, "r")$p.value <= pcrit
    return(test)
  }
}


#################
#' Forward-Momentum Signal
#'
#' Returns:
#' Logical of significant trend found
#'
#' Parameter List:
#'
#' inset
#' -----
#' Input data frame with the following columns:
#' dt.ym = Date format with each row being a unique month
#' count = Number of complaints
#'
#' eval.date
#' ---------
#' Which Date-format date in inset$dt.ym to run algorithm on
#'
#' pcrit
#' -----
#' Critical p-value to evaluate test
#'
#' zerorate
#' --------
#' Maximum proportion of months containing zeroes for this algorithm to run
#################

Fwd.Momentum <- function(inset, eval.date, pcrit=0.05, zerorate=(3/4)){
  this.mo <- which(inset$dt.ym == eval.date)
  if(length(this.mo) == 0) stop("Date not found")

  ctrl.period <- inset$count[max(1, this.mo - 11):this.mo]

  if(length(ctrl.period) < 12){
    warning("Forward-Momentum did not run - 11 months history required")
    return(NA)
  } else if(sum(ctrl.period == 0) / length(ctrl.period) >= zerorate){
    warning("Forward-Momentum did not run - too many zero months")
    return(NA)
  } else{
    obs <- c(sum(ctrl.period[10:12]), sum(ctrl.period[1:9]))
    exp <- c(round(.25 * sum(ctrl.period)), round(.75 * sum(ctrl.period)))
    pval <- stats::fisher.test(matrix(c(obs, exp), ncol=2, byrow=T),
                        alternative="greater")$p.value
    if(pval <= pcrit){
      return(TRUE)
    } else{
      return(FALSE)
      }
  }
}


#################
#' Changepoint Detection (both uptrend and downtrend)
#'
#' Returns:
#' Logical of significant trend found
#'
#' Parameter List:
#'
#' inset
#' -----
#' Input data frame with the following columns:
#' dt.ym = Date format with each row being a unique month
#' dpmo = Number of complaints
#'
#' eval.date
#' ---------
#' Which Date-format date in inset$dt.ym to run algorithm on
#'
#' pen.drop
#' -----
#' The minimum proportion drop in penalty to qualify a set of changepoints
#'
#' zerorate
#' --------
#' Maximum proportion of months containing zeroes for this algorithm to run
#'
#' n.mos
#' --------
#' Number of months back from the eval.date to detect a changepoint
#################

Changepoint <- function(inset, eval.date, pen.drop=0.2, zerorate=(3/4), n.mos=3){
  this.mo <- which(inset$dt.ym == eval.date)
  if(length(this.mo) == 0) stop("Date not found")
  # Control period starts at first non-zero DPMO month
  ctrl.period <- inset$dpmo[which.max(inset$dpmo > 0):this.mo]

  if(length(ctrl.period) < 4){
    warning("Changepoint did not run - 3 months history required")
    return(NA)
  } else if(sum(ctrl.period == 0) / length(ctrl.period) >= zerorate){
    warning("Changepoint did not run - too many zero months")
    return(NA)
  } else{
    # Add a nominal amount to zero DPMOs to avoid logical errors
    ctrl.period[ctrl.period == 0] <- 1e-5
    # Run BinSeg changepoint
    cpt <- tryCatch(changepoint::cpt.mean(ctrl.period, penalty="Manual", pen.value="2*log(n)",
                    Q=min(9.5, max(length(ctrl.period) / 2 + 1, 1)),
                    method="BinSeg", test.stat="CUSUM"),
                    error=function(x) cpt.mean(rep(10, 10), method="BinSeg"))
    # Find the maximum drop point
    max.drop <- abs(min(pen.value.full(cpt) - dplyr::lag(pen.value.full(cpt)), na.rm=T))
    pen.range <- max(pen.value.full(cpt))
    max.drop.rel <- max.drop / ifelse(pen.range > 0, pen.range, 1e-10)
    drop.lag <- abs(pen.value.full(cpt) - dplyr::lag(pen.value.full(cpt)))
    cpt.opt <- ifelse(max.drop.rel < pen.drop, 0,
                      which(drop.lag == max.drop)[1] - 1)
    # Changepoints
    cpts <- cpt@cpts.full[cpt.opt, 1:cpt.opt] + 1
    # Return the mean difference of a positive changepoint in the last n months
    # If there was one, else return 0
    if (length(cpts) == 0){
      # If there are no changepoints, then 0
      result <- 0
    } else if (max(cpts) > (length(ctrl.period) - n.mos) &
               ctrl.period[max(cpts)] > ctrl.period[max(cpts) - 1]){
      # If there was an up-changepoint in the last n months, then save the
      # mean difference
      # To get the mean difference, find the prior changepoint
      prior <- ifelse(length(cpts) == 1, 1,
                      cpts[order(cpts, decreasing=T)][2])
      result <- mean(ctrl.period[max(cpts):length(ctrl.period)]) -
        mean(ctrl.period[prior:(max(cpts) - 1)])
    } else result <- 0
    return(result)
    # # Return whether the eval.date contains an up-changepoint
    # return(ifelse(length(cpts) == 0, F, (max(cpts) == length(ctrl.period)) &
    #                 (ctrl.period[length(ctrl.period)] >
    #                    ctrl.period[length(ctrl.period) - 1])))
  }
}


#################
#' Uptrend Slope
#'
#' Returns:
#' Slope (dpmo/month) for last n months
#'
#' Parameter List:
#'
#' inset
#' -----
#' Input data frame with the following columns:
#' dt.ym = Date format with each row being a unique month
#' dpmo = Number of complaints
#'
#' eval.date
#' ---------
#' Which Date-format date in inset$dt.ym to run algorithm on
#'
#' n.mos
#' -----
#' The number of months to estimate the slope from
#################

Uptrend <- function(inset, eval.date, n.mos=3){
  this.mo <- which(inset$dt.ym == eval.date)
  if(length(this.mo) == 0) stop("Date not found")
  # Control period starts at first non-zero DPMO month
  ctrl.period <- inset$dpmo[which.max(inset$dpmo > 0):this.mo]
  len <- length(ctrl.period)
  if(len < n.mos){
    warning("Uptrend did not run - Minimum non-zero DPMO months required")
    return(NA)
  } else{
    return(round(lm(ctrl.period[(len - (n.mos - 1)):len] ~ c(1:n.mos))$coefficients[2], 1))
  }
}




#################
#' Return All Trending Signals
#'
#' Example: Trend.Signals(dataset, as.Date("2015-10-01"))
#'
#' Returns:
#' Logical vector of trends by algorithm
#'
#' Parameter List:
#'
#' this
#' -----
#' Input data frame with the following columns:
#' dt.ym = Date format with each row being a unique month
#' dpmo = Calculated complaint rate
#' count = Number of complaints
#'
#' mo
#' ---------
#' Which Date-format date from this$dt.ym to run algorithms on
#################

Trend.Signals <- function(this, mo, warnings=F){
  if(warnings==F){options(warn=-1)}
  sigs <- c(
    structure(Shewart(this, mo), names=c("Shw1", "Shw2", "Shw3", "Shw4")),
    structure(CUSUM(this, mo), names=c("CUS1", "CUS2", "CUS3", "CUS4")),
    structure(EWMA(this, mo), names=c("EWM1", "EWM2", "EWM3", "EWM4")),
    CoxS=Cox.Stuart(this, mo),
    FwdM=Fwd.Momentum(this, mo),
    PoiR=Poisson.Rare(this, mo),
    ChP=Changepoint(this, mo)
  )
  options(warn=0)
  return(sigs)
}

# Optimized Continuous Trend Signals
Trend.Signals.OC <- function(this, mo, warnings=F){
  if(warnings==F){options(warn=-1)}
  sigs <- c(
    Shw=Shewart.C(this, mo),
    EWM=EWMA.C(this, mo),
    PoiR=Poisson.Rare.C(this, mo),
    ChP=Changepoint(this, mo)
  )
  options(warn=0)
  names(sigs) <- c("Shw", "EWM", "PoiR", "ChP")
  return(sigs)
}


#' #################
#' #' Remove the middle numbers of consecutive integers in a vector
#' #################
#'
#' Non.Consec <- function(x){
#'   if(length(x) <= 2){
#'     return(x)
#'   } else{
#'     rle <- rbind(x, c(0, diff(x)))
#'     for(i in 2:(ncol(rle)-1)){
#'       if(rle[2, i] == 1 & rle[2, i+1] == 1){
#'         rle[2, i] <- NA
#'       }
#'     }
#'     x[!is.na(rle[2, ])]
#'   }
#' }


# ======
# Script
# ======

# N/A
