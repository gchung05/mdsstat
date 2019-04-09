#' Mean-Shift Changepoint
#'
#' Test on device-events using the mean-shift changepoint method
#' originally described in Xu, et al 2015.
#'
#' Function \code{cp_mean()} is an implementation of the mean-shift changepoint
#' method originally proposed by Xu, et al (2015) based on testing the
#' mean-centered absolute cumulative sum against a bootstrap null
#' distribution. This algorithm defines a signal as any changepoint found within
#' the last/most recent n=\code{min_seglen} measurements of \code{df}.
#'
#' The parameters in this implementation can be interpreted as
#' follows. Changepoints are detected at an \code{alpha} level based on
#' n=\code{bootstrap_iter} bootstrap iterations (with or without replacement
#' using \code{replace}) of the input time series
#' \code{df}. A minimum of n=\code{min_seglen} consecutive measurements without
#' a changepoint are required to test for an additional changepoint. Both
#' \code{epochs} and \code{cp_max} constrain the maximum possible number of
#' changepoints detectable as follows: within each epoch, each segment of
#' consecutive measurements at least n=\code{min_seglen} measurements long are
#' tested for a changepoint, until no additional changepoints are found.
#'
#' @return A named list of class \code{mdsstat_test} object, as follows:
#' \describe{
#'   \item{test_name}{Name of the test run}
#'   \item{analysis_of}{English description of what was analyzed}
#'   \item{status}{Named boolean of whether the test was run. The name contains
#'   the run status.}
#'   \item{result}{A standardized list of test run results: \code{statistic}
#'   for the test statistic, \code{lcl} and \code{ucl} for the 95%
#'   confidence bounds, \code{p} for the p-value, \code{signal} status, and
#'   \code{signal_threshold}.}
#'   \item{params}{The test parameters}
#'   \item{data}{The data on which the test was run}
#' }
#'
#' @param df Required input data frame of class \code{mds_ts} or, for generic
#' usage, any data frame with the following columns:
#' \describe{
#'   \item{time}{Unique times of class \code{Date}}
#'   \item{event}{Either the event count or rate of class \code{numeric}}
#' }
#' @param ts_event Required if \code{df} is of class \code{mds_ts}. Named string
#' indicating the variable corresponding to the event count or rate. Rate must
#' be calculated in a separate column in \code{df} as it is not calculated by
#' default. The name of the string is an English description of what was
#' analyzed.
#'
#' Default: \code{c("Count"="nA")} corresponding to the event count column in
#' \code{mds_ts} objects. Name is generated from \code{mds_ts} metadata.
#'
#' Example: \code{c("Rate of Bone Filler Events in Canada"="rate")}
#'
#' @param analysis_of Optional string indicating the English description of what
#' was analyzed. If specified, this will override the name of the
#' \code{ts_event} string parameter.
#'
#' Default: \code{NA} indicates no English description for plain \code{df}
#' data frames, or \code{ts_event} English description for \code{df} data frames
#' of class \code{mds_ts}.
#'
#' Example: \code{"Rate of bone cement leakage"}
#'
#' @param eval_period Optional positive integer indicating the number of unique
#' times counting in reverse chronological order to assess. This will be used to
#' establish the process mean and moving range.
#'
#' Default: \code{NULL} considers all times in \code{df}.
#'
#' @param alpha Alpha or Type-I error rate for detection of a changepoint, in
#' the range (0, 1).
#'
#' Default: \code{0.05} detects a changepoint at an alpha level of 0.05 or 5\%.
#'
#' @param cp_max Maximum number of changepoints detectable. This supersedes the
#' theoretical max set by \code{epochs}.
#'
#' Default: \code{100} detects up to a maximum of 100 changepoints.
#'
#' @param min_seglen Minimum required length of consecutive measurements without
#' a changepoint in order to test for an additional changepoint within.
#'
#' Default: \code{6} requires a minimum of 6 consecutive measurements.
#'
#' @param epochs Maximum number of epochs allowed in the iterative search for
#' changepoints, where \code{2^epochs} is the theoretical max changepoints
#' findable. Within each epoch, all measurement segments with a minimum of
#' \code{min_seglen} measurements are tested for a changepoint until no
#' additional changepoints are found.
#'
#' Default: \code{NULL} estimates max epochs from the number of observations or
#' measurements in \code{df} and \code{min_seglen}.
#'
#' @param bootstrap_iter Number of bootstrap iterations for constructing the
#' null distribution of means. Lowest recommended is 1000. Increasing iterations
#' also increases p-value precision.
#'
#' Default: \code{1000} uses 1000 bootstrap iterations.
#'
#' @param replace When sampling for the bootstrap, perform sampling with or
#' without replacement. Unless your \code{df} contains many measurements, and
#' definitely more than \code{bootstrap_iter}, it makes the most sense to set
#' this to \code{TRUE}.
#'
#' Default: \code{T} constructs bootstrap samples with replacement.
#'
#' @param zero_rate Required maximum proportion of \code{event}s in \code{df}
#' (constrained by \code{eval_period}) containing zeroes for this algorithm to
#' run. Because mean-shift changepoint does not perform well on time series with
#' many 0 values, a value >0 is recommended.
#'
#' Default: \code{1/3} requires no more than 1/3 zeros in \code{event}s in
#' \code{df} in order to run.
#'
#' @param ... Further arguments passed onto \code{cp_mean} methods
#'
#' @examples
#' # Basic Example
#' data <- data.frame(time=c(1:25), event=as.integer(stats::rnorm(25, 100, 25)))
#' a1 <- cp_mean(data)
#' # Example using an mds_ts object
#' a2 <- cp_mean(mds_ts[[3]])
#' # Example using a derived rate as the "event"
#' data <- mds_ts[[3]]
#' data$rate <- ifelse(is.na(data$nA), 0, data$nA) / data$exposure
#' a3 <- cp_mean(data, c(Rate="rate"))
#'
#' @references
#' Xu, Zhiheng, et al. "Signal detection using change point analysis in postmarket surveillance." Pharmacoepidemiology and Drug Safety 24.6 (2015): 663-668.
#' @export
cp_mean <- function (df, ...) {
  UseMethod("cp_mean", df)
}

#' @describeIn cp_mean Mean-shift changepoint on mds_ts data
#' @export
cp_mean.mds_ts <- function(
  df,
  ts_event=c("Count"="nA"),
  analysis_of=NA,
  ...
){
  input_param_checker(ts_event, check_names=df, max_length=1)
  if (is.null(names(ts_event))) stop("ts_event must be named")

  # Set NA counts to 0 for "nA" default
  df$nA <- ifelse(is.na(df$nA), 0, df$nA)

  # Set analysis_of
  if (is.na(analysis_of)){
    name <- paste(names(ts_event), "of",
                  paste0(attributes(df)$analysis$device_level_source, " ",
                         attributes(df)$analysis$device_level, ":",
                         attributes(df)$analysis$event_level_source, " ",
                         attributes(df)$analysis$event_level))
  } else name <- analysis_of

  out <- data.frame(time=df$time,
                    event=df[[ts_event]])
  cp_mean.default(out, analysis_of=name, ...)
}

#' @describeIn cp_mean Mean-shift changepoint on general data
#' @export
cp_mean.default <- function(
  df,
  analysis_of=NA,
  eval_period=NULL,
  alpha=0.05,
  cp_max=100,
  min_seglen=6,
  epochs=NULL,
  bootstrap_iter=1000,
  replace=T,
  zero_rate=1/3,
  ...
){
  input_param_checker(df, "data.frame")
  input_param_checker(c("time", "event"), check_names=df)
  input_param_checker(eval_period, "numeric", null_ok=T, max_length=1)
  input_param_checker(alpha, "numeric", null_ok=F, max_length=1)
  input_param_checker(cp_max, "numeric", null_ok=F, max_length=1)
  input_param_checker(min_seglen, "numeric", null_ok=F, max_length=1)
  input_param_checker(epochs, "numeric", null_ok=T, max_length=1)
  input_param_checker(bootstrap_iter, "numeric", null_ok=F, max_length=1)
  input_param_checker(replace, "logical", null_ok=F, max_length=1)
  input_param_checker(zero_rate, "numeric", null_ok=F, max_length=1)
  if (!is.null(eval_period)){
    if (eval_period %% 1 != 0) stop("eval_period must be an integer")
  }
  if (alpha <= 0 | alpha >= 1) stop("alpha must be in range (0, 1)")
  if (cp_max %% 1 != 0) stop("cp_max must be an integer")
  if (cp_max <= 0 ) stop("cp_max must be positive")
  if (min_seglen %% 1 != 0) stop("min_seglen must be an integer")
  if (min_seglen < 4) stop("min_seglen must be >=4")
  if (!is.null(epochs)){
    if (epochs %% 1 != 0) stop("epochs must be an integer")
    if (epochs <= 0 ) stop("epochs must be positive")
  }
  if (bootstrap_iter %% 1 != 0) stop("bootstrap_iter must be an integer")
  if (bootstrap_iter < 1000) stop("bootstrap_iter should be >=1000")
  if (zero_rate < 0 | zero_rate > 1) stop("zero_rate must be in range [0, 1]")

  # Order by time
  df <- df[order(df$time), ]
  # Restrict to eval_period
  if (!is.null(eval_period)){
    if (eval_period > nrow(df)){
      stop("eval_period cannot be greater than df rows")
    } else if (eval_period < 1){
      stop("eval_period must be greater than 0")
    } else df <- df[c((nrow(df) - eval_period + 1):nrow(df)), ]
  } else eval_period <- nrow(df)
  # Return data
  tlen <- nrow(df)
  rd <- list(reference_time=range(df$time),
             data=df)

  # Check for non-runnable conditions
  hyp <- "Not run"
  if(nrow(df) < 4){
    rr <- NA
    rs <- stats::setNames(F, ">3 time periods required")
  } else if(sum(df$time != 0) < 2){
    rr <- NA
    rs <- stats::setNames(F, "2 or more non-zero events required")
  } else if(sum(df$event == 0) / nrow(df) > zero_rate){
    rr <- NA
    rs <- stats::setNames(F, paste("Maximum zero_rate of", zero_rate, "exceeded"))
  } else{
    # If all conditions are met, run CP test

    # Single changepoint based on mean shift
    # Original code by Xu, Kass-Hout 2015
    # Updated by Chung 2019
    # --------------------------------------
    cp_mean_single <- function(
      df,
      bootstrap_iter,
      alpha,
      replace
    ){
      # Changepoint estimate based on max absolute CUSUM difference
      x_ctr <- cumsum(df$event - mean(df$event))
      cp <- which.max(abs(x_ctr)) + 1
      diff <- max(x_ctr) - min(x_ctr)
      # Bootstrap null distribution
      n <- ifelse(length(df$event) > 7, bootstrap_iter,
                  factorial(length(df$event)))
      cpwo <- mat.or.vec(n, 1)
      cpdiff <- mat.or.vec(n, 1)
      for (n1 in 1:n){
        x_random <- sample(df$event, replace=replace)
        x_random_ctr <- cumsum(x_random - mean(x_random))
        cpwo[n1] <- which.max(abs(x_random_ctr)) + 1
        cpdiff[n1] <- max(x_random_ctr) - min(x_random_ctr)
      }
      conflevel <- sum(diff > cpdiff) / n
      # Assemble output
      out <- data.frame(
        cp=df$time[cp],
        cpindex=cp,
        cpfound=(1 - conflevel) <= alpha,

        alpha=alpha,
        pvalue=(1 - conflevel))
      return(out)
    }

    # Changepoint search strategy
    # Inspired by code written by Xu, Kass-Hout 2015
    # Rewritten search start/stop strategy
    # ----------------------------------------------
    stat <- p <- NA
    sig <- F
    # Initial changepoint estimate
    a <- cp_mean_single(
      df=df,
      bootstrap_iter=bootstrap_iter,
      alpha=alpha,
      replace=replace)
    if (a$cpfound){
      cpct <- epochct <- 1
      # Max epochs based on number of measurements
      if (is.null(epochs)) epochs <- floor(log2(nrow(df) / (min_seglen / 2)))
      # Constrain maximum number of findable changepoints
      while (cpct < cp_max & epochct <= epochs){
        cpct_now <- cpct
        a <- a[order(a$cpindex), ]
        # Identify current segments
        cps_now <- c(1, a$cpindex, length(df$event) + 1)
        begin <- cps_now[1]
        bhold <- data.frame()
        # Iterate over current segments
        for (j in c(2:(length(cps_now)))){
          end <- cps_now[j] - 1
          if ((end - begin + 1) >= min_seglen){
            subx <- df[begin:end, ]
            b <- cp_mean_single(
              df=subx,
              bootstrap_iter=bootstrap_iter,
              alpha=alpha,
              replace=replace)
            b$cpindex[1] <- b$cpindex[1] + begin - 1
            if (b$cpfound){
              cpct_now <- cpct_now + 1
              bhold <- rbind(bhold, b)
            }
          }
          begin <- cps_now[j]
        }
        # If no changepoints are found in current epoch, stop search
        if (cpct_now == cpct){
          break
        } else{
          # If this epoch identifies more changepoints than max allowable
          # Keep up to max by ranked p-value
          if (cpct_now > cp_max){
            bhold <- bhold[order(bhold$pvalue), ][c(1:(cp_max - nrow(a))), ]
          }
          a <- rbind(a, bhold)
          cpct <- nrow(a)
          epochct <- epochct + 1
        }
      }
      a <- a[order(a$cpindex, decreasing=T), ]
      # Search for changepoint within the last min_seglen measurements
      cps_found <- a$cpindex[a$cpfound]
      sig_period <- c((tlen - min_seglen + 1):tlen)
      if (any(cps_found %in% sig_period)) sig <- T
      stat <- df$time[cps_found]
      p <- a$pvalue[a$cpfound]
    }

    # Save output parameters
    thresh <- lcl <- ucl <- NA
    hyp <- paste0("CP > bootstrap null CP at alpha=", alpha)
    rr <- list(statistic=stats::setNames(stat, rep("CP time", length(stat))),
               lcl=NA,
               ucl=NA,
               p=stats::setNames(p, rep("CP p-value", length(p))),
               signal=sig,
               signal_threshold=stats::setNames(alpha, "alpha"))
    rs <- stats::setNames(T, "Success")
  }

  # Return test
  out <- list(test_name="Mean-shift changepoint",
              analysis_of=analysis_of,
              status=rs,
              result=rr,
              params=list(test_hyp=hyp,
                          eval_period=eval_period,
                          zero_rate=zero_rate,
                          alpha=alpha,
                          cp_max=cp_max,
                          min_seglen=min_seglen,
                          epochs=epochs,
                          bootstrap_iter=bootstrap_iter,
                          replace=replace),
              data=rd)
  class(out) <- append(class(out), "mdsstat_test")
  return(out)
}
