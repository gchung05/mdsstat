#' Cumuluative Sum (CUSUM)
#'
#' Test on device-events using the tabular CUSUM (CUmulative SUM) method.
#'
#' Function \code{cusum()} is a standard implementation of the tabular CUSUM
#' method originally proposed by E.S. Page. CUSUM is part of the family of
#' statistical process control tests.
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
#' @param delta Required number of sigmas at which to detect a mean shift.
#'
#' Default: \code{1} detects a mean shift of one sigma.
#'
#' @param H Optional positive number representing the decision interval bound.
#' Lower values will result in a more sensitive test.
#'
#' Default: \code{NULL} uses a value of 5 times the estimated sigma.
#'
#' @param zero_rate Required maximum proportion of \code{event}s in \code{df}
#' (constrained by \code{eval_period}) containing zeroes for this algorithm to
#' run. Because CUSUM does not perform well on time series with many 0 values,
#' a value >0 is recommended.
#'
#' Default: \code{1/3} requires no more than 1/3 zeros in \code{event}s in
#' \code{df} in order to run.
#'
#' @param ... Further arguments passed onto \code{cusum} methods
#'
#' @examples
#' # Basic Example
#' data <- data.frame(time=c(1:25), event=as.integer(stats::rnorm(25, 100, 25)))
#' a1 <- cusum(data)
#' # Example using an mds_ts object
#' a2 <- cusum(mds_ts[[3]])
#' # Example using a derived rate as the "event"
#' data <- mds_ts[[3]]
#' data$rate <- ifelse(is.na(data$nA), 0, data$nA) / data$exposure
#' a3 <- cusum(data, c(Rate="rate"))
#'
#' @references
#' Page, E. S. (June 1954). "Continuous Inspection Scheme". Biometrika. 41 (1/2): 100–115. doi:10.1093/biomet/41.1-2.100. JSTOR 2333009.
#' @export
cusum <- function (df, ...) {
  UseMethod("cusum", df)
}

#' @describeIn cusum CUSUM on mds_ts data
#' @export
cusum.mds_ts <- function(
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
  cusum.default(out, analysis_of=name, ...)
}

#' @describeIn cusum CUSUM on general data
#' @export
cusum.default <- function(
  df,
  analysis_of=NA,
  eval_period=NULL,
  delta=1,
  H=NULL,
  zero_rate=1/3,
  ...
){
  input_param_checker(df, "data.frame")
  input_param_checker(c("time", "event"), check_names=df)
  input_param_checker(eval_period, "numeric", null_ok=T, max_length=1)
  input_param_checker(delta, "numeric", null_ok=F, max_length=1)
  input_param_checker(H, "numeric", null_ok=F, max_length=1)
  input_param_checker(zero_rate, "numeric", null_ok=F, max_length=1)
  if (!is.null(eval_period)){
    if (eval_period %% 1 != 0) stop("eval_period must be an integer")
  }
  if (delta <= 0) stop("delta must be >0")
  if (!is.null(H)){ if(H <= 0) stop("H must be >0")}
  if (zero_rate < 0 | zero_rate > 1) stop("zero_rate must be in range [0, 1]")

  # d2 is an anti-biasing constant used to estimate sigma. This d2 is set at a
  # subgroup size of 2
  d2 <- 2 / sqrt(pi)

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
    # If all conditions are met, run CUSUM test
    ctrl_period <- df$event[1:(nrow(df) - 1)]
    mu <- mean(ctrl_period)
    sigma <- mean(abs(diff(ctrl_period))) / d2
    K <- sigma * delta / 2
    H <- ifelse(is.null(H), 5 * sigma, H)
    Cplus <- Cminus <- c(0)
    for (i in 1:length(df$event)){
      Cplus[i] <- max(0, df$event[i] - (mu + K) + Cplus[i - 1])
      Cminus[i] <- max(0, (mu - K) - df$event[i] + Cminus[i - 1])
    }
    # Save output parameters
    stat <- Cplus
    # Cminus currently not used
    thresh <- H
    sig <- Cplus[length(Cplus)] > H
    hyp <- "CUSUM > decision interval"

    rr <- list(statistic=stats::setNames(stat, rep("CUSUM", length(stat))),
               lcl=rep(NA, length(stat)),
               ucl=rep(NA, length(stat)),
               p=stats::setNames(NA, "p-value"),
               signal=sig,
               signal_threshold=stats::setNames(H, "Decision Interval"),
               mu=mu,
               sigma=sigma,
               K=K)
    rs <- stats::setNames(T, "Success")
  }

  # Return test
  out <- list(test_name="CUSUM",
              analysis_of=analysis_of,
              status=rs,
              result=rr,
              params=list(test_hyp=hyp,
                          eval_period=eval_period,
                          zero_rate=zero_rate,
                          delta=delta,
                          H=H),
              data=rd)
  class(out) <- append(class(out), "mdsstat_test")
  return(out)
}
