#' Shewhart x-bar Control Chart
#'
#' Test on device-events using the Shewhart x-bar control chart. Includes
#' the first 4 Western Electric rules common to statistical process control.
#'
#' Function \code{shewhart()} is a standard implementation of the x-bar
#' Control Chart test from the family of statistical process control tests
#' originally proposed by Walter Shewhart.
#'
#' \code{we_rule} has four possible values: \code{1} is one point over the
#' 3-sigma limit. \code{2} is two out of three consecutive points over the
#' 2-sigma limit. \code{3} is four of five consecutive points over the 1-sigma
#' limit. \code{4} is nine consecutive points over the process mean.
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
#' @param zero_rate Required maximum proportion of \code{event}s in \code{df}
#' (constrained by \code{eval_period}) containing zeroes for this algorithm to
#' run. Because Shewhart does not perform well on time series with many 0 values,
#' a value >0 is recommended.
#'
#' Default: \code{1/3} requires no more than 1/3 zeros in \code{event}s in
#' \code{df} in order to run.
#'
#' @param we_rule Required integer from \code{1} to \code{4} representing the
#' Western Electric rule to use. See details for descriptions.
#'
#' Default: \code{1} represents the first Western Electric rule of one point
#' over the 3-sigma limit.
#'
#' @param ... Further arguments passed onto \code{shewhart} methods
#'
#' @examples
#' # Basic Example
#' data <- data.frame(time=c(1:25), event=as.integer(stats::rnorm(25, 100, 25)))
#' a1 <- shewhart(data)
#' # Example using an mds_ts object
#' a2 <- shewhart(mds_ts[[3]])
#' # Example using a derived rate as the "event"
#' data <- mds_ts[[3]]
#' data$rate <- ifelse(is.na(data$nA), 0, data$nA) / data$exposure
#' a3 <- shewhart(data, c(Rate="rate"))
#'
#' @references
#' Montgomery, Douglas C. Introduction to Statistical Quality Control by Douglas C. Montgomery, 5th Edition: Study Guide. Cram101, 2013.
#' @export
shewhart <- function (df, ...) {
  UseMethod("shewhart", df)
}

#' @describeIn shewhart Shewhart on mds_ts data
#' @export
shewhart.mds_ts <- function(
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
  shewhart.default(out, analysis_of=name, ...)
}

#' @describeIn shewhart Shewhart on general data
#' @export
shewhart.default <- function(
  df,
  analysis_of=NA,
  eval_period=NULL,
  zero_rate=1/3,
  we_rule=1L,
  ...
){
  input_param_checker(df, "data.frame")
  input_param_checker(c("time", "event"), check_names=df)
  input_param_checker(zero_rate, "numeric", null_ok=F, max_length=1)
  input_param_checker(we_rule, "integer", null_ok=F, max_length=1)
  input_param_checker(eval_period, "numeric", null_ok=F, max_length=1)
  if (eval_period %% 1 != 0) stop("eval_period must be an integer")
  if (zero_rate < 0 | zero_rate > 1) stop("zero_rate must be in range [0, 1]")
  if (we_rule < 1 | we_rule > 4) stop("we_rule must be in range [1L, 4L]")

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
  } else if (we_rule == 2L & (nrow(df) < 6)){
    rr <- NA
    rs <- stats::setNames(F, ">5 time periods required for WE rule 2")
  } else if (we_rule == 3L & (nrow(df) < 8)){
    rr <- NA
    rs <- stats::setNames(F, ">7 time periods required for WE rule 3")
  } else if (we_rule == 4L & (nrow(df) < 11)){
    rr <- NA
    rs <- stats::setNames(F, ">10 time periods required for WE rule 4")
  } else if(sum(df$time != 0) < 2){
    rr <- NA
    rs <- stats::setNames(F, "2 or more non-zero events required")
  } else if(sum(df$event == 0) / nrow(df) > zero_rate){
    rr <- NA
    rs <- stats::setNames(F, paste("Maximum zero_rate of", zero_rate, "exceeded"))
  } else{
    # If all conditions are met, run Shewhart test
    ctrl_period <- df$event[1:(nrow(df) - 1)]
    mu <- mean(ctrl_period)
    sigma <- mean(abs(diff(ctrl_period))) / d2
    nsigma <- (df$event - mu) / sigma
    # Calculate WE trend rules
    if (we_rule == 1L){
      stat <- nsigma[tlen]
      thresh <- (mu + 3 * sigma)
      sig <- stat > 3
      hyp <- "1 point > 3-sigma limit"
    } else if (we_rule == 2L){
      stat <- nsigma[(tlen - 2):tlen]
      thresh <- mu + 2 * sigma
      sig <- (nsigma[tlen] > 2) &
        any(nsigma[(tlen - 2):(tlen - 1)] > 2)
      hyp <- "2 of 3 points > 2-sigma limit"
    } else if (we_rule == 3L){
      stat <- nsigma[(tlen - 4):tlen]
      thresh <- mu + sigma
      sig <- (nsigma[tlen] > 1) &
        (sum(nsigma[(tlen - 4):(tlen - 1)] > 1) >= 3)
      hyp <- "4 of 5 points > 1-sigma limit"
    } else if (we_rule == 4L){
      stat <- nsigma[(tlen - 8):tlen]
      thresh <- mu
      sig <- sum(stat > 0) == 9
      hyp <- "9 points > mean"
    }

    rr <- list(statistic=stats::setNames(stat, rep("N-Sigmas", length(stat))),
               lcl=rep(mu + stats::qnorm(0.025) * sigma, length(stat)),
               ucl=rep(mu + stats::qnorm(0.975) * sigma, length(stat)),
               p=stats::setNames(NA, "p-value"),
               signal=sig,
               signal_threshold=stats::setNames(
                 rep(thresh, length(stat)),
                 rep("UCL", length(stat))),
               mu=mu,
               sigma=sigma)
    rs <- stats::setNames(T, "Success")
  }

  # Return test
  out <- list(test_name=paste("Shewhart x-bar Western Electric Rule", we_rule),
              analysis_of=analysis_of,
              status=rs,
              result=rr,
              params=list(test_hyp=hyp,
                          eval_period=eval_period,
                          zero_rate=zero_rate,
                          we_rule=we_rule),
              data=rd)
  class(out) <- append(class(out), "mdsstat_test")
  return(out)
}
