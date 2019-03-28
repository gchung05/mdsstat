#' Exponentially Weighted Moving Average (EWMA)
#'
#' Test on device-events using the EWMA method based on the normal distribution.
#'
#' Function \code{ewma()} is an implementation of the EWMA
#' method originally proposed by S.W. Roberts based on the normal distribution.
#' EWMA is part of the family of statistical process control tests.
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
#' @param delta Required number of sigmas at which to detect a mean shift. Sigma
#' in this context refers to the estimated standard deviation of the EWMA
#' statistic.
#'
#' Default: \code{3} detects a EWMA shift of three sigmas.
#'
#' @param lambda Required "memory" parameter in the range (0, 1]. Lower values
#' assign more weight to prior measurements. A value of 1 indicates no memory
#' and only considers the current measurement.
#'
#' Default: \code{0.85} assigns higher weight to the current measurement with
#' "light memory" of prior measurements.
#'
#' @param zero_rate Required maximum proportion of \code{event}s in \code{df}
#' (constrained by \code{eval_period}) containing zeroes for this algorithm to
#' run. Because EWMA does not perform well on time series with many 0 values,
#' a value >0 is recommended.
#'
#' Default: \code{1/3} requires no more than 1/3 zeros in \code{event}s in
#' \code{df} in order to run.
#'
#' @param mu Optional value of the in-control process mean, typically measured
#' from historical data.
#'
#' Default: \code{NULL} estimates the in-control process mean from timepoints
#' prior to the most recent timepoint in the time series. The most recent
#' measurement is then tested using this estimate.
#'
#' @param sigma Optional value of the in-control process standard deviation,
#' typically measured from historical data.
#'
#' Default: \code{NULL} estimates the in-control process standard deviation
#' from timepoints prior to the most recent timepoint in the time series using
#' a moving range calculation assuming an n=2 sampling approach. The most recent
#' measurement is then tested using this estimate.
#'
#' @param ... Further arguments passed onto \code{ewma} methods
#'
#' @examples
#' # Basic Example
#' data <- data.frame(time=c(1:25), event=as.integer(stats::rnorm(25, 100, 25)))
#' a1 <- ewma(data)
#' # Example using an mds_ts object
#' a2 <- ewma(mds_ts[[3]])
#' # Example using a derived rate as the "event"
#' data <- mds_ts[[3]]
#' data$rate <- ifelse(is.na(data$nA), 0, data$nA) / data$exposure
#' a3 <- ewma(data, c(Rate="rate"))
#'
#' @references
#' S. W. Roberts (1959) Control Chart Tests Based on Geometric Moving Averages, Technometrics, 1:3, 239-250, DOI: 10.1080/00401706.1959.10489860
#' @export
ewma <- function (df, ...) {
  UseMethod("ewma", df)
}

#' @describeIn ewma EWMA on mds_ts data
#' @export
ewma.mds_ts <- function(
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
  ewma.default(out, analysis_of=name, ...)
}

#' @describeIn ewma EWMA on general data
#' @export
ewma.default <- function(
  df,
  analysis_of=NA,
  eval_period=NULL,
  delta=3,
  lambda=0.85,
  zero_rate=1/3,
  mu=NULL,
  sigma=NULL,
  ...
){
  input_param_checker(df, "data.frame")
  input_param_checker(c("time", "event"), check_names=df)
  input_param_checker(eval_period, "numeric", null_ok=T, max_length=1)
  input_param_checker(delta, "numeric", null_ok=F, max_length=1)
  input_param_checker(lambda, "numeric", null_ok=F, max_length=1)
  input_param_checker(zero_rate, "numeric", null_ok=F, max_length=1)
  input_param_checker(mu, "numeric", null_ok=T, max_length=1)
  input_param_checker(sigma, "numeric", null_ok=T, max_length=1)
  if (!is.null(eval_period)){
    if (eval_period %% 1 != 0) stop("eval_period must be an integer")
  }
  if (delta <= 0) stop("delta must be >0")
  if (lambda <= 0 | lambda > 1) stop("lambda must be in range (0, 1]")
  if (zero_rate < 0 | zero_rate > 1) stop("zero_rate must be in range [0, 1]")
  if (!is.null(mu)){ if(mu <= 0) stop("mu must be >0")}
  if (!is.null(sigma)){ if(sigma <= 0) stop("sigma must be >0")}

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
    # If all conditions are met, run EWMA test
    # Set control limits
    ctrl_period <- df$event[1:(nrow(df) - 1)]
    if (is.null(mu)) mu <- mean(ctrl_period)
    if (is.null(sigma)) sigma <- mean(abs(diff(ctrl_period))) / d2
    ewma_sigma <- sigma * sqrt(lambda / (2 - lambda) *
                                 (1 - (1 - lambda)) * (2 * nrow(df)))
    lcl <- mu - delta * ewma_sigma
    ucl <- mu + delta * ewma_sigma
    # Calculate EWMA Values
    ewma <- vector()
    for(i in 1:length(df$event)){
      if(i==1){
        ewma[1] <- lambda * df$event[1]
      } else{
        ewma[i] <- lambda * df$event[i] + (1 - lambda) * ewma[i - 1]
      }
    }

    # Save output parameters
    stat <- ewma
    thresh <- ucl
    sig <- stat[length(stat)] > thresh
    hyp <- paste0("EWMA > ", delta, "-sigma UCL")

    rr <- list(statistic=stats::setNames(stat, rep("EWMA", length(stat))),
               lcl=rep(lcl, length(stat)),
               ucl=rep(ucl, length(stat)),
               p=stats::setNames(NA, "p-value"),
               signal=sig,
               signal_threshold=stats::setNames(ucl,
                                                paste0(delta, "-sigma UCL")),
               mu=mu,
               sigma=sigma,
               ewma_sigma=ewma_sigma)
    rs <- stats::setNames(T, "Success")
  }

  # Return test
  out <- list(test_name="EWMA",
              analysis_of=analysis_of,
              status=rs,
              result=rr,
              params=list(test_hyp=hyp,
                          eval_period=eval_period,
                          zero_rate=zero_rate,
                          delta=delta,
                          lambda=lambda),
              data=rd)
  class(out) <- append(class(out), "mdsstat_test")
  return(out)
}
