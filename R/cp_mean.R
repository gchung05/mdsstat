#' Mean-Shift Changepoint
#'
#' Test on device-events using the mean-shift changepoint method
#' originally described in Xu, et al 2015.
#'
#' Function \code{cp_mean()} is an implementation of the mean-shift changepoint
#' method originally proposed by Xu, et al (2015) based on a bootstrap null
#' distribution. The parameters in this implementation can be interpreted as
#' follows. Changepoints are detected at an \code{alpha} level based on
#' n=\code{bootstrap_iter} bootstrap iterations (with or without replacement
#' using \code{replace}) of the input time series
#' \code{df}. A minimum of n=\code{min_seglen} consecutive measurements without
#' a changepoint are required to test for an additional changepoint. Both
#' \code{epochs} and \code{cpmax} constrain the maximum possible number of
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
#' the range (0, 1].
#'
#' Default: \code{0.05} detects a changepoint at an alpha level of 0.05 or 5%.
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
#' Xu, Zhiheng, et al. "Signal detection using change point analysis in postmarket surveillance." pharmacoepidemiology and drug safety 24.6 (2015): 663-668.
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
  cpmax=100,
  min_seglen=6,
  epochs=NULL,
  bootstrap_iter=1000,
  replace=T,
  ...
){
  input_param_checker(df, "data.frame")
  input_param_checker(c("time", "event"), check_names=df)
  input_param_checker(eval_period, "numeric", null_ok=T, max_length=1)
  input_param_checker(alpha, "numeric", null_ok=F, max_length=1)
  input_param_checker(cpmax, "numeric", null_ok=F, max_length=1)
  input_param_checker(min_seglen, "numeric", null_ok=F, max_length=1)
  input_param_checker(epochs, "numeric", null_ok=T, max_length=1)
  input_param_checker(bootstrap_iter, "numeric", null_ok=F, max_length=1)
  input_param_checker(replace, "logical", null_ok=F, max_length=1)
  if (!is.null(eval_period)){
    if (eval_period %% 1 != 0) stop("eval_period must be an integer")
  }
  if (alpha <= 0 | alpha >= 1) stop("alpha must be in range (0, 1)")
  if (cpmax %% 1 != 0) stop("cpmax must be an integer")
  if (cpmax <= 0 ) stop("cpmax must be positive")
  if (min_seglen %% 1 != 0) stop("min_seglen must be an integer")
  if (min_seglen < 4) stop("min_seglen must be >=4")
  if (epochs %% 1 != 0) stop("epochs must be an integer")
  if (epochs <= 0 ) stop("epochs must be positive")
  if (bootstrap_iter %% 1 != 0) stop("bootstrap_iter must be an integer")
  if (bootstrap_iter < 1000) stop("bootstrap_iter should be >=1000")

  # I AM HERE!!! READY TO WRITE IN THE FUNC

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
    # If all conditions are met, run CP test
    # Set control limits
    ctrl_period <- df$event[1:(nrow(df) - 1)]
    if (is.null(mu)) mu <- mean(ctrl_period)
    if (is.null(sigma)) sigma <- mean(abs(diff(ctrl_period))) / d2
    ewma_sigma <- sigma * sqrt(lambda / (2 - lambda) *
                                 (1 - (1 - lambda)) * (2 * nrow(df)))
    lcl <- mu - delta * ewma_sigma
    ucl <- mu + delta * ewma_sigma
    # Calculate CP Values
    cp_mean <- vector()
    for(i in 1:length(df$event)){
      if(i==1){
        cp_mean[1] <- lambda * df$event[1]
      } else{
        cp_mean[i] <- lambda * df$event[i] + (1 - lambda) * cp_mean[i - 1]
      }
    }

    # Save output parameters
    stat <- cp_mean
    thresh <- ucl
    sig <- stat[length(stat)] > thresh
    hyp <- paste0("CP > ", delta, "-sigma UCL")

    rr <- list(statistic=stats::setNames(stat, rep("CP", length(stat))),
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
  out <- list(test_name="Mean-shift changepoint",
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
