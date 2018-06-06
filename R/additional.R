#' Poisson for Rare Events
#'
#' Test on rare events using an exact test on the Poisson distribution rate
#' parameter (\code{stats::poisson.test()}).
#'
#' \code{p_rate} default of \code{0.2} is a suggested null value for
#' the Poisson rate parameter. However this value is highly advised to be set
#' based on known priors and/or your specific application.
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
#' times counting in reverse chronological order to assess.
#'
#' Default: \code{NULL} considers all times in \code{df}.
#'
#' @param zero_rate Required minimum proportion of \code{event}s in \code{df}
#' (constrained by \code{eval_period}) containing zeroes for this algorithm to
#' run.
#'
#' Default: \code{2/3} requires a minimum of 2/3 zeros in \code{event}s in
#' \code{df}.
#'
#' @param p_rate Hypothesized Poisson rate parameter null value at which the
#' Poisson test is performed (null vs. greater). See details for more.
#'
#' Default: \code{0.2}
#'
#' @param p_crit Critical p-value for the Poisson test..
#'
#' Default: \code{0.05}
#'
#' @param ... Further arguments passed onto \code{poisson_rare} methods
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
#' @examples
#' # Basic Example
#' data <- data.frame(time=c(1:8), event=c(rep(0, 6), stats::rpois(2, 4)))
#' a1 <- poisson_rare(data)
#' # Example using an mds_ts object
#' a2 <- poisson_rare(mds_ts[[1]])
#' # Example using a derived rate as the "event"
#' data <- mds_ts[[1]]
#' data$rate <- ifelse(is.na(data$nA), 0, data$nA) / data$exposure
#' a3 <- poisson_rare(data, c("Rate"="rate"))
#' @export
poisson_rare <- function (df, ...) {
  UseMethod("poisson_rare", df)
}

#' @describeIn poisson_rare Poisson on mds_ts data
#' @export
poisson_rare.mds_ts <- function(
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
                  paste(attributes(df)$nLabels$nA, collapse=" for "))
  } else name <- analysis_of

  out <- data.frame(time=df$time,
                    event=df[[ts_event]])
  poisson_rare.default(out, analysis_of=name, ...)
}

#' @describeIn poisson_rare Poisson on general data
#' @export
poisson_rare.default <- function(
  df,
  analysis_of=NA,
  eval_period=NULL,
  zero_rate=2/3,
  p_rate=0.2,
  p_crit=0.05,
  ...
){
  input_param_checker(df, "data.frame")
  input_param_checker(c("time", "event"), check_names=df)
  input_param_checker(eval_period, "integer")
  input_param_checker(zero_rate, "numeric", null_ok=F, max_length=1)
  input_param_checker(p_rate, "numeric", null_ok=F, max_length=1)
  input_param_checker(p_crit, "numeric", null_ok=F, max_length=1)
  if (zero_rate < 0 | zero_rate > 1) stop("zero_rate must be in range [0, 1]")
  if (p_crit < 0 | p_crit > 1) stop("p_crit must be in range [0, 1]")

  # Order by time
  df <- df[order(df$time), ]
  # Restrict to eval_period
  if (!is.null(eval_period)){
    if (eval_period > nrow(df)){
      stop("eval_period cannot be greater than df rows")
    } else if (eval_period < 1){
      stop("eval_period must be greater than 0")
    } else df <- df[c((nrow(df) - eval_period + 1):nrow(df)), ]
  }
  # Return data
  rd <- list(reference_time=range(df$time),
             data=df)
  # Set Poisson test type
  h_alternative <- "greater"
  # Check for non-runnable conditions
  if(nrow(df) < 4){
    rr <- NA
    rs <- stats::setNames(F, ">3 time periods required")
  } else if(sum(df$time != 0) < 2){
    rr <- NA
    rs <- stats::setNames(F, "2 or more non-zero events required")
  } else if(sum(df$event == 0) / nrow(df) < zero_rate){
    rr <- NA
    rs <- stats::setNames(F, paste("Minimum zero_rate of", zero_rate, "not met"))
  } else{
    # If all conditions are met, run Poisson test
    run <- stats::poisson.test(round(sum(df$event)), nrow(df), r=p_rate,
                              alternative=h_alternative)
    rr <- list(statistic=stats::setNames(run$estimate, "Rate"),
               lcl=run$conf.int[1],
               ucl=run$conf.int[2],
               p=stats::setNames(run$p.value, "p-value"),
               signal=(run$p.value <= p_crit),
               signal_threshold=stats::setNames(p_crit, "critical p-value"))
    rs <- stats::setNames(T, "Success")
  }

  # Return test
  out <- list(test_name="Poisson Rare",
              analysis_of=analysis_of,
              status=rs,
              result=rr,
              params=list(test_hyp=paste0("Poisson test p-value <=", p_crit),
                          eval_period=eval_period,
                          zero_rate=zero_rate,
                          p_rate=p_rate,
                          p_crit=p_crit,
                          h_alternative=h_alternative),
              data=rd)
  class(out) <- append(class(out), "mdsstat_test")
  return(out)
}
