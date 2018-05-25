#' Proportional Reporting Ratio
#'
#' Test on device-events using the proportional reporting ratio (PRR). From
#' the family of disproportionality analyses (DPA) used to generate signals of
#' disproportionate reporting (SDRs).
#'
#' @param df Required input data frame of class \code{mds_ts} or, for generic
#' usage, any data frame with the following columns:
#' \describe{
#'   \item{time}{Unique times of class \code{Date}}
#'   \item{nA}{Cell A count (class \code{numeric}) of the 2x2 table:
#'   device/event of interest.}
#'   \item{nB}{Cell B count (class \code{numeric}) of the 2x2 table:
#'   device/non-event of interest.}
#'   \item{nC}{Cell C count (class \code{numeric}) of the 2x2 table:
#'   non-device/event of interest.}
#'   \item{nD}{Cell D count (class \code{numeric}) of the 2x2 table:
#'   non-device/non-event of interest.}
#' }
#' @param ts_event Required if \code{df} is of class \code{mds_ts}. Named string
#' indicating the variable corresponding to the event count (cell A in the 2x2
#' contingency table). In most cases, the default (\code{"nA"}) is the appropriate
#' setting. Otherwise specify the name of the alternate variable containing
#' event counts for cell A. The name of the string is an English description of
#' what was analyzed.
#'
#' Default: \code{"nA"} corresponding to the event count column in \code{mds_ts}
#' objects. Name is generated from \code{mds_ts} metadata.
#'
#' Example: \code{stats::setNames("eventCount", "Count of Bone Cement Leakages")}
#'
#' @param analysis_of Optional string indicating the English description of what
#' was analyzed. If specified, this will override the name of the
#' \code{ts_event} string parameter.
#'
#' Default: \code{NA} indicates no English description.
#'
#' Example: \code{"Count of bone cement leakages"}
#'
#' @param eval_period Required positive integer indicating the number of unique
#' times counting in reverse chronological order to sum over to create the 2x2
#' contingency table.
#'
#' Default: \code{1} considers only the most recent time in \code{df}.
#'
#' EXample: \code{12} sums over the last 12 time periods to create the 2x2
#' contingency table.
#'
#' @return A named list of class \code{mdsstat_test} object, as follows:
#' \describe{
#'   \item{test_name}{Name of the test run}
#'   \item{analysis_of}{English description of what was analyzed}
#'   \item{status}{Named boolean of whether the test was run. The name contains
#'   the run status.}
#'   \item{result}{A standardized list of test run results: \code{statistic}
#'   for the test statistic, \code{ll95} and \code{ul95} for the 95%
#'   confidence bounds, \code{p} for the p-value, \code{signal} status, and
#'   \code{signal_threshold}.}
#'   \item{params}{The test parameters}
#'   \item{data}{The data on which the test was run}
#' }
#'
#' @examples
#' # Basic Example
#' data <- data.frame(time=c(1:8), event=c(rep(0, 6), rpois(2, 4)))
#' a1 <- poisson_rare(data)
#' # Example using an mds_ts object
#' a2 <- poisson_rare(mds_ts[[1]])
#' # Example using a derived rate as the "event"
#' data <- mds_ts[[1]]
#' data$rate <- ifelse(is.na(data$nA), 0, data$nA) / data$exposure
#' a3 <- poisson_rare(data, "rate")
#'
#' @references
#' Evans, S. J. W., Waller, P. C., & Davis, S. (2001). Use of proportional reporting ratios (PRRs) for signal generation from spontaneous adverse drug reaction reports. Pharmacoepidemiology and Drug Safety, 10(6), 483–486. https://doi.org/10.1002/pds.677
#' @export
prr <- function (df, ...) {
  UseMethod("prr", df)
}

prr.mds_ts <- function(
  df,
  ts_event="nA",
  analysis_of=NA,
  ...
){
  input_param_checker(ts_event, check_names=df, max_length=1)

  # Set NA counts to 0 for "nA" default
  df$nA <- ifelse(is.na(df$nA), 0, df$nA)

  # Set analysis_of
  if (ts_event == "nA" & is.na(analysis_of)){
    name <- paste("Count of", paste(sapply(attributes(df)$nA, function(x){
      paste0(names(x), ":", x)}), collapse=" for "))
  } else if (is.na(analysis_of)){
    name <- names(ts_event)
  } else name <- analysis_of

  out <- data.frame(time=df$time,
                    event=df[[ts_event]])
  shewhart.default(out, analysis_of=name, ...)
}

prr.default <- function(
  df,
  analysis_of=NA,
  eval_period=1
){
  input_param_checker(df, "data.frame")
  input_param_checker(eval_period, "integer")
  input_param_checker(zero_rate, "numeric", null_ok=F, max_length=1)
  input_param_checker(we_rule, "integer", null_ok=F, max_length=1)
  if (!all(c("time", "event") %in% names(df))){
    stop("df must contain columns named time and event")
  }
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
  }
  # Return data
  tlen <- nrow(df)
  rd <- list(reference_time=df$time[tlen],
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
    # If all conditions are met, run Shewhart test
    mu <- mean(df$event)
    sigma <- mean(abs(diff(df$event))) / d2
    df$event[length(df$event)]
    nsigma <- (df$event - mu) / sigma
    if (we_rule == 1L){
      stat <- nsigma[tlen]
      thresh <- (mu + 3 * sigma)
      sig <- stat > thresh
      hyp <- "1 point > 3-sigma limit"
    } else if (we_rule == 2L){
      stat <- nsigma[(tlen - 2):tlen]
      thresh <- mu + 2 * sigma
      sig <- (nsigma[tlen] > thresh) &
        any(nsigma[(tlen - 2):(tlen - 1)] > thresh)
      hyp <- "2 of 3 points > 2-sigma limit"
    } else if (we_rule == 3L){
      stat <- nsigma[(tlen - 4):tlen]
      thresh <- mu + sigma
      sig <- (nsigma[tlen] > thresh) &
        (sum(nsigma[(tlen - 4):(tlen - 1)] > thresh) >= 3)
      hyp <- "4 of 5 points > 1-sigma limit"
    } else if (we_rule == 4L){
      stat <- nsigma[(tlen - 8):tlen]
      thresh <- mu
      sig <- sum(stat > thresh) == 9
      hyp <- "9 points > mean"
    }

    rr <- list(statistic=stat,
               ll95=rep(mu + qnorm(0.025) * sigma, length(stat)),
               ul95=rep(mu + qnorm(0.975) * sigma, length(stat)),
               p=stats::setNames(NA, "p-value"),
               signal=sig,
               signal_threshold=stats::setNames(
                 rep(thresh, length(stat)),
                 rep("UCL", length(stat))),
               mu=mu,
               sigma=sigma)
    rm <- "Run success"
    rs <- stats::setNames(T, "Success")
  }

  # Return test
  out <- list(test_name=paste("Shewhart x-bar Western Electric Rule", we_rule),
              analysis_of=analysis_of,
              status=rs,
              result=rr,
              params=list(test_hyp=hyp,
                          zero_rate=zero_rate,
                          we_rule=we_rule),
              data=rd)
  class(out) <- append(class(out), "mdsstat_test")
  return(out)
}
