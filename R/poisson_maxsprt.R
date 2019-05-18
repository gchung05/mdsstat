#' Poisson MaxSPRT
#'
#' Test on device-events using the Poisson Maximized Sequential Probability
#' Ratio Test, or Poisson MaxSPRT, originally proposed by Kulldorff et al
#' (2011). Requires a call to Silva & Kulldorff's \code{Sequential} package.
#'
#' Runs a variant of Kulldorff's Poisson MaxSPRT where the mean of the null
#' Poisson distribution, \code{u_t}, reflecting the baseline risk of the event,
#' can be inferred from the data as the evaluation period, \code{eval_period},
#' less the observation period, \code{obs_period}.
#'
#' This null distribution is used in conjunction with \code{alpha}, \code{M},
#' and other parameters passed into \code{Sequential::SampleSize.Poisson()} to
#' establish the critical MaxSPRT log likelihood ratios. Additional parameters
#' eligible to be passed via \code{...} are \code{power}, \code{D}, and
#' \code{RR}.
#'
#' The following constraints are on this test: 1) If \code{u_t} is being
#' inferred from the data, the total number of events used to construct
#' \code{u_t} must be at least
#' equal to the required sample size in order to control Type I error, per
#' Li 2009.
#' 2) The total number of observed events cannot be greater than the length
#' of surveillance (required sample size).
#'
#' For parameter \code{ts_event}, in the uncommon case where the
#' device-event count (Cell A) variable is not \code{"nA"}, the name of the
#' variable may be specified here. A named character
#' vector may be used where the name is the English description of what was
#' analyzed. Note that if the parameter \code{analysis_of} is specified, it will
#' override this name. Example: \code{ts_event=c("Count of Bone Cement
#' Leakages"="event_count")}
#'
#' @param df Required input data frame of class \code{mds_ts} or, for generic
#' usage, any data frame with the following columns:
#' \describe{
#'   \item{time}{Unique times of class \code{Date}}
#'   \item{event}{Either the event count or rate of class \code{numeric}}
#' }
#' @param ts_event Required if \code{df} is of class \code{mds_ts}. Named string
#' indicating the variable corresponding to the event count (cell A in the 2x2
#' contingency table). In most cases, the default is the appropriate setting.
#' See details for alternative options.
#'
#' Default: \code{c("Count"="nA")} corresponding to the event count column in
#' \code{mds_ts} objects. Name is generated from \code{mds_ts} metadata.
#'
#' @param analysis_of Optional string indicating the English description of what
#' was analyzed. If specified, this will override the name of the
#' \code{ts_event} string parameter.
#'
#' Default: \code{NA} indicates no English description for plain \code{df}
#' data frames, or \code{ts_event} English description for \code{df} data frames
#' of class \code{mds_ts}.
#'
#' Example: \code{"Count of bone cement leakages"}
#'
#' @param eval_period Positive integer indicating the number of unique
#' times counting in reverse chronological order to evaluate. If \code{u_t} is
#' not specified, \code{eval_period} less \code{obs_period} will be used
#' to establish an empirical null Poisson distribution, \code{u_t}.
#'
#' Default: \code{NULL} considers all times in \code{df}.
#'
#' @param obs_period Positive integer indicating the number of unique
#' times in reverse chronological order to observe and test against the null
#' hypothesis, described by \code{u_t}. This cannot be greater than
#' \code{eval_period}.
#'
#' Default: \code{3} indicates the last 3 times from \code{df} constitute the
#' observation period. It is not \code{1} to emphasize its use on rarely
#' occurring events. The user is encouraged to explicitly set \code{obs_period}.
#'
#' @param alpha Required Type I error probability in the range (0, 0.5].
#'
#' Default: \code{0.05} is a significance level of 5\%.
#'
#' @param M Required minimum number of events needed during the observation
#' period \code{obs_period} before the null hypothesis can be
#' rejected. Passed into \code{Sequential::SampleSize.Poisson()}.
#'
#' Default: \code{4} is a minimum of 4 events, per Kulldorff and Silva, 2015.
#'
#' @param u_t Optional value for the mean of the null Poisson distribution used
#' to test the hypothesis.
#'
#' Default: \code{NA} will infer the null distribution from the data using the
#' \code{eval_period} less the \code{obs_period}.
#'
#' @param ... Further arguments passed onto \code{poisson_maxsprt} methods,
#' which are currently select parameters for
#' \code{Sequential::SampleSize.Poisson()} as follows: \code{RR}, \code{power},
#' and \code{D}.
#'
#' @return A named list of class \code{mdsstat_test} object, as follows:
#' \describe{
#'   \item{test_name}{Name of the test run}
#'   \item{analysis_of}{English description of what was analyzed}
#'   \item{status}{Named boolean of whether the test was run. The name contains
#'   the run status.}
#'   \item{result}{A standardized list of test run results: \code{statistic}
#'   for the test statistic, \code{lcl} and \code{ucl} for the set
#'   confidence bounds, \code{p} for the p-value, \code{signal} status,
#'   \code{signal_threshold}, and other MaxSPRT test values.}
#'   \item{params}{The test parameters}
#'   \item{data}{The data on which the test was run}
#' }
#'
#' @examples
#' # Basic Example
#' data <- data.frame(time=c(1:25), event=as.integer(stats::rnorm(25, 100, 25)))
#' a1 <- poisson_maxsprt(data)
#' # Example using an mds_ts object
#' a2 <- poisson_maxsprt(mds_ts[[3]])
#'
#' @references
#' Martin Kulldorff, Robert L. Davis, Margarette Kolczak, Edwin Lewis, Tracy Lieu & Richard Platt (2011) A Maximized Sequential Probability Ratio Test for Drug and Vaccine Safety Surveillance, Sequential Analysis, 30:1, 58-78.
#'
#' Ivair Ramos Silva, Martin Kulldorff (2019). Sequential: Exact Sequential Analysis for Poisson and Binomial Data. R package version 3.0.1. https://CRAN.R-project.org/package=Sequential
#'
#' Lingling Li, Martin Kulldorff (2009). A conditional maximized sequential probability ratio test for Pharmacovigilance. Statistics in Medicine. doi:10.1002/sim.3780
#'
#' @export
poisson_maxsprt <- function (df, ...) {
  UseMethod("poisson_maxsprt", df)
}

#' @describeIn poisson_maxsprt Poisson MaxSPRT on mds_ts data
#' @export
poisson_maxsprt.mds_ts <- function(
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
  poisson_maxsprt.default(out, analysis_of=name, ...)
}

#' @describeIn poisson_maxsprt Poisson MaxSPRT on general data
#' @export
poisson_maxsprt.default <- function(
  df,
  analysis_of=NA,
  eval_period=NULL,
  obs_period=3,
  alpha=0.05,
  M=4,
  u_t=NA,
  ...
){
  input_param_checker(df, "data.frame")
  input_param_checker(c("time", "event"), check_names=df)
  input_param_checker(eval_period, "numeric", null_ok=T, max_length=1)
  input_param_checker(obs_period, "numeric", null_ok=F, max_length=1)
  input_param_checker(alpha, "numeric", null_ok=F, max_length=1)
  input_param_checker(M, "numeric", null_ok=F, max_length=1)
  input_param_checker(u_t, "numeric", null_ok=F, na_ok=T, max_length=1)
  if (!is.null(eval_period)){
    if (eval_period < 2) stop("eval_period must be 2 or greater")
    if (eval_period %% 1 != 0) stop("eval_period must be an integer")
    if (obs_period > eval_period){
      stop("obs_period cannot be greater than eval_period")
    }
    if (is.na(u_t) & obs_period == eval_period){
      stop("obs_period cannot equal eval_period unless u_t is declared")
    }
  }
  if (is.na(u_t) & obs_period == nrow(df)){
    stop("obs_period cannot equal df rows unless u_t is declared")
  }
  if (obs_period %% 1 != 0 | obs_period <= 0){
    stop("obs_period must be a positive integer")
  }
  if (obs_period > nrow(df)){
    stop("obs_period cannot be longer than rows in df (number of times)")}
  if (alpha <= 0 | alpha > 0.5) stop("alpha must be in range (0, 0.5]")
  if (M %% 1 != 0) stop("M must be an integer")
  if (M < 1) stop("D must be 1 or greater")
  if (!is.na(u_t)){ if (u_t <= 0) stop("u_t must be >0")}
  # Ellipsis parameters for Sequential::SampleSize.Poisson()
  dots <- list(...)
  RR <- ifelse(is.null(dots$RR),
               formals(Sequential::SampleSize.Poisson)$RR, dots$RR)
  power <- ifelse(is.null(dots$power),
                  formals(Sequential::SampleSize.Poisson)$power, dots$power)
  D <- ifelse(is.null(dots$D),
              formals(Sequential::SampleSize.Poisson)$D, dots$D)

  u_t_input <- u_t

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
  rr <- NA
  if(nrow(df) < 2){
    rs <- stats::setNames(F, ">=2 time periods required")
  } else if(is.null(u_t) & (obs_period == eval_period)){
    rs <- stats::setNames(F, paste("null distribution (u_t) cannot be inferred",
                                   "from data if all data is in observation",
                                   "period"))
  } else{
    # If all conditions are met, run Poisson MaxSPRT test
    n_obs <- sum(df$event[c((nrow(df) - obs_period + 1):nrow(df))])
    n_eval <- sum(df$event[c((nrow(df) - eval_period + 1)):
                             c((nrow(df) - obs_period))])
    # Estimation of sample size and critical value
    sspoi <- Sequential::SampleSize.Poisson(alpha=alpha, M=M,
                                            RR=RR, power=power, D=D,
                                            precision=0.01)
    sspoi_ss <- sspoi[, "Sample Size"]
    sspoi_cv <- sspoi[, "Critical value"]
    if (is.na(u_t) & (n_eval < sspoi_ss)){
      # If the data are being used to construct the null, only run test if
      # total events under the null are at least than "length of surveillance",
      # as defined in Kulldorff et al (2011), in order to control Type I error.
      rs <- stats::setNames(F, paste0(
        "since u_t is being inferred from the data, ",
        "total events during eval_period less obs_period (n_eval=", n_eval,
        ") must be at least equal to the ",
        "length of surveillance (sspoi_ss=", round(sspoi_ss, 1),
        " events), which is a function of parameters: ",
        "alpha, M, D, power, RR"))
    } else if (n_obs > sspoi_ss){
      rs <- stats::setNames(F, paste0(
        "total events during observation period (n_obs=", n_eval,
        ") exceeds the length of surveillance (sspoi_ss=", round(sspoi_ss, 1)))
    } else if (n_obs < M){
      # Only run test if minimum number of events needed, M, is met
      rs <- stats::setNames(F, paste0("minimum number of events M=", M,
                                      " was not exceeded during observation ",
                                      "period (n=", n_obs, ")"))
    } else if ((is.na(u_t) & n_eval == 0) | (!is.na(u_t) & u_t == 0)){
      # Only run test if the null Poisson estimate is > 0
      rs <- stats::setNames(F, paste("null distribution (u_t) cannot be",
                                     "inferred from data if no events have",
                                     "occurred"))
    } else{
      # Calculate MaxSPRT
      c_t <- n_obs / obs_period
      # Estimate u_t from evaluation period excluding observation period
      if (is.na(u_t)) u_t <- n_eval / (eval_period - obs_period)
      if (c_t < u_t){
        llr <- 0
        sig <- F
      } else{
        llr <- (u_t - c_t) + c_t * log(c_t / u_t)
        sig <- as.logical(llr > sspoi_cv)
      }
      hyp <- paste("Observed", round(c_t, 1),
                   "events greater than expected null", round(u_t, 1),
                   "events based on a length of surveillance of",
                   round(sspoi_ss, 1),
                   "events")

      rr <- list(statistic=stats::setNames(llr, "Log-Likelihood Ratio"),
                 lcl=NA,
                 ucl=NA,
                 p=stats::setNames(NA, "p-value"),
                 signal=sig,
                 signal_threshold=sspoi_cv,
                 sample_threshold=sspoi_ss,
                 c_t=c_t,
                 u_t=u_t)
      rs <- stats::setNames(T, "Success")
    }

    # Return test
    out <- list(test_name="Poisson MaxSPRT",
                analysis_of=analysis_of,
                status=rs,
                result=rr,
                params=list(test_hyp=hyp,
                            eval_period=eval_period,
                            obs_period=obs_period,
                            alpha=alpha,
                            M=M,
                            u_t=u_t_input),
                data=rd)
    class(out) <- append(class(out), "mdsstat_test")
    return(out)
  }
}
