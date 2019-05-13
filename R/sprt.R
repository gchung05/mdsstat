#' Sequential Probability Ratio Test
#'
#' Test on device-events using Wald's Sequential Probability Ratio Test (SPRT).
#' Supports both Poisson (default) and normal distribution functions.
#'
#' Runs Wald's SPRT where the null hypothesis \code{h0} and alternative
#' hypothesis \code{h1} represent event occurrence in a single time period.
#' Event occurrence in Wald's context is the number of events in a time period.
#' However, at the user's discretion, this function allows event occurrence to
#' be any continuous number (such as event rate).
#'
#' In typical medical device surveillance, \code{h1} is greater than \code{h0}
#' or \code{relative_risk} is greater than 1, and \code{h1_type="greater"}. This is
#' because we wish to detect elevated occurrences of an undesirable event.
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
#' @param eval_period Optional positive integer indicating the number of unique
#' times counting in reverse chronological order to evaluate. This will be used
#' to establish the default null hypothesis \code{h0}.
#'
#' Default: \code{NULL} considers all times in \code{df}.
#'
#' @param obs_period Required positive integer indicating the number of unique times
#' in reverse chronological order to observe and test against the null
#' hypothesis \code{h0}. This cannot be greater than \code{eval_period}. Used
#' with \code{eval_period} to establish the default null hypothesis \code{h0}.
#'
#' Default: \code{1} indicates only the latest time in \code{df} constitutes the
#' observation period.
#'
#' Example: \code{3} indicates the last three times in \code{df} constitute the
#' observation period.
#'
#' @param h0 Optional numeric value representing the null hypothesis. See
#' details for more.
#'
#' Default: \code{NULL} estimates the null hypothesis from the evaluation period
#' \code{eval_period} less the number of observation periods \code{obs_period}.
#'
#' @param h1 Optional numeric value representing the test/alternative
#' hypothesis. Either \code{h1} or \code{relative_risk} must be specified. If
#' both are specified, \code{relative_risk} takes priority. See details for more.
#'
#' Default: \code{NULL} assumes that \code{relative_risk} is being used to infer
#' the alternative hypothesis.
#'
#' @param relative_risk Optional numeric value representing the relative risk used to
#' infer the test/alternative hypothesis as follows: \code{h1=relative_risk * h0}.
#' Either \code{h1} or \code{relative_risk} must be specified. If both are
#' specified, \code{relative_risk} takes priority. See details for more.
#'
#' Default: \code{NULL} assumes that \code{h1} is defining the alternative
#' hypothesis.
#'
#' Example: \code{1.2} tests if event occurrence is 1.2 times what is specified
#' by the null hypothesis \code{h0}.
#'
#' @param distribution Required distribution to estimate. Must be either
#' \code{"poisson"} or \code{"normal"}.
#'
#' Default: \code{"poisson"}
#'
#' @param alpha Required Type I error probability between 0 and 1.
#'
#' Default: \code{0.05}
#'
#' @param beta Required Type II error probability between 0 and 1.
#'
#' Default: \code{0.20}
#'
#' @param h1_type Required type of alternative hypothesis. Must be any of three
#' values: \code{"greater"}, \code{"less"}, or \code{"not equal"}. For a
#' two-sided test, set to \code{"not equal"}.
#'
#' Default: \code{"greater"}
#'
#' @param ... Further arguments passed onto \code{sprt} methods
#'
#' @return A named list of class \code{mdsstat_test} object, as follows:
#' \describe{
#'   \item{test_name}{Name of the test run}
#'   \item{analysis_of}{English description of what was analyzed}
#'   \item{status}{Named boolean of whether the test was run. The name contains
#'   the run status.}
#'   \item{result}{A standardized list of test run results: \code{statistic}
#'   for the test statistic, \code{lcl} and \code{ucl} for the set
#'   confidence bounds, \code{p} for the p-value, \code{signal} status, and
#'   \code{signal_threshold}.}
#'   \item{params}{The test parameters}
#'   \item{data}{The data on which the test was run}
#' }
#'
#' @examples
#' # At minimum, the df parameter and either the h1 or relative_risk parameter
#' # must be specified.
#' # Basic Example
#' data <- data.frame(time=c(1:25), event=as.integer(stats::rnorm(25, 100, 25)))
#' a1 <- sprt(data, h1=110)
#' # Example using an mds_ts object
#' a2 <- sprt(mds_ts[[3]], relative_risk=1.2)
#'
#' @references
#' Wald, Abraham (June 1945). "Sequential Tests of Statistical Hypotheses". Annals of Mathematical Statistics. 16 (2): 117-186.
#'
#' Martin Kulldorff, Robert L. Davis, Margarette Kolczak, Edwin Lewis, Tracy Lieu & Richard Platt (2011) A Maximized Sequential Probability Ratio Test for Drug and Vaccine Safety Surveillance, Sequential Analysis, 30:1, 58-78.
#'
#' Stephane Mikael Bottine (2015). SPRT: Wald's Sequential Probability Ratio Test. R package version 1.0. https://CRAN.R-project.org/package=SPRT
#'
#' @export
sprt <- function (df, ...) {
  UseMethod("sprt", df)
}

#' @describeIn sprt SPRT on mds_ts data
#' @export
sprt.mds_ts <- function(
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
  sprt.default(out, analysis_of=name, ...)
}

#' @describeIn sprt SPRT on general data
#' @export
sprt.default <- function(
  df,
  analysis_of=NA,
  eval_period=NULL,
  obs_period=1,
  h0=NULL,
  h1=NULL,
  relative_risk=NULL,
  distribution="poisson",
  alpha=0.05,
  beta=0.20,
  h1_type="greater",
  ...
){
  input_param_checker(df, "data.frame")
  input_param_checker(c("time", "event"), check_names=df)
  input_param_checker(eval_period, "numeric", null_ok=T, max_length=1)
  input_param_checker(obs_period, "numeric", null_ok=F, max_length=1)
  input_param_checker(h0, "numeric", null_ok=T, max_length=1)
  input_param_checker(h1, "numeric", null_ok=T, max_length=1)
  input_param_checker(relative_risk, "numeric", null_ok=T, max_length=1)
  input_param_checker(distribution, "character", null_ok=F, max_length=1)
  input_param_checker(alpha, "numeric", null_ok=F, max_length=1)
  input_param_checker(beta, "numeric", null_ok=F, max_length=1)
  input_param_checker(h1_type, "character", null_ok=F, max_length=1)
  if (!is.null(eval_period)){
    if (eval_period < 2) stop("eval_period must be 2 or greater")
    if (eval_period %% 1 != 0) stop("eval_period must be an integer")
    if (obs_period > eval_period){
      stop("obs_period cannot be greater than eval_period")}
    if (is.null(h0) & obs_period == eval_period){
      stop("obs_period cannot equal eval_period unless h0 is declared")}
  }
  if (obs_period %% 1 != 0 | obs_period <= 0){
    stop("obs_period must be a positive integer")
  }
  if (obs_period > nrow(df)){
    stop("obs_period cannot be longer than rows in df (number of times)")}
  if (!is.null(h0)){ if(h0 <=0) stop("h0 must be >0")}
  if (is.null(h0) & obs_period == nrow(df)){
    stop("obs_period cannot equal rows in df unless h0 is declared.")
  }
  if (!is.null(h1)){ if (h1 <= 0) stop("h1 must be >0")}
  if (!is.null(relative_risk)){
    if (relative_risk <= 0) stop("relative_risk must be >0")
    if (relative_risk == 1) stop("relative_risk cannot be 1")
  }
  if (is.null(h1) & is.null(relative_risk)){
    stop("Either h1 or relative_risk must be specified")}
  if (!distribution %in% c("poisson", "normal")){
    stop("distribution must be either \'poisson\' or \'normal\'")}
  if (alpha <= 0 | alpha >= 1) stop("alpha must be in range (0, 1)")
  if (beta <= 0 | beta >= 1) stop("alpha must be in range (0, 1)")
  if (!h1_type %in% c("greater", "less", "not equal")){
    stop("distribution must be either \'greater\', \'less\', or \'not equal\'")}

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
  if(nrow(df) < 2){
    rr <- NA
    rs <- stats::setNames(F, ">=2 time periods required")
  } else if(is.null(h0) & (obs_period == eval_period)){
    rr <- NA
    rs <- stats::setNames(F, paste("h0 cannot be inferred from data if",
                                   "all data is in observation period"))
  } else{
    # If all conditions are met, run SPRT test
    # Calculate n and k
    n <- obs_period
    k <- sum(df$event[(tlen - n + 1):tlen])
    # Calculate h0 if null
    if (is.null(h0)){
      h0 <- mean(df$event[1:(tlen - n)])
    }
    # Calculate h1 if relative_risk is present
    if (!is.null(relative_risk)) h1 <- relative_risk * h0
    # Wald Boundaries
    bU <- log((1 - beta) / alpha)
    bL <- log(beta / (1 - alpha))
    # Log-likelihood function coefficients
    if (distribution == "poisson"){
      kc <- log(h1) - log(h0)
      nc <- -1 * (h1 - h0)
    } else if (distribution == "normal"){
      kc <- h1 - h0
      nc <- -1 * (((h1 ^ 2) / 2) - ((h0 ^ 2) / 2))
    }
    # Generalized log-likelihood function
    gllf <- function(nc, kc) {
      function(n, k) {
        n * nc + k * kc
      }
    }
    # Log-likelihood function
    llf <- gllf(nc, kc)
    # SPRT
    llr <- llf(n, k)

    # Save output parameters
    if (h1_type == "greater"){
      sig <- llr >= bU
      thresh <- bU
      hyp <- "LLR >= upper rejection level"
    } else if (h1_type == "less"){
      sig <- llr <= bL
      thresh <- bL
      hyp <- "LLR <= lower rejection level"
    } else if (h1_type == "not equal"){
      sig <- (llr >= bU) | (llr <= bL)
      thresh <- c(bL, bU)
      hyp <- "LLR exceeds indifference region"
    }

    rr <- list(statistic=llr,
               lcl=bL,
               ucl=bU,
               p=stats::setNames(NA, "p-value"),
               signal=sig,
               signal_threshold=thresh,
               h0=h0,
               h1=h1)
    rs <- stats::setNames(T, "Success")

    capitalize <- function(x){
      paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
    }
    # Return test
    out <- list(test_name=paste(capitalize(distribution), "SPRT"),
                analysis_of=analysis_of,
                status=rs,
                result=rr,
                params=list(test_hyp=hyp,
                            eval_period=eval_period,
                            obs_period=obs_period,
                            h1_source=ifelse(is.null(relative_risk), "h1",
                                             "relative_risk"),
                            alpha=alpha,
                            beta=beta),
                data=rd)
    class(out) <- append(class(out), "mdsstat_test")
    return(out)
  }
}
