#' Bayesian Confidence Propagation Neural Network
#'
#' Test on device-events using a one-layer BCPNN as proposed by Bate et al
#' (1998), which assumes beta-distributed probabilities and a normal
#' approximation of the variance of the information component (IC). From
#' the family of disproportionality analyses (DPA) used to generate signals of
#' disproportionate reporting (SDRs).
#'
#' \code{null_ratio} and \code{conf_interval} are used together to establish the
#' signal criteria. The \code{null_ratio} is conceptually similar to the
#' relative reporting ratio under a null hypothesis of no signal. Common values
#' are \code{1} and, more conservatively (fewer false signals), \code{2}. The
#' \code{conf_interval} is the IC confidence interval used to test for a
#' signal. A value of \code{0.90} returns the 5% and 95% confidence bounds and
#' tests if the lower bound exceeds \code{null_ratio}. Effectively,
#' \code{conf_interval=0.90} conducts a one-sided test at the conventional 0.05
#' alpha level.
#'
#' \code{cont_adj} provides the option to allow \code{bcpnn()} to proceed running,
#' however this is done at the user's discretion because there are adverse
#' effects of adding a positive number to every cell of the contingency table.
#' By default, \code{bcpnn()} allows 0 in all cells except D.
#' It has been suggested that 0.5 may be an appropriate value. However, the
#' user is cautioned that interpretation may be compromised by adding continuity
#' adjustments.
#'
#' For parameter \code{ts_event}, in the uncommon case where the
#' device-event count (Cell A) variable is not \code{"nA"}, the name of the
#' variable may be specified here. Note that the remaining 3 cells of the 2x2
#' contingency table (Cells B, C, D) must be the variables \code{"nB"},
#' \code{"nC"}, and \code{"nD"} respectively in \code{df}. A named character
#' vector may be used where the name is the English description of what was
#' analyzed. Note that if the parameter \code{analysis_of} is specified, it will
#' override this name. Example: \code{ts_event=c("Count of Bone Cement
#' Leakages"="event_count")}
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
#' @param eval_period Required positive integer indicating the number of unique
#' times counting in reverse chronological order to sum over to create the 2x2
#' contingency table.
#'
#' Default: \code{1} considers only the most recent time in \code{df}.
#'
#' Example: \code{12} sums over the last 12 time periods to create the 2x2
#' contingency table.
#'
#' @param null_ratio Numeric value representing the null relative reporting
#' ratio (RR), used with \code{conf_interval} to establish the signal status.
#' This \code{null_ratio} is saved in the output as the signal threshold. See
#' details for more.
#'
#' Default: \code{1} indicates a null RR of 1 and tests if the lower bound of
#' the IC \code{conf_interval} exceeds \code{1}.
#'
#' @param conf_interval Numeric value between 0 and 1 representing the width of
#' the IC confidence interval, where the lower bound of the
#' interval is assessed against the \code{null_ratio}. The interval bounds are
#' returned as the lcl and ucl. See details for more.
#'
#' Default: \code{0.90} indicates a 90\% confidence interval with bounds at 5\%
#' and 95\%. The signal test is against the lower 5\% bound, effectively
#' creating a one-sided test at the conventional 0.05 alpha level.
#'
#' @param quantiles Vector of quantiles between 0 and 1. \code{bcpnn()} will
#' return an equal length vector of Gaussian-approximated quantiles from the
#' posterior distribution of the IC. Specify \code{quantiles=NULL} if no
#' quantiles are desired.
#'
#' Default: \code{c(.05, .95)} corresponds to the 5\% and 95\% quantiles of the
#' IC.
#'
#' @param cont_adj Non-negative number representing the continuity
#' adjustment to be added to each cell of the 2x2 contingency table. A value
#' greater than 0 allows for contingency tables with a 0 in the D cell to run
#' the algorithm. Adding a continuity adjustment will adversely affect the
#' algorithm estimates, user discretion is advised. See details for more.
#'
#' Default: \code{0} adds zero to each cell, thus an unadjusted table.
#'
#' @param ... Further arguments passed onto \code{bcpnn} methods
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
#' # Basic Example
#' data <- data.frame(time=c(1:25),
#'                    nA=as.integer(stats::rnorm(25, 25, 5)),
#'                    nB=as.integer(stats::rnorm(25, 50, 5)),
#'                    nC=as.integer(stats::rnorm(25, 100, 25)),
#'                    nD=as.integer(stats::rnorm(25, 200, 25)))
#' a1 <- bcpnn(data)
#' # Example using an mds_ts object
#' a2 <- bcpnn(mds_ts[[3]])
#'
#' @references
#' Bate A, Lindquist M, et al. A Bayesian Neural Network Method for Adverse Drug Reaction Signal Generation. European Journal of Clinical Pharmacology, 1998, 54, 315-321.
#'
#' Ahmed I, Poncet A. PhViD: PharmacoVigilance Signal Detection, 2016. R package version 1.0.8.
#'
#' Lansner A, Ekeberg Ö. A one-layer feedback artificial neural network with a bayesian learning rule. Int. J. Neural Syst., 1989, 1, 77-87.
#'
#' @export
bcpnn <- function (df, ...) {
  UseMethod("bcpnn", df)
}

#' @describeIn bcpnn BCPNN on mds_ts data
#' @export
bcpnn.mds_ts <- function(
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

  if (attributes(df)$dpa){
    out <- data.frame(time=df$time,
                      nA=df[[ts_event]],
                      nB=df$nB,
                      nC=df$nC,
                      nD=df$nD)
  } else{
    stop("Input mds_ts df does not contain data for disproportionality analysis.")
  }
  bcpnn.default(out, analysis_of=name, ...)
}

#' @describeIn bcpnn BCPNN on general data
#' @export
bcpnn.default <- function(
  df,
  analysis_of=NA,
  eval_period=1,
  null_ratio=1,
  conf_interval=0.90,
  quantiles=c(.05, .95),
  cont_adj=0,
  ...
){
  # Contingency table primary variables
  c2x2 <- c("nA", "nB", "nC", "nD")

  input_param_checker(df, "data.frame")
  input_param_checker(c("time", c2x2), check_names=df)
  input_param_checker(eval_period, "numeric", null_ok=F, max_length=1)
  input_param_checker(null_ratio, "numeric", null_ok=F, max_length=1)
  input_param_checker(conf_interval, "numeric", null_ok=F, max_length=1)
  input_param_checker(quantiles, "numeric", null_ok=T)
  input_param_checker(cont_adj, "numeric", null_ok=F, max_length=1)
  if (eval_period %% 1 != 0) stop("eval_period must be an integer")
  if (null_ratio < 1) stop("null_ratio must be 1 or greater")
  if (conf_interval <= 0 | conf_interval >= 1){
    stop("conf_interval must be in range (0, 1)")}
  if (any(quantiles <= 0) | any(quantiles >= 1)){
    stop("quantiles must be in range (0, 1)")
  }
  if (cont_adj < 0){
    stop("cont_adj must be 0 or greater")}

  # Order by time
  df <- df[order(df$time), ]
  # Restrict to eval_period
  if (!is.null(eval_period)){
    if (eval_period > nrow(df)){
      stop("eval_period cannot be greater than df rows")
    } else if (eval_period < 1){
      stop("eval_period must be greater than 0")
    } else{
      df <- df[c((nrow(df) - eval_period + 1):nrow(df)), ]
      # Sum over eval_period
      timeRange <- range(df$time)
      df <- cbind(data.frame(time_start=timeRange[1],
                             time_end=timeRange[2]),
                  data.frame(t(colSums(df[, c2x2], na.rm=T))))
      # Apply continuity adjustment
      df[, c2x2] <- df[, c2x2] + cont_adj
    }
  }
  # Return data
  tlen <- nrow(df)
  rd <- list(reference_time=timeRange,
             data=df)

  # Check for non-runnable conditions
  hyp <- "Not run"
  if (any(df[, c("nD")] == 0)){
    rr <- NA
    rs <- stats::setNames(F, "contingency cell D is zero")
  } else{
    # If all conditions are met, run BCPNN test
    # Observed and expected
    N <- unlist(df[, c2x2])
    E <- E2x2(N)
    # Row, column, and table marginals
    nd <- c(rep(sum(N[c(1, 2)]), 2), rep(sum(N[c(3, 4)]), 2))
    ne <- rep(c(sum(N[c(1, 3)]), sum(N[c(2, 4)])), 2)
    n <- sum(N)
    # Interim variables
    p1 <- nd + 1
    p2 <- n - nd + 1
    q1 <- ne + 1
    q2 <- n - ne + 1
    r1 <- N + 1
    r2b <- n - N - 1 + (2 + n) ^ 2 / (q1 * p1)
    # Expectation (eIC) and variance (vIC) of the information component
    eIC <- log(2) ^ -1 * (digamma(r1) - digamma(r1 + r2b) -
                            (digamma(p1) - digamma(p1 + p2) +
                               digamma(q1) - digamma(q1 + q2)))
    vIC <- log(2) ^ -2 * (trigamma(r1) - trigamma(r1 + r2b) +
                            (trigamma(p1) - trigamma(p1 + p2) +
                               trigamma(q1) - trigamma(q1 + q2)))
    # Posterior IC estimate
    IC <- exp(eIC)
    # p-value
    p <- stats::pnorm(log(null_ratio), eIC, sqrt(vIC))[1]
    # Confidence interval limits
    cl_l <- round((1 - conf_interval) / 2, 3)
    cl_u <- 1 - cl_l
    cl <- sapply(c(cl_l, cl_u),
                 function(x) exp(stats::qnorm(x, eIC[1], sqrt(vIC[1]))))
    sig <- p <= cl_l
    hyp <- paste0(1e2 * cl_l, "% quantile of the posterior distribution > ",
                  null_ratio)
    # Estimate quantiles
    qEst <- numeric()
    if (length(quantiles) > 0){
      if (all(!is.na(quantiles))){
        qEst <- sapply(quantiles,
                       function(x) exp(stats::qnorm(x, eIC[1], sqrt(vIC[1]))))
        stats::setNames(qEst, quantiles)
      }
    }

    rr <- list(statistic=stats::setNames(IC[1], "IC"),
               lcl=cl[1],
               ucl=cl[2],
               p=stats::setNames(p, "p-value"),
               signal=sig,
               signal_threshold=stats::setNames(null_ratio, "null ratio"),
               quantiles=qEst)
    rs <- stats::setNames(T, "Success")
  }

  # Return test
  out <- list(test_name="BCPNN",
              analysis_of=analysis_of,
              status=rs,
              result=rr,
              params=list(test_hyp=hyp,
                          eval_period=eval_period,
                          null_ratio=null_ratio,
                          conf_interval=conf_interval,
                          cont_adj=cont_adj),
              data=rd)
  class(out) <- append(class(out), "mdsstat_test")
  return(out)
}
