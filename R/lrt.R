#' Likelihood Ratio Rest
#'
#' Test on device-events using the Likelihood Ratio Test, originally proposed by
#' Huang & Tiwari (2012). From
#' the family of disproportionality analyses (DPA) used to generate signals of
#' disproportionate reporting (SDRs).
#'
#' This is an implementation of the "Regular LRT" per the revised 2019
#' article by Huang & Tiwari. It assumes a test on a single event of interest
#' where all other events & devices are collapsed, effectively testing a 2x2
#' table only. Therefore this is a test on the significance of the likelihood
#' ratio instead of the maximum likelihood over \code{n} events.
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
#' @param alpha Alpha or Type-I error rate in the range (0, 1), used to
#' determine signal status.
#' It is the threshold for determining if the observed reporting rate is
#' greater than the expected based on Monte Carlo simulations of the null.
#'
#' Default: \code{0.05} is an alpha level of 0.05 or 5\%.
#'
#' @param mc_sample Number of Monte Carlo samples for constructing the
#' null distribution based on empirical data. Lowest recommended is 1000.
#' Increasing iterations also increases p-value precision.
#'
#' Default: \code{10000} uses 10000 bootstrap iterations.
#'
#' @param ... Further arguments passed onto \code{lrt} methods
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
#' a1 <- lrt(data)
#' # Example using an mds_ts object
#' a2 <- lrt(mds_ts[[3]])
#'
#' @references
#' Huang L, Zalkikar J, Tiwari RC. A Likelihood Ratio Test Based Method for Signal Detection With Application to FDA’s Drug Safety Data. Journal of the American Statistical Association, 2011, Volume 106, Issue 496, 1230-1241.
#'
#' Huang L, Zalkikar J, Tiwari RC. Likelihood-Ratio-Test Methods for Drug Safety Signal Detection from Multiple Clinical Datasets. Comput Math Methods Med. 2019, PMC6399568.
#'
#' @export
lrt <- function (df, ...) {
  UseMethod("lrt", df)
}

#' @describeIn lrt LRT on mds_ts data
#' @export
lrt.mds_ts <- function(
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
                  attributes(df)$dpa_detail$nA)
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
  lrt.default(out, analysis_of=name, ...)
}

#' @describeIn lrt LRT on general data
#' @export
lrt.default <- function(
  df,
  analysis_of=NA,
  eval_period=1,
  alpha=0.05,
  mc_sample=10000,
  ...
){
  # Contingency table primary variables
  c2x2 <- c("nA", "nB", "nC", "nD")

  input_param_checker(df, "data.frame")
  input_param_checker(c("time", c2x2), check_names=df)
  input_param_checker(eval_period, "numeric", null_ok=F, max_length=1)
  input_param_checker(alpha, "numeric", null_ok=F, max_length=1)
  input_param_checker(mc_sample, "numeric", null_ok=F, max_length=1)
  if (eval_period %% 1 != 0) stop("eval_period must be an integer")
  if (alpha <= 0 | alpha >= 1) stop("alpha must be in range (0, 1)")
  if (mc_sample %% 1 != 0) stop("mc_sample must be an integer")
  if (mc_sample < 1000) stop("mc_sample should be >=1000")

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
    # Observed and expected
    N <- unlist(df[, c2x2])
    E <- E2x2(N)
    # If all conditions are met, run LRT test
    nij <- N[1]
    Eij <- E[1]
    n_j <- sum(N[c(1, 3)])
    P_ <- sum(N) # Exposure availability not assumed, use n__ instead
    pij <- N[1] / sum(N[c(1, 2)])
    qij <- N[3] / sum(N[c(3, 4)])
    # Log likelihood ratio
    llr <- function(nij, Eij, n_j, P_){
      nij * (log(nij) - log(Eij)) +
        (n_j - nij) * (log(n_j - nij) - log(n_j - Eij)) -
        n_j * (log(n_j) - log(P_))
    }
    llr_hat <- llr(nij, Eij, n_j, P_)
    # Monte Carlo null distribution
    mc_llr <- mat.or.vec(mc_sample, 1)
    for(i in 1:mc_sample){
      samp <- stats::rmultinom(1, n_j, c(N[1] / P_, N[3] / P_))
      samp_N <- c(samp[1], N[2], samp[2], N[4])
      samp_E <- E2x2(samp_N)
      mc_llr[i] <- llr(samp[1], samp_E[1], n_j, sum(samp_N))
    }
    pvalue <- sum(llr_hat < mc_llr) / mc_sample
    sig <- pvalue <= alpha
    hyp <- paste0("Device-event reporting rate (RR) (",
                  round(pij, 3), ") > comparator RR (",
                  round(qij, 3), ") at alpha=", alpha)

    # List of output values
    rr <- list(statistic=stats::setNames(llr_hat, "Log-Likelihood Ratio"),
               lcl=NA,
               ucl=NA,
               p=stats::setNames(pvalue, "p-value"),
               signal=sig,
               signal_threshold=stats::setNames(alpha, "alpha"))
    rs <- stats::setNames(T, "Success")
  }

  # Return test
  out <- list(test_name="LRT",
              analysis_of=analysis_of,
              status=rs,
              result=rr,
              params=list(test_hyp=hyp,
                          eval_period=eval_period,
                          alpha=alpha,
                          mc_sample=mc_sample),
              data=rd)
  class(out) <- append(class(out), "mdsstat_test")
  return(out)
}
