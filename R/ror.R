#' Reporting Odds Ratio
#'
#' Test on device-events using the reporting odds ratio (ROR). From
#' the family of disproportionality analyses (DPA) used to generate signals of
#' disproportionate reporting (SDRs).
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
#' Default: \code{1L} considers only the most recent time in \code{df}.
#'
#' Example: \code{12L} sums over the last 12 time periods to create the 2x2
#' contingency table.
#'
#' @param null_ratio Numeric ROR value representing the null hypothesis, used
#' with \code{alpha} to establish the signal status and the p-value.
#'
#' Default: \code{1} indicates a null hypothesis of ROR=1 and tests if the
#' actual ROR is greater than 1.
#'
#' @param alpha Numeric value representing the statistical alpha used to
#' establish the signal status.
#'
#' Default: \code{0.05} corresponds to the standard alpha value of 5\%.
#'
#' @param cont_adj Numeric value 0 or greater representing the continuity
#' adjustment to be added to each cell of the 2x2 contingency table. A value
#' greater than 0 allows for contingency tables with 0 cells to run the
#' algorithm. A typical non-zero value is 0.5.
#'
#' Default: \code{0} adds zero to each cell, thus an unadjusted table. If any
#' cell of the 2x2 is 0, the algorithm will not run.
#'
#' @param ... Further arguments passed onto \code{ror} methods
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
#' a1 <- ror(data)
#' # Example using an mds_ts object
#' a2 <- ror(mds_ts[[3]])
#'
#' @references
#' Stricker BH, Tijssen JG. Serum sickness-like reactions to cefaclor. J Clin Epidemiol. 1992;45(10):1177-84.
#'
#' Bohm R, Klein H.-J. (v2018-10-16). Primer on Disportionality Analysis. OpenVigil http://openvigil.sourcefourge.net/doc/DPA.pdf
#'
#' @export
ror <- function (df, ...) {
  UseMethod("ror", df)
}

#' @describeIn ror ROR on mds_ts data
#' @export
ror.mds_ts <- function(
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
  ror.default(out, analysis_of=name, ...)
}

#' @describeIn ror ROR on general data
#' @export
ror.default <- function(
  df,
  analysis_of=NA,
  eval_period=1L,
  null_ratio=1,
  alpha=0.05,
  cont_adj=0,
  ...
){
  # Contingency table primary variables
  c2x2 <- c("nA", "nB", "nC", "nD")

  input_param_checker(df, "data.frame")
  input_param_checker(c("time", c2x2), check_names=df)
  input_param_checker(eval_period, "integer")
  input_param_checker(null_ratio, "numeric")
  input_param_checker(alpha, "numeric")
  input_param_checker(cont_adj, "numeric")
  if (null_ratio < 1) stop("null_ratio must be 1 or greater")
  if (alpha <= 0 | alpha >= 1) stop("alpha must be in range (0, 1)")
  if (cont_adj < 0) stop("cont_adj must be 0 or greater")

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
  if(any(df[, c2x2] == 0)){
    rr <- NA
    rs <- stats::setNames(F, "contingency table has zero counts")
  } else{
    # If all conditions are met, run ROR test
    # Calculate ROR
    stat <- (df$nA / df$nB) / (df$nC / df$nD)
    s <- sqrt((1 / df$nA) + (1 / df$nB) + (1 / df$nC) + (1 / df$nD))
    # Establish confidence limits
    z <- stats::qnorm(1 - (alpha / 2))
    cl <- c(exp(log(stat) - z * s), exp(log(stat) + z * s))
    p <- min(stats::pnorm((log(null_ratio) - log(stat)) / s) * 2, 1)
    # Determine signal & hypothesis
    sig <- p <= alpha
    hyp <- paste0("Two-sided test at alpha=", alpha, " of ROR > ", null_ratio)

    rr <- list(statistic=stats::setNames(stat, "ROR"),
               lcl=cl[1],
               ucl=cl[2],
               p=p,
               signal=sig,
               signal_threshold=stats::setNames(alpha, "critical p-value"),
               sigma=exp(s))
    rs <- stats::setNames(T, "Success")
  }

  # Return test
  out <- list(test_name="Reporting Odds Ratio",
              analysis_of=analysis_of,
              status=rs,
              result=rr,
              params=list(test_hyp=hyp,
                          eval_period=eval_period,
                          null_ratio=null_ratio,
                          alpha=alpha,
                          cont_adj=cont_adj),
              data=rd)
  class(out) <- append(class(out), "mdsstat_test")
  return(out)
}
