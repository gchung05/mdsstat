#' Empirical Bayes Gamma-Poisson Shrinker
#'
#' Test on device-events using William DuMouchel's Empirical Bayes
#' Gamma-Poisson Shrinker. From
#' the family of disproportionality analyses (DPA) used to generate signals of
#' disproportionate reporting (SDRs).
#'
#' \code{null_ratio} and \code{cred_interval} are used together to establish the
#' signal criteria. The \code{null_ratio} is conceptually similar to the
#' relative reporting ratio under a null hypothesis of no signal. Common values
#' are \code{1} and, more conservatively (fewer false signals), \code{2}. The
#' \code{cred_interval} is the posterior credibility interval used to test for a
#' signal. A value of \code{0.90} returns the 5% and 95% credible bounds and
#' tests if the lower bound exceeds \code{null_ratio}. Effectively,
#' \code{cred_interval=0.90} conducts the well-known EB05 test.
#'
#' \code{init_prior} specifies the initial guess for the 5 parameters of the
#' prior gamma mixture distribution as described in DuMouchel (1999, Eqs. 4, 7)
#' in the sequence: alpha1, beta1, alpha2, beta2, p. \code{gamma_lower}
#' specifies the optimization lower bound for the two alphas and betas.
#' \code{gamma_upper} specifies similarly the upper bound. The initial guess,
#' upper and lower bounds are fed into PORT optimization using the
#' \code{stats::nlminb()} routine.
#'
#' \code{cont_adj} provides the option to allow \code{gps()} to proceed running,
#' however this is done at the user's discretion because there are adverse
#' effects of adding a positive integer to every cell of the contingency table.
#' By default, \code{gps()} runs with 0 in the C cell only, but not in A, B, or
#' D. It has been suggested that 0.5 may be an appropriate value. However,
#' values <1 have been shown to be unstable using box-constrained PORT
#' optimization, which is the only optimization considered in this release.
#' Overall, posterior distribution estimates have been shown to be unstable with
#' very low or 0 count cells.
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
#' ratio (RR), used with \code{cred_interval} to establish the signal status.
#' This \code{null_ratio} is saved in the output as the signal threshold. See
#' details for more.
#'
#' Default: \code{1} indicates a null RR of 1 and tests if the lower bound of
#' the \code{cred_interval} exceeds \code{1}.
#'
#' @param cred_interval Numeric value between 0 and 1 representing the width of
#' the Bayesian posterior credible interval, where the lower bound of the
#' interval is assessed against the \code{null_ratio}. The interval bounds are
#' returned as the lcl and ucl. See details for more.
#'
#' Default: \code{0.90} indicates a 90\% credible interval with bounds at 5\% and
#' 95\%. The signal test is against the lower 5\% bound, effectively creating the
#' EB05 test.
#'
#' @param init_prior A numeric vector of length 5 representing the
#' initialization parameters for the prior gamma mixture distribution in this
#' order: \code{alpha1, beta1, alpha2, beta2, p}. See details for more.
#'
#' Default: \code{c(.2, .02, 2, 4, 1/3)} as suggested in openEBGM package
#' v0.7.0.
#'
#' @param gamma_lower Positive mumeric value representing the lower bound for
#' the two alphas and betas of the prior during PORT optimization.
#'
#' Default: \code{1e-5} is a value suggested in openEBGM package v0.7.0.
#'
#' @param gamma_upper Positive mumeric value representing the upper bound for
#' the two alphas and betas of the prior during PORT optimization.
#'
#' Default: \code{20} is a value suggested in openEBGM package v0.7.0.
#'
#' @param quantiles Vector of quantiles between 0 and 1. \code{gps()} will
#' return an equal length vector of estimated empirical Bayes quantiles from the
#' posterior distribution. Specify \code{quantiles=NULL} if no quantiles are
#' desired.
#'
#' Default: \code{c(.05, .95)} corresponds to the 5\% (EB05) and 95\% (EB95)
#' quantiles.
#'
#' @param cont_adj Positive integer representing the continuity
#' adjustment to be added to each cell of the 2x2 contingency table. A value
#' greater than 0 allows for contingency tables with 0 cells to run the
#' algorithm. Adding a continuity adjustment will adversely affect the algorithm
#' estimates, user discretion is advised. See details for more.
#'
#' Default: \code{0} adds zero to each cell, thus an unadjusted table.
#'
#' @param ... Further arguments passed onto \code{gps} methods
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
#' a1 <- gps(data)
#' # Example using an mds_ts object
#' a2 <- gps(mds_ts[[3]])
#'
#' @references
#' DuMouchel W. Bayesian data mining in large frequency tables, with an application to the FDA spontaneous reporting system. The American Statistician, 53(3):177-190, August 1999.
#'
#' Ahmed I, Poncet A. PhViD: PharmacoVigilance Signal Detection, 2016. R package version 1.0.8.
#'
#' Ihrie J, Canida T. openEBGM: EBGM Scores for Mining Large Contingency Tables, 2018. R package version 0.7.0.
#'
#' @export
gps <- function (df, ...) {
  UseMethod("gps", df)
}

#' @describeIn gps GPS on mds_ts data
#' @export
gps.mds_ts <- function(
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
  gps.default(out, analysis_of=name, ...)
}

#' @describeIn gps GPS on general data
#' @export
gps.default <- function(
  df,
  analysis_of=NA,
  eval_period=1,
  null_ratio=1,
  cred_interval=0.90,
  init_prior=c(.2, .02, 2, 4, 1/3),
  gamma_lower=1e-5,
  gamma_upper=20,
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
  input_param_checker(cred_interval, "numeric", null_ok=F, max_length=1)
  input_param_checker(init_prior, "numeric", null_ok=F, max_length=5)
  input_param_checker(gamma_lower, "numeric", null_ok=F, max_length=1)
  input_param_checker(gamma_upper, "numeric", null_ok=F, max_length=1)
  input_param_checker(quantiles, "numeric", null_ok=T)
  input_param_checker(cont_adj, "numeric", null_ok=F, max_length=1)
  if (eval_period %% 1 != 0) stop("eval_period must be an integer")
  if (null_ratio < 1) stop("null_ratio must be 1 or greater")
  if (cred_interval <= 0 | cred_interval >= 1){
    stop("cred_interval must be in range (0, 1)")}
  if (length(init_prior) != 5) stop("init_prior must be of length 5")
  if (init_prior[5] <= 0 | init_prior[5] >= 1){
    stop("init_prior p must be in range (0, 1)")}
  if (gamma_lower >= gamma_upper){
    stop("gamma_lower must be less than gamma_upper")}
  if (gamma_lower <= 0 | gamma_upper <= 0){
    stop("gamma_lower and gamma_upper must be positive")
  }
  if (any(init_prior[1:4] < gamma_lower)){
    stop("init_prior alphas and betas cannot be lower than gamma_lower")}
  if (any(init_prior[1:4] > gamma_upper)){
    stop("init_prior alphas and betas cannot be greater than gamma_upper")}
  if (any(quantiles <= 0) | any(quantiles >= 1)){
    stop("quantiles must be in range (0, 1)")
  }
  if (cont_adj < 0 | cont_adj %% 1 != 0){
    stop("cont_adj must be a positive integer")}

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
  if (any(df[, c("nB", "nD")] == 0)){
    rr <- NA
    rs <- stats::setNames(F, "contingency cells B or D have zero counts")
  } else if (any(df[, c("nA")] == 0)){
    rr <- NA
    rs <- stats::setNames(F, "contingency cell A is zero")
  } else{
    # If all conditions are met, run GPS test
    # Observed and expected
    N <- unlist(df[, c2x2])
    E <- E2x2(N)
    # Prior mixture function
    prior_f <- function(p, N, E){
      sum(-log((p[5] * stats::dnbinom(N,
                                      size=p[1],
                                      prob=p[2] / (p[2] + E)) +
                  (1 - p[5]) * stats::dnbinom(N,
                                              size=p[3],
                                              prob=p[4] / (p[4] + E)))))
    }
    # Prior estimation optimization using PORT
    optim <- suppressWarnings(
      stats::nlminb(start=init_prior, objective=prior_f,
                    lower=c(rep(gamma_lower, 5)),
                    upper=c(rep(gamma_upper, 4), .9999),
                    N=N, E=E))
    prior <- optim$par
    # Posterior mixture Q
    p_gamma1 <- prior[5] * stats::dnbinom(N[1],
                                          size=prior[1],
                                          prob=prior[2] / (prior[2] + E[1]))
    p_gamma_2 <- (1 - prior[5]) * stats::dnbinom(N[1],
                                                 size=prior[3],
                                                 prob=prior[4] / (prior[4] + E[1]))
    Q <- p_gamma1 / (p_gamma1 + p_gamma_2)
    # Posterior expectation of lambda
    pe <- (Q * (digamma(prior[1] + N[1]) - log(prior[2] + E[1])) +
             (1 - Q) * (digamma(prior[3] + N[1]) - log(prior[4] + E[1]))) / log(2)
    ebgm <- 2 ^ pe
    # Quantile estimation functions
    quantEst <- function(thresh, Q, a1, b1, a2, b2){
      l <- rep(-1e5, length(Q))
      u <- rep(1e5, length(Q))
      x <- rep(1, length(Q))
      Cout <- fCout(x, thresh, Q, a1, b1, a2, b2)
      while (max(round(Cout * 1e4)) != 0){
        S <- sign(Cout)
        xnow <- (1 + S) / 2 * ((x + l) / 2) + (1 - S) / 2 * ((x + u) / 2)
        u <- (1 + S) / 2 * x + (1 - S) / 2 * u
        l <- (1 + S) / 2 * l + (1 - S) / 2 * x
        x <- xnow
        Cout <- fCout(x, thresh, Q, a1, b1, a2, b2)
      }
      x
    }
    fCout <- function(p, thresh, Q, a1, b1, a2, b2){
      Q * stats::pgamma(p, shape=a1, rate=b1) +
        (1 - Q) * stats::pgamma(p, shape=a2, rate=b2) - thresh
    }
    # Estimate quantiles
    qEst <- numeric()
    if (length(quantiles) > 0){
      if (all(!is.na(quantiles))){
        for (i in c(1:length(quantiles))){
          qEst[i] <- quantEst(quantiles[i], Q,
                              prior[1] + N[1],
                              prior[2] + E[1],
                              prior[3] + N[1],
                              prior[4] + E[1])
        }
        stats::setNames(qEst, quantiles)
      }
    }
    # Credibility interval limits
    cl_l <- round((1 - cred_interval) / 2, 3)
    cl_u <- 1 - cl_l
    cl <- c(quantEst(cl_l, Q,
                   prior[1] + N[1],
                   prior[2] + E[1],
                   prior[3] + N[1],
                   prior[4] + E[1]),
            quantEst(cl_u, Q,
                     prior[1] + N[1],
                     prior[2] + E[1],
                     prior[3] + N[1],
                     prior[4] + E[1]))
    sig <- cl[1] > null_ratio
    hyp <- paste0(1e2 * cl_l, "% quantile of the posterior distribution > ",
                  null_ratio)

    rr <- list(statistic=stats::setNames(ebgm, "EBGM"),
               lcl=cl[1],
               ucl=cl[2],
               p=stats::setNames(NA, "p-value"),
               signal=sig,
               signal_threshold=stats::setNames(null_ratio, "null ratio"),
               quantiles=qEst)
    rs <- stats::setNames(T, "Success")
  }

  # Return test
  out <- list(test_name="Gamma Poisson Shrinker",
              analysis_of=analysis_of,
              status=rs,
              result=rr,
              params=list(test_hyp=hyp,
                          eval_period=eval_period,
                          null_ratio=null_ratio,
                          cred_interval=cred_interval,
                          init_prior=init_prior,
                          gamma_lower=gamma_lower,
                          gamma_upper=gamma_upper,
                          cont_adj=cont_adj),
              data=rd)
  class(out) <- append(class(out), "mdsstat_test")
  return(out)
}
