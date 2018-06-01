#' Coerce mdsstat Test to 1-Row Data Frame
#'
#' Coerce an mdsstat test (class \code{mdsstat_test}) to a 1-row data frame.
#'
#' @param df Required input object of class \code{mdsstat_test}
#' @return 1-row data frame (class \code{mdsstat_df}) summarizing the test.
#' @examples
#' test2row(prr(mds_ts[[3]]))
#' @export
test2row <- function(
  df
){
  input_param_checker(df, "mdsstat_test")

  out <- data.frame(test_name=df$test_name,
                    analysis_of=df$analysis_of,
                    run_status=df$status,
                    run_msg=names(df$status),
                    reference_time=as.character(df$data$reference_time),
                    eval_period=df$params$eval_period,
                    test_hyp=df$params$test_hyp,
                    test_params=I(list(df$params[!names(df$params) %in%
                                                   c("test_hyp",
                                                     "eval_period")])))
  if (all(is.na(df$result))){
    out <- cbind(out,
                 data.frame(signal=NA,
                            signal_threshold=I(list(NA)),
                            stat=I(list(NA)),
                            stat_name=I(list(NA)),
                            stat_lcl=I(list(NA)),
                            stat_ucl=I(list(NA)),
                            p_value=NA,
                            stat_addtl=I(list(NA))))
  } else{
    out <- cbind(out,
                 data.frame(signal=df$result$signal,
                            signal_threshold=I(list(df$result$signal_threshold)),
                            stat=I(list(df$result$statistic)),
                            stat_name=I(list(names(df$result$statistic))),
                            stat_lcl=I(list(names(df$result$lcl))),
                            stat_ucl=I(list(names(df$result$ucl))),
                            p_value=df$result$p,
                            stat_addtl=I(list(df$result[!names(df$result) %in%
                                                          c("signal",
                                                            "signal_threshold",
                                                            "statistic",
                                                            "lcl", "ucl", "p")]))))
  }
  rownames(out) <- c()
  class(out) <- append(class(out), "mdsstat_df")
  return(out)
}


#' Set List of Algorithms to Run
#'
#' Define any number of algorithms with various parameter settings and save as a
#' reusable set of instructions.
#'
#' @param algos Required named list of \code{mdsstat} algorithms to run. Each
#' named list element must be a single list of parameter values for the
#' algorithm named. The list of parameters may be an empty list (indicating the
#' default values) and must not contain the first parameter \code{df}. See
#' details and examples for more.
#' @return Validated list of instructions that may be used in the
#' \code{\link{run_algos}} function.
#' @details Valid list names (algorithm function names in \code{mdsstat}) are
#' currently \code{\link{shewhart}}, \code{\link{prr}}, and
#' \code{\link{poisson_rare}}. Each algorithm may be named multiple times (to
#' allow running of multiple parameter settings). Do not specify the \code{df}
#' parameter.
#' @examples
#' x <- list(prr=list(),
#'   shewhart=list(),
#'   shewhart=list(ts_event=c(Rate="rate"), we_rule=2L),
#'   poisson_rare=list(p_rate=0.3))
#' define_algos(x)
#' @export
define_algos <- function(
  algos
){
  input_param_checker(algos, "list")

  # Each list item must be an algorithm in mdsstat
  algolist <- ls("package:mdsstat")[grepl("\\.mds_ts$", ls("package:mdsstat"))]
  algolist <- gsub("\\.mds_ts$", "", algolist)
  if (!all(names(algos) %in% algolist)){
    notalgo <- names(algos)[!names(algos) %in% algolist]
    stop(paste(paste(notalgo, collapse=", "),
               ifelse(length(notalgo) > 1, "are", "is"), "not",
               ifelse(length(notalgo) > 1, "mdsstat", "an mdsstat"),
               ifelse(length(notalgo) > 1, "algorithms", "algorithm")))
  }

  # Each algorithm list element must itself be a list
  if (!all(sapply(algos, class) == "list")){
    notlist <- names(algos)[sapply(algos, class) != "list"]
    stop(paste(paste(notlist, collapse=", "), "must be",
               ifelse(length(notlist) > 1, "lists", "a list")))
  }

  # All list elements in each algorithm must be parameter names and not be df
  for (i in 1:length(algos)){
    args <- unique(unlist(
      sapply(c(names(algos)[i], suppressWarnings(methods(names(algos)[i]))), formalArgs)))
    args <- args[!args %in% c("df", "...")]
    if (any(names(algos[[i]]) %in% "df")){
      stop(paste("Do not specify df parameter in", names(algos)[i]))
    } else if (!all(names(algos[[i]]) %in% args)){
      notargs <- names(algos[[i]])[!names(algos[[i]]) %in% args]
      stop(paste(paste(notargs, collapse=", "),
                 ifelse(length(notargs) > 1, "are", "is"),
                 ifelse(length(notargs) > 1, "not", "not a"),
                 names(algos)[i],
                 ifelse(length(notargs) > 1, "parameters", "parameter")))
    }
  }

  class(algos) <- append(class(algos), "mdsstat_da")
  return(algos)
}


#' Run Multiple Algorithms
#'
#' Run a multiple number of \code{mdsstat} algorithms on a single input dataset.
#'
#' @param data Required input dataset. Note that the dataset must satisfy the
#' dataset (parameter \code{df}) requirements for each algorithm specified.
#' @param algos Input list of algorithms to run. Must be a list generated by
#' \code{\link{define_algos}}.
#' @param dataframe Logical on whether to output results as a
#' \code{mdsstat_tests} data frame. If \code{FALSE}, will output as a list of
#' \code{mdsstat_test} lists.
#'
#' Default: \code{TRUE}
#'
#' @return A \code{mdsstat_tests} data frame or list of \code{mdsstat_test}
#' lists with the results of the algorithm runs.
#' @examples
#' data <- data.frame(time=c(1:25),
#'                    nA=as.integer(stats::rnorm(25, 25, 5)),
#'                    nB=as.integer(stats::rnorm(25, 50, 5)),
#'                    nC=as.integer(stats::rnorm(25, 100, 25)),
#'                    nD=as.integer(stats::rnorm(25, 200, 25)))
#' x <- list(prr=list(),
#'   shewhart=list(),
#'   shewhart=list(ts_event=c(Rate="rate"), we_rule=2L),
#'   poisson_rare=list(p_rate=0.3))
#' algos <- define_algos(x)
#' run_algos(data, algos)
#' run_algos(data, algos, F)
#' @export
run_algos <- function(
  data,
  algos,
  dataframe=T
){
  input_param_checker(data, "data.frame")
  input_param_checker(algos, "mdsstat_da")
  input_param_checker(dataframe, "logical")

  # I AM HERE!!!

  if (dataframe) class(stats) <- append(class(stats), "mdsstat_tests")
  return(stats)
}

# do.call(runif,c(list(n=100),params))

# x <- list(prr=list(),
#           shewhart=list(ts_event=c(Rate="rate"),
#                         we_rule=2L),
#           poisson_rare=list(p_rate=0.3))
