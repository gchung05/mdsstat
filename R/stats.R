#' Coerce mdsstat Test to 1-Row Data Frame
#'
#' Coerce an mdsstat test (class \code{mdsstat_test}) to a 1-row data frame.
#'
#' @param df Required input object of class \code{mdsstat_test}
#' @return 1-row data frame (class \code{mdsstat_df}) summarizing the test.
#' @examples
#' test_as_row(prr(mds_ts[[1]]))
#' @export
test_as_row <- function(
  df
){
  input_param_checker(df, "mdsstat_test")
  # Set eval period to length of time series as default behavior
  eval_period <- ifelse(is.null(df$params$eval_period), nrow(df$data$data),
                        df$params$eval_period)

  out <- data.frame(test_name=df$test_name,
                    analysis_of=df$analysis_of,
                    run_status=df$status,
                    run_msg=names(df$status),
                    ref_time_start=df$data$reference_time[1],
                    ref_time_end=df$data$reference_time[2],
                    eval_period=eval_period,
                    test_hyp=df$params$test_hyp,
                    test_params=I(list(df$params[!names(df$params) %in%
                                                   c("test_hyp",
                                                     "eval_period")])))
  if (all(is.na(df$result))){
    out <- cbind(out,
                 data.frame(signal=NA,
                            signal_threshold=I(list(NA)),
                            stat=I(list(NA)),
                            stat_lcl=I(list(NA)),
                            stat_ucl=I(list(NA)),
                            p_value=NA,
                            stat_addtl=I(list(NA))))
  } else{
    out <- cbind(out,
                 data.frame(signal=df$result$signal,
                            signal_threshold=I(list(df$result$signal_threshold)),
                            stat=I(list(df$result$statistic)),
                            stat_lcl=I(list(df$result$lcl)),
                            stat_ucl=I(list(df$result$ucl)),
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
#' @details Each algorithm may be named multiple times (to
#' allow running of multiple parameter settings). Do not specify the \code{df}
#' parameter.
#' @examples
#' x <- list(prr=list(),
#'   xbar=list(),
#'   xbar=list(ts_event=c(Rate="rate"), we_rule=2),
#'   poisson_rare=list(p_rate=0.3))
#' define_algos(x)
#' @export
define_algos <- function(
  algos
){
  input_param_checker(algos, "list")

  # Each list item must be an algorithm in mdsstat
  algolist <- c("poisson_rare", "prr", "shewhart", "xbar",
                "cusum", "ror",
                "sprt", "gps", "bcpnn",
                "ewma", "cp_mean",
                "lrt", "poisson_maxsprt")
  # algolist <- ls("package:mdsstat")[grepl("\\.mds_ts$", ls("package:mdsstat"))]
  # algolist <- gsub("\\.mds_ts$", "", algolist)
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
    tname <- names(algos)[i]
    tmeth <- c(tname, as.character(utils::methods(tname)))
    args <- unique(unlist(sapply(tmeth, function(x) names(formals(x)))))
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
#' dataset requirements for each algorithm specified (parameter \code{df}). An
#' \code{mds} times series object (class \code{mds_ts}) is a natural fit.
#' @param algos Input list of algorithms to run. Must be a list generated by
#' \code{\link{define_algos}}.
#' @param dataframe Logical on whether to output results as a
#' \code{mdsstat_tests} data frame. If \code{FALSE}, will output as a list of
#' \code{mdsstat_test} lists.
#'
#' Default: \code{TRUE}
#'
#' @param non_dpa What to do when input \code{data} is not prepared for
#' disproportionality analysis (DPA) data. Three values are accepted:
#' \code{"skip"}, \code{"warn"}, and \code{"stop"}. \code{"skip"} skips the
#' DPA test without warnings or errors. \code{"warn"} outputs a warning and
#' then skips the DPA test. \code{"stop"} stops the function call.
#'
#' Default: \code{"skip"}
#'
#' @param ... Further arguments for future work.
#' @return A \code{mdsstat_tests} data frame or list of \code{mdsstat_test}
#' lists with the results of the algorithm runs.
#' @examples
#' data <- mds_ts[[1]]
#' data$rate <- data$nA / data$exposure
#' x <- list(prr=list(),
#'   xbar=list(),
#'   xbar=list(ts_event=c(Rate="rate"), we_rule=2),
#'   poisson_rare=list(p_rate=0.3))
#' algos <- define_algos(x)
#' run_algos(data, algos)
#' run_algos(data, algos, FALSE)
#' @export
run_algos <- function(
  data,
  algos,
  dataframe=T,
  non_dpa="skip",
  ...
){
  UseMethod("run_algos", data)
}

#' @describeIn run_algos Run algorithms on a list of time series
#' @export
run_algos.list <- function(
  data,
  algos,
  dataframe=T,
  non_dpa="skip",
  ...
){
  dots <- list(...)
  if (dataframe){
    out <- data.frame()
  } else out <- list()
  for(i in 1:length(data)){
    this <- run_algos(data=data[[i]],
                      algos=algos,
                      dataframe=dataframe,
                      non_dpa=non_dpa,
                      ...)
    if (dataframe){
      out <- rbind(out, this)
    } else{
      out[[i]] <- this
    }
  }
  return(out)
}

#' @describeIn run_algos Run algorithms on a single time series
#' @export
run_algos.default <- function(
  data,
  algos,
  dataframe=T,
  non_dpa="skip",
  ...
){
  input_param_checker(data, "data.frame")
  input_param_checker(algos, "mdsstat_da")
  input_param_checker(dataframe, "logical")
  if (!non_dpa %in% c("skip", "stop", "warn")){
    stop(paste("non_dpa must have one of the following values:",
               "skip, stop, warn"))
  }
  # Define DPA algorithms currently in mdsstat
  dpaalgos <- c("prr", "ror", "gps", "bcpnn")

  if (dataframe){
    stats <- data.frame()
  } else stats <- list()
  for (i in 1:length(algos)){
    algo <- eval(parse(text=names(algos)[i]))
    # Special handling for DPA algorithms
    flag <- F
    if (names(algos)[i] %in% dpaalgos){
      if (!all(c("nA", "nB", "nC", "nD") %in% names(data))){
        if (non_dpa == "stop"){
          stop(paste("data is not in the required format for", names(algos)[i],
                     "analysis"))
        } else if (non_dpa == "warn"){
          warning(paste("data is not in the required format for", names(algos)[i],
                        "analysis. Skipping analysis."))
        }
        flag <- T
      }
    }
    # Run algorithm
    if (!flag){
      test <- do.call(algo, c(list(df=data), algos[[i]]))
      test$ts_id <- attributes(data)$analysis$id
      if (dataframe){
        ts_row <- test_as_row(test)
        ts_row$ts_id <- test$ts_id
        stats <- rbind(stats, ts_row)
      } else stats[[i]] <- test
    }
  }

  if (dataframe) class(stats) <- append(class(stats), "mdsstat_tests")

  # # Final check for invalid input data
  # if (is.data.frame(stats)){
  #   if (nrow(stats) == 0) stop("Input dataset is invalid format.")
  # } else if (is.list(stats)){
  #   if (length(stats[[1]]) == 0) stop("Input dataset is invalid format.")
  # }

  return(stats)
}

