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
