context("DPA Algorithms")

# Reference example
data <- data.frame(time=c(1:25), event=as.integer(rnorm(25, 100, 25)))
a1 <- shewhart(data)

# Basic
# -----

# Return behavior
test_that("function returns the correct class", {
  expect_is(a1, "list")
  expect_is(a1, "mdsstat_test")
})
test_that("function returns core mdsstat_test components", {
  expect_true(all(names(a1) %in% c("test_name",
                                   "analysis_of",
                                   "status",
                                   "result",
                                   "params",
                                   "data")))
})
test_that("outputs are as expected", {
  expect_equal(a1$test_name, "Shewhart x-bar Western Electric Rule 1")
  expect_equal(a1$analysis_of, NA)
  expect_true(a1$status, T)
  expect_true(all(names(a1$result) %in% c("statistic",
                                          "ll95",
                                          "ul95",
                                          "p",
                                          "signal",
                                          "signal_threshold",
                                          "mu",
                                          "sigma")))
  expect_true(abs(a1$result$statistic) > 0)
  expect_true(abs(a1$result$ll95) > 0)
  expect_true(abs(a1$result$ul95) > 0)
  expect_true(a1$result$ul95 > a1$result$ll95)
  expect_true(is.na(a1$result$p))
  expect_is(a1$result$signal, "logical")
  expect_true(a1$result$signal_threshold > 0)
  expect_true(all(names(a1$params) %in% c("test_hyp",
                                          "zero_rate",
                                          "we_rule")))
  expect_is(a1$params$test_hyp, "character")
  expect_equal(a1$params$zero_rate, 1/3)
  expect_equal(a1$params$we_rule, 1L)
  expect_true(all(names(a1$data) %in% c("reference_time",
                                          "data")))
  expect_equal(a1$data$reference_time, 25)
  expect_equal(a1$data$data, data)
})

# Reference example
data <- data.frame(time=c(1:8), event=c(rep(0, 6), rpois(2, 4)))
a1a <- shewhart(data)

test_that("test does not run on rare events", {
  expect_true(isFALSE(a1a$status))
})

# Parameter checks
# ----------------
a2d <- mds_ts[[3]]
a2 <- shewhart(a2d)
test_that("df parameter functions as expected", {
  expect_is(a2, "list")
  expect_is(a2, "mdsstat_test")
  expect_true(all(names(a2) %in% c("test_name",
                                   "analysis_of",
                                   "status",
                                   "result",
                                   "params",
                                   "data")))
  expect_equal(a2$test_name, "Shewhart x-bar Western Electric Rule 1")
  expect_match(a2$analysis_of, "Count of .+")
  expect_true(a2$status)
  expect_true(all(names(a2$result) %in% c("statistic",
                                          "ll95",
                                          "ul95",
                                          "p",
                                          "signal",
                                          "signal_threshold",
                                          "mu",
                                          "sigma")))
  expect_true(abs(a2$result$statistic) > 0)
  expect_true(abs(a2$result$ll95) > 0)
  expect_true(abs(a2$result$ul95) > 0)
  expect_true(a2$result$ul95 > a2$result$ll95)
  expect_true(is.na(a1$result$p))
  expect_is(a2$result$signal, "logical")
  expect_true(a2$result$signal_threshold > 0)
  expect_true(all(names(a2$params) %in% c("test_hyp",
                                          "zero_rate",
                                          "we_rule")))
  expect_is(a2$params$test_hyp, "character")
  expect_equal(a2$params$zero_rate, 1/3)
  expect_equal(a2$params$we_rule, 1L)
  expect_true(all(names(a2$data) %in% c("reference_time",
                                        "data")))
  expect_equal(a2$data$reference_time, max(a2d$time))
  expect_equal(a2$data$data[[1]], a2d[[1]])
  expect_equal(a2$data$data[[2]], ifelse(is.na(a2d$nA), 0, a2d$nA))
})

a2d <- mds_ts[[3]]
a2d$rate <- ifelse(is.na(a2d$nA), 0, a2d$nA)
a2d$rate <- a2d$rate / a2d$exposure
a2 <- shewhart(a2d, ts_event="rate")
test_that("ts_event parameter functions as expected", {
  expect_equal(a2$data$reference_time, max(a2d$time))
  expect_equal(a2$data$data[[1]], a2d$time)
  expect_equal(a2$data$data[[2]], a2d$rate)
})

a2 <- shewhart(a2d, analysis_of="Testing")
test_that("ts_event parameter functions as expected", {
  expect_equal(a2$analysis_of, "Testing")
})

a2 <- shewhart(a2d, eval_period=3L)
test_that("eval_period parameter functions as expected", {
  expect_equal(names(a2$status), ">3 time periods required")
  expect_equal(nrow(a2$data$data), 3L)
  expect_error(eval_period(a2d, eval_period=2))
  expect_error(eval_period(a2d, eval_period=0))
  expect_error(eval_period(a2d, eval_period=nrow(a2d) + 1))
})

test_that("zero_rate parameter functions as expected", {
  expect_is(shewhart(a2d, zero_rate=0), "mdsstat_test")
  expect_is(shewhart(a2d, zero_rate=1), "mdsstat_test")
  expect_error(shewhart(a2d, zero_rate=1.1))
  expect_error(shewhart(a2d, zero_rate=-1))
})

a22L <- shewhart(a2d, we_rule=2L)
a23L <- shewhart(a2d, we_rule=3L)
a24L <- shewhart(a2d, we_rule=4L)
test_that("we_rule parameter functions as expected", {
  expect_is(a22L, "mdsstat_test")
  expect_is(a23L, "mdsstat_test")
  expect_is(a24L, "mdsstat_test")
  expect_error(shewhart(a2d, we_rule=1))
  expect_error(shewhart(a2d, we_rule=5L))
  expect_error(shewhart(a2d, we_rule=5))
  expect_error(shewhart(a2d, we_rule=0L))
})
test_that("we_rule 2 functions as expected", {
  expect_equal(a22L$test_name, "Shewhart x-bar Western Electric Rule 2")
  expect_equal(length(a22L$result$statistic), 3)
  expect_equal(length(a22L$result$ll95), 3)
  expect_equal(length(a22L$result$ul95), 3)
  expect_equal(length(a22L$result$signal_threshold), 3)
})
test_that("we_rule 3 functions as expected", {
  expect_equal(a23L$test_name, "Shewhart x-bar Western Electric Rule 3")
  expect_equal(length(a23L$result$statistic), 5)
  expect_equal(length(a23L$result$ll95), 5)
  expect_equal(length(a23L$result$ul95), 5)
  expect_equal(length(a23L$result$signal_threshold), 5)
})
test_that("we_rule 4 functions as expected", {
  expect_equal(a24L$test_name, "Shewhart x-bar Western Electric Rule 4")
  expect_equal(length(a24L$result$statistic), 9)
  expect_equal(length(a24L$result$ll95), 9)
  expect_equal(length(a24L$result$ul95), 9)
  expect_equal(length(a24L$result$signal_threshold), 9)
})
