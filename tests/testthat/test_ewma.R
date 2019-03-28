context("EWMA Algorithm")

# Reference example
data <- data.frame(time=c(1:25), event=as.integer(stats::rnorm(25, 100, 25)))
a1 <- ewma(data)

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
  expect_equal(a1$test_name, "EWMA")
  expect_equal(a1$analysis_of, NA)
  expect_true(a1$status)
  expect_true(all(names(a1$result) %in% c("statistic",
                                          "lcl",
                                          "ucl",
                                          "p",
                                          "signal",
                                          "signal_threshold",
                                          "mu",
                                          "sigma",
                                          "ewma_sigma")))
  expect_true(all(a1$result$statistic >= 0))
  expect_true(all(abs(a1$result$lcl) >= 0))
  expect_true(all(a1$result$ucl >= 0))
  expect_equal(length(a1$result$statistic), nrow(data))
  expect_equal(length(a1$result$lcl), nrow(data))
  expect_equal(length(a1$result$ucl), nrow(data))
  expect_true(is.na(a1$result$p))
  expect_is(a1$result$signal, "logical")
  expect_true(a1$result$signal_threshold > 0)
  expect_true(all(names(a1$params) %in% c("test_hyp",
                                          "eval_period",
                                          "zero_rate",
                                          "delta",
                                          "lambda")))
  expect_is(a1$params$test_hyp, "character")
  expect_equal(a1$params$zero_rate, 1/3)
  expect_equal(a1$params$delta, 3)
  expect_is(a1$params$lambda, "numeric")
  expect_equal(a1$params$lambda, 0.85)
  expect_true(all(names(a1$data) %in% c("reference_time",
                                          "data")))
  expect_equal(a1$data$reference_time, range(a1$data$data$time))
  expect_equal(a1$data$data, data)
})

# Reference example
data <- data.frame(time=c(1:8), event=c(rep(0, 6), rpois(2, 4)))
a1a <- ewma(data)

test_that("test does not run on rare events", {
  expect_true(!a1a$status)
})

# Parameter checks
# ----------------
data <- data.frame(time=c(1:25), event=as.integer(stats::rnorm(25, 100, 25)))
a2d <- data
a2 <- ewma(a2d)
test_that("df parameter functions as expected", {
  expect_is(a2, "list")
  expect_is(a2, "mdsstat_test")
  expect_true(all(names(a2) %in% c("test_name",
                                   "analysis_of",
                                   "status",
                                   "result",
                                   "params",
                                   "data")))
  expect_equal(a2$test_name, "EWMA")
  expect_equal(a2$analysis_of, NA)
  expect_true(a2$status)
  expect_true(all(names(a2$result) %in% c("statistic",
                                          "lcl",
                                          "ucl",
                                          "p",
                                          "signal",
                                          "signal_threshold",
                                          "mu",
                                          "sigma",
                                          "ewma_sigma")))
  expect_true(all(a2$result$statistic >= 0))
  expect_true(all(abs(a2$result$lcl) >= 0))
  expect_true(all(a2$result$ucl >= 0))
  expect_equal(length(a2$result$statistic), nrow(a2d))
  expect_equal(length(a2$result$lcl), nrow(a2d))
  expect_equal(length(a2$result$ucl), nrow(a2d))
  expect_true(is.na(a2$result$p))
  expect_is(a2$result$signal, "logical")
  expect_true(a2$result$signal_threshold >= 0)
  expect_true(all(names(a2$params) %in% c("test_hyp",
                                          "eval_period",
                                          "zero_rate",
                                          "delta",
                                          "lambda")))
  expect_is(a2$params$test_hyp, "character")
  expect_equal(a2$params$zero_rate, 1/3)
  expect_equal(a2$params$delta, 3)
  expect_is(a2$params$lambda, "numeric")
  expect_equal(a2$params$lambda, 0.85)
  expect_true(all(names(a2$data) %in% c("reference_time",
                                        "data")))
  expect_equal(a2$data$reference_time, range(a2d$time))
  expect_equal(a2$data$data[[1]], a2d[[1]])
  expect_equal(a2$data$data[[2]], ifelse(is.na(a2d$event), 0, a2d$event))
})

a2d <- data
a2d$rate <- ifelse(is.na(a2d$event), 0, a2d$event)
a2d$exposure <- stats::rnorm(25, 50, 5)
a2d$rate <- a2d$rate / a2d$exposure
a2 <- ewma(a2d, ts_event=c("Rate"="rate"))
test_that("ts_event parameter functions as expected", {
  expect_equal(a2$data$reference_time, range(a2d$time))
  expect_equal(a2$data$data[[1]], a2d$time)
  expect_equal(a2$data$data$rate, a2d$rate)
})

a2 <- ewma(a2d, analysis_of="Testing")
test_that("ts_event parameter functions as expected", {
  expect_equal(a2$analysis_of, "Testing")
})

a2 <- ewma(a2d, eval_period=3)
test_that("eval_period parameter functions as expected", {
  expect_equal(names(a2$status), ">3 time periods required")
  expect_equal(nrow(a2$data$data), 3)
  expect_equal(names(ewma(data, eval_period=2)$status),
               ">3 time periods required")
  expect_error(ewma(a2d, eval_period=0))
  expect_error(ewma(a2d, eval_period=nrow(a2d) + 1))
})

test_that("zero_rate parameter functions as expected", {
  expect_is(ewma(a2d, zero_rate=0), "mdsstat_test")
  expect_is(ewma(a2d, zero_rate=1), "mdsstat_test")
  expect_error(ewma(a2d, zero_rate=1.1))
  expect_error(ewma(a2d, zero_rate=-1))
})


a22 <- ewma(a2d, delta=2)
a245 <- ewma(a2d, delta=4.5)
test_that("delta parameter functions as expected", {
  expect_is(a22, "mdsstat_test")
  expect_equal(a22$params$delta, 2)
  expect_equal(a245$params$delta, 4.5)
  expect_error(ewma(a2d, delta=0))
  expect_error(ewma(a2d, delta=-1))
})

a21 <- ewma(a2d, lambda=1)
a203 <- ewma(a2d, lambda=0.3)
test_that("lambda parameter functions as expected", {
  expect_is(a21, "mdsstat_test")
  expect_equal(a21$params$lambda, 1)
  expect_equal(a203$params$lambda, 0.3)
  expect_error(ewma(a2d, lambda=2))
  expect_error(ewma(a2d, lambda=0))
  expect_error(ewma(a2d, lambda=-1))
})

