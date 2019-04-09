context("CP-Mean Algorithm")

# Reference example
data <- data.frame(time=seq(1:49),
                   event=c(rnorm(12, 35, 3), #13 ideal changepoint
                           rnorm(4, 45, 5), #17
                           rnorm(8, 53, 8), #25
                           rnorm(15, 30, 3), #40
                           rnorm(5, 40, 5), #45
                           rnorm(5, 80, 2))) #46
a1 <- cp_mean(data)

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
  expect_equal(a1$test_name, "Mean-shift changepoint")
  expect_equal(a1$analysis_of, NA)
  expect_true(a1$status)
  expect_true(all(names(a1$result) %in% c("statistic",
                                          "lcl",
                                          "ucl",
                                          "p",
                                          "signal",
                                          "signal_threshold")))
  expect_true(all(a1$result$statistic > 0))
  expect_true(is.na(a1$result$lcl))
  expect_true(is.na(a1$result$ucl))
  expect_true(all(a1$result$p >= 0))
  expect_equal(length(a1$result$p), length(a1$result$statistic))
  expect_is(a1$result$signal, "logical")
  expect_equal(as.numeric(a1$result$signal_threshold), a1$params$alpha)
  expect_true(all(names(a1$params) %in% c("test_hyp",
                                          "eval_period",
                                          "zero_rate",
                                          "alpha",
                                          "cp_max",
                                          "min_seglen",
                                          "epochs",
                                          "bootstrap_iter",
                                          "replace")))
  expect_is(a1$params$test_hyp, "character")
  expect_equal(a1$params$zero_rate, 1/3)
  expect_equal(a1$params$alpha, 0.05)
  expect_equal(a1$params$cp_max, 100)
  expect_equal(a1$params$min_seglen, 6)
  expect_equal(a1$params$epochs, 4)
  expect_equal(a1$params$bootstrap_iter, 1000)
  expect_equal(a1$params$replace, TRUE)
  expect_true(all(names(a1$data) %in% c("reference_time",
                                          "data")))
  expect_equal(a1$data$reference_time, range(a1$data$data$time))
  expect_equal(a1$data$data, data)
})

# Reference example
data2 <- data.frame(time=c(1:8), event=c(rep(0, 6), rpois(2, 4)))
a1a <- cp_mean(data2)

test_that("test does not run on rare events", {
  expect_true(!a1a$status)
})

# Parameter checks
# ----------------
a2d <- data

a22 <- cp_mean(a2d, alpha=.1)
a245 <- cp_mean(a2d, alpha=.2)
test_that("alpha parameter functions as expected", {
  expect_is(a22, "mdsstat_test")
  expect_equal(a22$params$alpha, .1)
  expect_equal(a245$params$alpha, .2)
  expect_error(cp_mean(a2d, alpha=0))
  expect_error(cp_mean(a2d, alpha=1))
  expect_error(cp_mean(a2d, alpha=-1))
})

a21 <- cp_mean(a2d, cp_max=10)
a203 <- cp_mean(a2d, cp_max=1)
test_that("cp_max parameter functions as expected", {
  expect_is(a21, "mdsstat_test")
  expect_equal(a21$params$cp_max, 10)
  expect_equal(a203$params$cp_max, 1)
  expect_equal(length(a203$result$statistic), 1)
  expect_error(cp_mean(a2d, cp_max=0))
  expect_error(cp_mean(a2d, cp_max=2.5))
  expect_error(cp_mean(a2d, cp_max=-1))
})

a21 <- cp_mean(a2d, min_seglen=10)
a203 <- cp_mean(a2d, min_seglen=12)
test_that("min_seglen parameter functions as expected", {
  expect_is(a21, "mdsstat_test")
  expect_equal(a21$params$min_seglen, 10)
  expect_equal(a203$params$min_seglen, 12)
  expect_error(cp_mean(a2d, min_seglen=0))
  expect_error(cp_mean(a2d, min_seglen=3))
  expect_error(cp_mean(a2d, min_seglen=3.5))
  expect_error(cp_mean(a2d, min_seglen=-1))
})

a21 <- cp_mean(a2d, epochs=10)
a203 <- cp_mean(a2d, epochs=12)
test_that("epochs parameter functions as expected", {
  expect_is(a21, "mdsstat_test")
  expect_equal(a21$params$epochs, 10)
  expect_equal(a203$params$epochs, 12)
  expect_error(cp_mean(a2d, epochs=0))
  expect_error(cp_mean(a2d, epochs=0.5))
  expect_error(cp_mean(a2d, epochs=-1))
})

a21 <- cp_mean(a2d, bootstrap_iter=1001)
a203 <- cp_mean(a2d, bootstrap_iter=1002)
test_that("bootstrap_iter parameter functions as expected", {
  expect_is(a21, "mdsstat_test")
  expect_equal(a21$params$bootstrap_iter, 1001)
  expect_equal(a203$params$bootstrap_iter, 1002)
  expect_error(cp_mean(a2d, bootstrap_iter=500))
  expect_error(cp_mean(a2d, bootstrap_iter=1000.5))
  expect_error(cp_mean(a2d, bootstrap_iter=-1))
})

a21 <- cp_mean(a2d, replace=F)
test_that("replace parameter functions as expected", {
  expect_is(a21, "mdsstat_test")
  expect_equal(a21$params$replace, F)
  expect_error(cp_mean(a2d, replace="no"))
  expect_error(cp_mean(a2d, replace=1))
  expect_error(cp_mean(a2d, replace="true"))
})

a21 <- cp_mean(a2d, zero_rate=1/2)
a203 <- cp_mean(a2d, zero_rate=0)
test_that("zero_rate parameter functions as expected", {
  expect_is(a21, "mdsstat_test")
  expect_equal(a21$params$zero_rate, 1/2)
  expect_equal(a203$params$zero_rate, 0)
  expect_error(cp_mean(a2d, zero_rate=1.1))
  expect_error(cp_mean(a2d, zero_rate=-1))
  expect_error(cp_mean(a2d, zero_rate=2))
})
