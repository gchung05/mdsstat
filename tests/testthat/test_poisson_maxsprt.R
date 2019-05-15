context("Poisson MaxSPRT Algorithm")

# Reference example
data <- data.frame(time=c(1:100), event=as.integer(stats::rpois(100, 1)))
data$event[100] <- 6
a1 <- poisson_maxsprt(data)

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
  expect_equal(a1$test_name, "Poisson MaxSPRT")
  expect_equal(a1$analysis_of, NA)
  expect_true(a1$status)
  expect_true(all(names(a1$result) %in% c("statistic",
                                          "lcl",
                                          "ucl",
                                          "p",
                                          "signal",
                                          "signal_threshold",
                                          "sample_threshold",
                                          "c_t",
                                          "u_t")))
  expect_true(abs(a1$result$statistic) > 0)
  expect_true(is.na(a1$result$lcl))
  expect_true(is.na(a1$result$ucl))
  expect_true(is.na(a1$result$p))
  expect_is(a1$result$signal, "logical")
  expect_true(a1$result$signal_threshold > 0)
  expect_true(a1$result$sample_threshold > 0)
  expect_true(a1$result$c_t > 0)
  expect_true(a1$result$u_t > 0)
  expect_true(a1$result$c_t > a1$result$u_t)
  expect_true(all(names(a1$params) %in% c("test_hyp",
                                          "eval_period",
                                          "obs_period",
                                          "alpha",
                                          "M",
                                          "u_t")))
  expect_is(a1$params$test_hyp, "character")
  expect_equal(a1$params$eval_period, 100)
  expect_equal(a1$params$obs_period, 3)
  expect_equal(a1$params$alpha, 0.05)
  expect_equal(a1$params$M, 4)
  expect_true(is.na(a1$params$u_t))
  expect_true(all(names(a1$data) %in% c("reference_time",
                                          "data")))
  expect_equal(a1$data$reference_time, range(a1$data$data$time))
  expect_equal(a1$data$data, data)
})

# Parameter checks
# ----------------
a2d <- data
a2d$rate <- ifelse(is.na(a2d$event), 0, a2d$event)
a2d$exposure <- stats::rnorm(100, 50, 5)
a2d$rate <- a2d$rate / a2d$exposure
a2 <- poisson_maxsprt(a2d, ts_event=c("Rate"="rate"))
test_that("ts_event parameter functions as expected", {
  expect_equal(a2$data$reference_time, range(a2d$time))
  expect_equal(a2$data$data[[1]], a2d$time)
  expect_equal(a2$data$data$rate, a2d$rate)
  skip_on_cran()
  expect_equal(poisson_maxsprt(a2d, analysis_of="Testing")$analysis_of,
               "Testing")
})

test_that("eval_period parameter functions as expected", {
  expect_error(poisson_maxsprt(a2d, eval_period=0))
  expect_error(poisson_maxsprt(a2d, eval_period=5.5))
  expect_error(poisson_maxsprt(a2d, eval_period=101))
  expect_error(poisson_maxsprt(a2d, eval_period=100, obs_period=101))
})

test_that("obs_period parameter functions as expected", {
  expect_error(poisson_maxsprt(a2d, obs_period=0))
  expect_error(poisson_maxsprt(a2d, obs_period=5.5))
  expect_error(poisson_maxsprt(a2d, obs_period=100))
  expect_error(poisson_maxsprt(a2d, obs_period=101))
  expect_error(poisson_maxsprt(a2d, obs_period="now"))
  skip_on_cran()
  expect_is(poisson_maxsprt(a2d, obs_period=1), "mdsstat_test")
  expect_is(poisson_maxsprt(a2d, obs_period=9), "mdsstat_test")
})

a2i <- poisson_maxsprt(a2d, alpha=0.01)
test_that("alpha parameter functions as expected", {
  expect_equal(a2i$params$alpha, 0.01)
  expect_is(a2i, "mdsstat_test")
  expect_error(poisson_maxsprt(a2d, alpha=0))
  expect_error(poisson_maxsprt(a2d, alpha=1))
  expect_error(poisson_maxsprt(a2d, alpha=.51))
  expect_error(poisson_maxsprt(a2d, alpha=-.01))
  expect_error(poisson_maxsprt(a2d, alpha=1.01))
})

a2i <- poisson_maxsprt(a2d, M=1)
test_that("M parameter functions as expected", {
  expect_equal(a2i$params$M, 1)
  expect_is(a2i, "mdsstat_test")
  expect_error(poisson_maxsprt(a2d, M=0))
  expect_error(poisson_maxsprt(a2d, M=1.5))
  expect_error(poisson_maxsprt(a2d, M=-1))
  expect_error(poisson_maxsprt(a2d, M=-.01))
  expect_error(poisson_maxsprt(a2d, M=1.01))
})

a2i <- poisson_maxsprt(a2d, u_t=1)
test_that("u_t parameter functions as expected", {
  expect_equal(a2i$params$u_t, 1)
  expect_equal(a2i$result$u_t, 1)
  expect_is(a2i, "mdsstat_test")
  expect_error(poisson_maxsprt(a2d, u_t=0))
  expect_error(poisson_maxsprt(a2d, u_t=-1))
  expect_error(poisson_maxsprt(a2d, u_t=-.01))
  skip_on_cran()
  expect_equal(poisson_maxsprt(a2d, u_t=1.01)$params$u_t, 1.01)
  expect_equal(poisson_maxsprt(a2d, u_t=0.01)$params$u_t, 0.01)
})
