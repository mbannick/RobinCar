library(fastDummies)
library(dplyr)

test_that("sr", {
  set.seed(0)
  a <- car_sr(10000, p_trt=0.5)
  expect_equal(mean(a), 0.5, tolerance=1e-2)
})

test_that("pb", {
  set.seed(0)
  x <- runif(10000)
  z <- cut(x, breaks=c(0, 0.25, 0.5, 0.75, 1.0))
  z <- dummy_cols(z) %>%
       mutate(across(where(is.numeric), as.factor))

  a <- car_pb(z[, 2:5], c(0, 1, 2), trt_alc=c(1/4, 1/4, 2/4), blocksize=4L)
  expect_equal(mean(a == 0), 0.25, tolerance=1e-2)
  expect_equal(mean(a == 1), 0.25, tolerance=1e-2)
  expect_equal(mean(a == 2), 0.50, tolerance=1e-2)
})

test_that("ps", {
  set.seed(0)
  x <- runif(500)
  z <- cut(x, breaks=c(0, 0.25, 0.5, 0.75, 1.0))
  z <- dummy_cols(z)

  a <- car_ps(
    z=z[, 2:5],
    treat=c(0, 1, 2),
    ratio=c(1, 1, 1),
    imb_measure="Range"
  )
  expect_equal(mean(a==0), 0.33, tolerance=1e-1)
  expect_equal(mean(a==1), 0.33, tolerance=1e-1)
  expect_equal(mean(a==2), 0.33, tolerance=1e-1)
})
