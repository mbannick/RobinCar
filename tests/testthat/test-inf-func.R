library(survival)
library(tidyverse)

files <- Sys.glob("~/Documents/GitHub/RobinCar/R/*.R")
map(files, source)

DATA <- ovarian %>%
  rename(tte = futime, obs = fustat) %>%
  arrange(tte)

test_that("Logrank influence function error checking 6", {
  DATA$id <- DATA$ecog.ps
  expect_error(robincar_logrank(
    adj_method = "CL",
    df = DATA,
    treat_col = "rx",
    p_trt = 0.5,
    ref_arm = 1,
    response_col = "tte",
    event_col = "obs",
    return_influence = T,
    id_col = "idx"
  ), "Column idx not found")
})

test_that("Logrank influence function error checking 5", {
  DATA$id <- DATA$ecog.ps
  expect_no_error(robincar_logrank(
    adj_method = "CL",
    df = DATA,
    treat_col = "rx",
    p_trt = 0.5,
    covariate_col = "id",
    ref_arm = 1,
    response_col = "tte",
    event_col = "obs"
  ))
})

test_that("Logrank influence function error checking 4", {
  DATA$idx <- 1:nrow(DATA)
  DATA$id <- DATA$ecog.ps
  expect_error(robincar_logrank(
    adj_method = "CL",
    df = DATA,
    treat_col = "rx",
    p_trt = 0.5,
    covariate_col = "id",
    ref_arm = 1,
    response_col = "tte",
    event_col = "obs",
    return_influence = T,
    id_col = "idx"
  ), "if return_influence is TRUE, no covariate can be called 'id'.")
})

test_that("Logrank influence function error checking 3", {
  DATA$idx <- 1:nrow(DATA)
  expect_error(robincar_logrank(
    adj_method = "CL",
    df = DATA,
    treat_col = "rx",
    p_trt = 0.5,
    ref_arm = 1,
    response_col = "tte",
    event_col = "obs",
    return_influence = "cheese",
    id_col = "idx"
  ), "return_influence must be either TRUE or FALSE")
})

test_that("Logrank influence function error checking 2", {
  DATA$idx <- 5
  expect_error(robincar_logrank(
    adj_method = "CL",
    df = DATA,
    treat_col = "rx",
    p_trt = 0.5,
    ref_arm = 1,
    response_col = "tte",
    event_col = "obs",
    return_influence = T,
    id_col = "idx"
  ), "Id column must not have duplicated values")
})

test_that("Logrank influence function error checking 1", {
  expect_error(robincar_logrank(
    adj_method = "CL",
    df = DATA,
    treat_col = "rx",
    p_trt = 0.5,
    ref_arm = 1,
    response_col = "tte",
    event_col = "obs",
    return_influence = T
  ), "id_col must not be NULL if return_influence is TRUE")
})

test_that("Logrank influence function", {
  DATA$idx <- 1:nrow(DATA)
  RC1 <- robincar_logrank(
    adj_method = "CL",
    df = DATA,
    treat_col = "rx",
    p_trt = 0.5,
    ref_arm = 1,
    response_col = "tte",
    event_col = "obs",
    return_influence = T,
    id_col = "idx"
  )
  expect_equal(sum(RC1$result$inf_func$inf_func), RC1$result$strata_sum$U_SL_z[1])
})