library(dplyr)

test_that("test finite sample correction", {
  
  set.seed(10)
  n <- 100
  
  # Simulate data
  data.simu <- data_gen2(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case1"
  )
  
  data.simu$strata_2 <- factor(rbinom(100, 1, 0.5))
  data.simu$strata_3 <- factor(rbinom(100, 1, 0.5))
  
  cl_1_F <- robincar_logrank(
    df=data.simu,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols=c("strata_2", "strata_3"),
    covariate_cols=c("model_w3"),
    car_scheme="permuted-block",
    adj_method="CL",
    ref_arm=0,
    finite_sample = FALSE
  )
  
  cl_1_T <- robincar_logrank(
    df=data.simu,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols=c("strata_2", "strata_3"),
    covariate_cols=c("model_w3"),
    car_scheme="permuted-block",
    adj_method="CL",
    ref_arm=0,
    finite_sample = TRUE
  )

  expect_equal(round(cl_1_F$result$se, 6), 0.045337)
  expect_equal(round(cl_1_T$result$se, 6), 0.047828)
  
  cl_2_F <- robincar_logrank(
    df=data.simu,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols=NULL,
    covariate_cols=c("model_w3", "strata_2", "strata_3"),
    car_scheme="simple",
    adj_method="CL",
    ref_arm=0,
    finite_sample = FALSE
  )
  
  cl_2_T <- robincar_logrank(
    df=data.simu,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols=NULL,
    covariate_cols=c("model_w3", "strata_2", "strata_3"),
    car_scheme="simple",
    adj_method="CL",
    ref_arm=0,
    finite_sample = TRUE
  )
  
  expect_equal(round(cl_2_F$result$se, 6), 0.045343)
  expect_equal(round(cl_2_T$result$se, 6), 0.047302)
  
  
  csl_1_F <- robincar_logrank(
    df=data.simu,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols=c("strata_2", "strata_3"),
    covariate_cols=c("model_w3"),
    car_scheme="permuted-block",
    adj_method="CSL",
    ref_arm=0,
    finite_sample = FALSE
  )
  
  csl_1_T <- robincar_logrank(
    df=data.simu,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols=c("strata_2", "strata_3"),
    covariate_cols=c("model_w3"),
    car_scheme="permuted-block",
    adj_method="CSL",
    ref_arm=0,
    finite_sample = TRUE
  )
  
  expect_equal(round(csl_1_F$result$se, 6), 0.043468)
  expect_equal(round(csl_1_T$result$se, 6), 0.044374)
  
})
