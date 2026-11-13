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
    variance_correction="none"
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
    variance_correction="normal_approximation"
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
    variance_correction="none"
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
    variance_correction="normal_approximation"
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
    variance_correction="none"
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
    variance_correction="normal_approximation"
  )
  
  expect_equal(round(csl_1_F$result$se, 6), 0.043468)
  expect_equal(round(csl_1_T$result$se, 6), 0.044374)
  
})


test_that("hc3 correction inflates se relative to no correction (log-rank)", {

  set.seed(10)
  n <- 100

  data.simu <- data_gen2(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case1"
  )
  data.simu$strata_2 <- factor(rbinom(100, 1, 0.5))

  cl_none <- robincar_logrank(
    df=data.simu,
    treat_col="I1", response_col="t", event_col="delta",
    car_strata_cols=NULL,
    covariate_cols=c("model_w3", "strata_2"),
    car_scheme="simple",
    adj_method="CL", ref_arm=0,
    variance_correction="none"
  )

  cl_hc3 <- robincar_logrank(
    df=data.simu,
    treat_col="I1", response_col="t", event_col="delta",
    car_strata_cols=NULL,
    covariate_cols=c("model_w3", "strata_2"),
    car_scheme="simple",
    adj_method="CL", ref_arm=0,
    variance_correction="hc3"
  )

  # HC3 must inflate the SE
  expect_gt(cl_hc3$result$se, cl_none$result$se)

  # statistic must be recomputed consistently with the new SE
  expect_equal(cl_hc3$result$statistic,  cl_hc3$result$U  / cl_hc3$result$se)

})


test_that("hc3 correction inflates se_theta_CL relative to no correction (covhr)", {

  set.seed(10)
  n <- 100

  data.simu <- data_gen2(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case1"
  )
  data.simu$strata_2 <- factor(rbinom(100, 1, 0.5))

  hr_none <- robincar_covhr(
    df=data.simu,
    treat_col="I1", response_col="t", event_col="delta",
    covariate_cols=c("model_w3", "strata_2"),
    car_scheme="simple",
    adj_method="CL", ref_arm=0,
    variance_correction="none"
  )

  hr_hc3 <- robincar_covhr(
    df=data.simu,
    treat_col="I1", response_col="t", event_col="delta",
    covariate_cols=c("model_w3", "strata_2"),
    car_scheme="simple",
    adj_method="CL", ref_arm=0,
    variance_correction="hc3"
  )

  # HC3 must inflate the SE
  expect_gt(hr_hc3$result$se_theta_CL, hr_none$result$se_theta_CL)

  # theta_CL itself is unaffected by the variance correction
  expect_equal(hr_hc3$result$theta_CL, hr_none$result$theta_CL)

})


test_that("hc3 with no covariates leaves se unchanged", {

  set.seed(10)
  n <- 100

  data.simu <- data_gen2(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case1"
  )

  lr_none <- robincar_logrank(
    df=data.simu,
    treat_col="I1", response_col="t", event_col="delta",
    covariate_cols=NULL,
    car_scheme="simple",
    adj_method="CL", ref_arm=0,
    variance_correction="none"
  )

  lr_hc3 <- robincar_logrank(
    df=data.simu,
    treat_col="I1", response_col="t", event_col="delta",
    covariate_cols=NULL,
    car_scheme="simple",
    adj_method="CL", ref_arm=0,
    variance_correction="hc3"
  )

  expect_equal(lr_hc3$result$se, lr_none$result$se)

})
