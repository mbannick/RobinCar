library(survival)
library(dplyr)

test_that("No X no Z case under SR, CL = logrank test",{

  set.seed(0)

  n <- 1000
  data.simu0 <- data_gen(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case2"
  )

  data.simu0 <- data.simu0 %>% mutate(
    strata1=sample(letters[1:3],n,replace=TRUE),
    strata2=sample(LETTERS[4:5],n,replace=TRUE)
  )

  out <- robincar_covhr(df=data.simu0,
                        treat_col="I1",
                        p_trt=0.5,
                        ref_arm=0,
                        response_col="t",
                        event_col="delta",
                        car_strata_cols=NULL,
                        covariate_cols=NULL,
                        car_scheme=c("simple"),
                        adj_method=c("CL")
  )

  surv_out <- coxph(Surv(t, delta) ~ I1, data=data.simu0)
  surv_out.score <-
    sqrt(surv_out$score) * ifelse(surv_out$coefficients > 0, 1, -1)

  expect_equal(out$result$theta_L, out$result$theta_CL)
  # Need a tolerance here because possibly a slightly different optimization
  # algorithm to find the hazard ratio than coxph
  expect_equal(out$result$theta_L, unname(surv_out$coefficients), tolerance=1e-5)
  expect_equal(out$result$se_theta_L, as.numeric(sqrt(surv_out$var)))

})



test_that("X yes Z yes, case1, CL", {

  set.seed(0)
  n <- 100
  data.simu0 <- data_gen(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case1"
  )

  data.simu <- data.simu0 %>%
    tidyr::pivot_longer(cols=starts_with("car_strata"),
                        names_prefix = "car_strata",
                        names_to = "strt") %>%
    filter(value==1) %>% select(-value) %>%
    mutate(strt=forcats::as_factor(strt)) %>%
    select(t,strt) %>%
    left_join(data.simu0,.)

  test1 <- robincar_covhr(
    df=data.simu,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols="strt",
    covariate_cols=c("model_Z1", "model_Z21", "model_Z22"),
    car_scheme="permuted-block",
    adj_method=c("CL"),
    ref_arm=0
  )
  test1_ty <- covariate_adjusted_logrank_estimation(data.simu, p_trt=0.5)
  expect_equal(test1$result$theta_L,test1_ty$theta_L)
  expect_equal(test1$result$theta_CL,test1_ty$theta_CL)
  expect_equal(test1$result$se_theta_L,test1_ty$se.theta_L)
  expect_equal(test1$result$se_theta_CL, as.numeric(test1_ty$se.theta_CL))

  data.simu01 <- data_gen2(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case1"
  )

  data.simu1 <- data.simu01 %>%
    tidyr::pivot_longer(cols=starts_with("car_strata"),
                        names_prefix = "car_strata",
                        names_to = "strt") %>%
    filter(value==1) %>% select(-value) %>%
    mutate(strt=forcats::as_factor(strt)) %>%
    select(t,strt) %>%
    left_join(data.simu01,.)

  test12 <- robincar_covhr(
    df=data.simu1,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols="strt",
    covariate_cols=c("model_w3"),
    car_scheme="permuted-block",
    adj_method=c("CL"),
    ref_arm=0
  )

  test12_ty <- covariate_adjusted_logrank_estimation(data.simu1, p_trt=0.5)
  expect_equal(test12$result$theta_L,test12_ty$theta_L)
  expect_equal(test12$result$theta_CL,test12_ty$theta_CL)
  expect_equal(test12$result$se_theta_L,test12_ty$se.theta_L)
  expect_equal(test12$result$se_theta_CL,as.numeric(test12_ty$se.theta_CL))

})


test_that("X yes Z yes, case1, CSL",{

  n <- 100

  data.simu01 <- data_gen2(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case1"
  )

  data.simu1 <- data.simu01 %>%
    tidyr::pivot_longer(cols=starts_with("car_strata"),
                        names_prefix = "car_strata",
                        names_to = "strt") %>%
    filter(value==1) %>% select(-value) %>%
    mutate(strt=forcats::as_factor(strt)) %>%
    select(t,strt) %>%
    left_join(data.simu01,.)

  test12 <- robincar_covhr(
    df=data.simu1,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols="strt",
    covariate_cols=c("model_w3"),
    car_scheme="permuted-block",
    adj_method=c("CSL"),
    ref_arm=0
  )
  test12_ty <- covariate_adjusted_stratified_logrank_estimation(data.simu1,p_trt=0.5)
  expect_equal(test12$result$theta_L,test12_ty$theta_SL)
  expect_equal(test12$result$theta_CL,test12_ty$theta_CSL)
  expect_equal(test12$result$se_theta_L,test12_ty$se.theta_SL)
  expect_equal(test12$result$se_theta_CL,as.numeric(test12_ty$se.theta_CSL))

})
