library(survival)
library(dplyr)

test_that("No X no Z case under SR, CL = logrank test", {

  set.seed(0)
  n <- 100
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

  out <- robincar_logrank(df=data.simu0,
                          treat_col="I1",
                          p_trt=0.5,
                          ref_arm=0,
                          response_col="t",
                          event_col="delta",
                          car_scheme=c("simple"),
                          adj_method=c("CL")
  )

  surv_out <- coxph(Surv(t, delta) ~ I1, data = data.simu0)
  surv_out.score <-
    sqrt(surv_out$score) * ifelse(surv_out$coefficients > 0, 1, -1)

  expect_equal(out$result$statistic, unname(surv_out.score))

})

test_that("No X yes Z case under CAR, CSL = stratified logrank test",{

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
    tidyr::pivot_longer(
      cols=starts_with("car_strata"),
      names_prefix="car_strata",
      names_to="strt"
    ) %>%
    dplyr::filter(value == 1) %>%
    dplyr::select(-value) %>%
    dplyr::mutate(strt = forcats::as_factor(strt)) %>%
    dplyr::select(t, strt) %>%
    dplyr::left_join(data.simu0, .)

  out1 <- robincar_logrank(
    df=data.simu,
    treat_col="I1",
    p_trt=0.5,
    ref_arm=0,
    response_col="t",
    event_col="delta",
    car_strata_cols=c("strt"),
    car_scheme=c("permuted-block"),
    adj_method=c("CSL")
  )

  surv_out1 <- coxph(Surv(t, delta) ~ I1 + strata(strt), data = data.simu)
  surv_out.score1 <- sqrt(surv_out1$score) *
    ifelse(surv_out1$coefficients > 0, 1, -1)

  expect_equal(out1$result$statistic, unname(surv_out.score1))

})

test_that("X yes Z yes, case1, CL",{

  set.seed(0)
  n <- 100
  data.simu0 <- data_gen(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case1"
  )

  data.simu <-
    data.simu0 %>%
    tidyr::pivot_longer(
      cols=starts_with("car_strata"),
      names_prefix="car_strata",
      names_to="strt"
    ) %>%
    dplyr::filter(value == 1) %>%
    dplyr::select(-value) %>%
    dplyr::mutate(strt = forcats::as_factor(strt)) %>%
    dplyr::select(t, strt) %>%
    dplyr::left_join(data.simu0, .)

  test1 <- robincar_logrank(
    df=data.simu,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols=c("strt"),
    covariate_cols=c("model_Z1", "model_Z21", "model_Z22"),
    car_scheme="permuted-block",
    adj_method=c("CL"),
    ref_arm=0
  )
  test1_ty <- covariate_adjusted_logrank(data.simu, p_trt = 0.5)
  expect_equal(test1$result$U, test1_ty$U_CL)
  expect_equal(test1$result$se, as.numeric(test1_ty$se))

  data.simu01 <- data_gen2(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case1"
  )

  data.simu1 <- data.simu01 %>%
    tidyr::pivot_longer(
      cols=starts_with("car_strata"),
      names_prefix="car_strata",
      names_to="strt"
    ) %>%
    dplyr::filter(value == 1) %>%
    dplyr::select(-value) %>%
    dplyr::mutate(strt=forcats::as_factor(strt)) %>%
    dplyr::select(t, strt) %>%
    dplyr::left_join(data.simu01, .)

  test12 <- robincar_logrank(
    df=data.simu1,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols=c("strt"),
    covariate_cols=c("model_w3"),
    car_scheme="permuted-block",
    adj_method=c("CL"),
    ref_arm=0
  )
  test12_ty <- covariate_adjusted_logrank(data.simu1, p_trt = 0.5)
  expect_equal(test12$result$U, test12_ty$U_CL)
  expect_equal(test12$result$se, as.numeric(test12_ty$se))

})

test_that("X (continuous) yes Z no, case1, CL",{

  set.seed(0)
  n <- 100

  data.simu01 <- data_gen2(
    n=n,
    theta=0,
    randomization="SR",
    p_trt=0.5,
    case="case1"
  )

  # Add another continuous covariate because Ting's code requires
  data.simu1 <- data.simu01 %>%
    dplyr::mutate(strt=0)

  test12 <- robincar_logrank(
    df=data.simu1,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    covariate_cols=c("model_w3"),
    car_scheme="simple",
    adj_method=c("CL"),
    ref_arm=0
  )

  # Using the function with no Z, modified version of Ting's original code
  test12_ty <- covariate_adjusted_logrank(data.simu1, p_trt = 0.5, Z=FALSE)
  expect_equal(test12$result$U, test12_ty$U_CL)
  expect_equal(test12$result$se, as.numeric(test12_ty$se))

})

test_that("X yes Z yes, case1, CSL",{

  set.seed(100)
  n <- 100

  data.simu01 <- data_gen2(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case1"
  )

  data.simu1 <-
    data.simu01 %>% tidyr::pivot_longer(
      cols=starts_with("car_strata"),
      names_prefix="car_strata",
      names_to="strt"
    ) %>%
    dplyr::filter(value == 1) %>% select(-value) %>%
    dplyr::mutate(strt=forcats::as_factor(strt)) %>% select(t, strt) %>%
    dplyr::left_join(data.simu01, .)

  test12 <- robincar_logrank(
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
  test12_ty <- covariate_adjusted_stratified_logrank(data.simu1, p_trt=0.5)
  expect_equal(test12$result$U, test12_ty$U_CSL)
  expect_equal(test12$result$se, as.numeric(test12_ty$se))

})

