library(dplyr)
library(survival)

test_that("X yes Z yes, permuted-block, coxscore",{

  set.seed(10)
  n <- 100
  data.simu0 <- data_gen(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case2"
  )

  data.simu <-
    data.simu0 %>% tidyr::pivot_longer(
      cols=starts_with("car_strata"),
      names_prefix="car_strata",
      names_to="strt"
    ) %>%
    dplyr::filter(value == 1) %>%
    dplyr::select(-value) %>%
    dplyr::mutate(strt=forcats::as_factor(strt)) %>%
    dplyr::select(t, strt) %>%
    dplyr::left_join(data.simu0, .)

  test1 <- robincar_coxscore(
    df=data.simu,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols="strt",
    covariate_cols=c("model_z1", "model_z2"),
    car_scheme="permuted-block",
    ref_arm=0
  )
  # print(test1)
  test1_ty <- score_robust(data.simu,
                           randomization="permuted_block",
                           p_trt=0.5)
  expect_equal(test1$result$U, test1_ty$U)
  expect_equal(test1$result$se, as.numeric(test1_ty$se))
  expect_equal(test1$result$statistic, as.numeric(test1_ty$T_SR))

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
    dplyr::filter(value == 1) %>%
    dplyr::select(-value) %>%
    dplyr::mutate(strt=forcats::as_factor(strt)) %>%
    dplyr::select(t, strt) %>%
    dplyr::left_join(data.simu01, .)

  test2 <- robincar_coxscore(
    df=data.simu1,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols="strt",
    covariate_cols=c("model_w3"),
    car_scheme="permuted-block",
    ref_arm=0
  )
  test2_ty <- score_robust(data.simu1,
                           randomization="permuted_block",
                           p_trt=0.5)
  expect_equal(test2$result$U, test2_ty$U)
  expect_equal(test2$result$se, as.numeric(test2_ty$se))
  expect_equal(test2$result$statistic, as.numeric(test2_ty$T_SR))

})




test_that("No X no Z case under SR, coxscore",{

  n <- 100
  data.simu0 <- data_gen(
    n=n,
    theta=0,
    randomization="SR",
    p_trt=0.5,
    case="case2"
  )

  data.simu <-
    data.simu0 %>% tidyr::pivot_longer(
      cols=starts_with("car_strata"),
      names_prefix="car_strata",
      names_to="strt"
    ) %>%
    dplyr::filter(value == 1) %>%
    dplyr::select(-value) %>%
    dplyr::mutate(strt=forcats::as_factor(strt)) %>%
    dplyr::select(t, strt) %>%
    dplyr::left_join(data.simu0, .)

  test1 <- robincar_coxscore(
    df=data.simu,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols=NULL,
    covariate_cols=NULL,
    car_scheme="simple",
    ref_arm=0
  )

  expect_true(is.numeric(test1$result$"statistic"))

})
