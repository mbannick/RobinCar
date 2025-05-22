library(dplyr)

test_that("data preprocessing", {

  set.seed(0)
  n <- 100

  data.simu0 <- data_gen(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case2"
  ) %>% mutate(
    strata1=sample(letters[1:3], n, replace=TRUE),
    strata2=sample(LETTERS[4:5], n, replace=TRUE)
  )

  # CSL for LOGRANK ----------------------------------------
  # Should have a car_strata with multiple car_strata values

  data_csl <- robincar_logrank(
    df=data.simu0,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols=c("strata1", "strata2"),
    covariate_cols=c("model_z1", "model_z2"),
    car_scheme="permuted-block",
    adj_method="CSL"
  )

  df_csl <- create.tte.df(model=data_csl$settings, data=data_csl$data)
  expect_true(length(unique(df_csl$car_strata)) != 1)
  expect_true(is.factor(df_csl$car_strata))
  expect_true(is.factor(df_csl$carcov_z))

  data_cl <- robincar_logrank(
    df=data.simu0,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols=c("strata1", "strata2"),
    covariate_cols=c("model_z1", "model_z2"),
    car_scheme="permuted-block",
    adj_method="CL"
  )

  df_cl <- create.tte.df(model=data_cl$settings, data=data_cl$data)
  expect_true(length(unique(df_cl$car_strata)) == 1)
  expect_true(is.factor(df_cl$carcov_z))

  data_score <- robincar_coxscore(
    df=data.simu0,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols=c("strata1", "strata2"),
    covariate_cols=c("model_z1", "model_z2"),
    car_scheme="permuted-block",
  )

  df_score <- create.tte.df(model=data_score$settings, data=data_score$data)
  expect_true(length(unique(df_score$car_strata)) == 1)
  expect_true(is.factor(df_score$carcov_z))

  # If there are no covariates
  data_score1 <- robincar_coxscore(
    df=data.simu0,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols=c("strata1", "strata2"),
    car_scheme="permuted-block"
  )
  df_score1 <- create.tte.df(model=data_score$settings, data=data_score$data)
  expect_true(length(unique(df_score1$car_strata)) == 1)
  expect_true(is.factor(df_score1$carcov_z))
  expect_false(all(grepl("robcarx_", colnames(df_score1))))

})

