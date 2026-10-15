make_tte_df <- function(n=100){
  set.seed(0)
  df <- data_gen(
    n=n, theta=0, randomization="permuted_block", p_trt=0.5, case="case2"
  )
  df$strata1 <- sample(letters[1:3], n, replace=TRUE)
  df
}

run_tte <- function(df, covariate_cols, ...){
  robincar_logrank(
    df=df,
    treat_col="I1",
    response_col="t",
    event_col="delta",
    car_strata_cols="strata1",
    covariate_cols=covariate_cols,
    car_scheme="permuted-block",
    adj_method="CL",
    ...
  )
}

test_that("covariates named s0 or s1 raise an error", {
  df <- make_tte_df()
  for(nm in c("s0", "s1")){
    d <- df
    d[[nm]] <- d$model_z1
    expect_error(run_tte(d, c("model_z2", nm)), "reserved names")
  }
})

test_that("all reserved covariate names raise an error", {
  df <- make_tte_df()
  for(nm in .TTE.RESERVED.NAMES){
    d <- df
    d[[nm]] <- d$model_z1
    expect_error(run_tte(d, c("model_z2", nm)), "reserved names", info=nm)
  }
})

test_that("reserved prefixes raise an error", {
  df <- make_tte_df()
  df$xmat_a <- df$model_z1
  df$trt1_a <- df$model_z1
  expect_error(run_tte(df, "xmat_a"), "reserved names")
  expect_error(run_tte(df, "trt1_a"), "reserved names")
})

test_that("reserved names for strata and wrong roles raise an error", {
  df <- make_tte_df()
  df$s1 <- df$strata1
  expect_error(
    robincar_logrank(
      df=df, treat_col="I1", response_col="t", event_col="delta",
      car_strata_cols="s1", car_scheme="permuted-block", adj_method="CL"
    ), "reserved name")
  df$event <- df$delta
  expect_error(
    robincar_logrank(
      df=df, treat_col="I1", response_col="event", event_col="delta",
      car_strata_cols="strata1", car_scheme="permuted-block", adj_method="CL"
    ), "reserved name")
})

test_that("valid covariate names do not error", {
  df <- make_tte_df()
  expect_error(run_tte(df, c("model_z1", "model_z2")), NA)
})
