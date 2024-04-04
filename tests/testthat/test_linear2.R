library(MASS)

DATA <- RobinCar:::data_sim
DATA$A <- as.factor(DATA$A)

# LINEAR V. GLM WITHOUT Z ----------------------------------- #

test_that("GLM full function -- linear (ANOVA)", {
  lin <- robincar_linear(
    df=DATA,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_scheme="simple",
    adj_method="ANOVA",
    covariate_to_include_strata=FALSE)
  non <- robincar_linear2(
    df=DATA,
    response_col="y",
    treat_col="A",
    car_scheme="simple",
    adj_method="ANOVA")
  expect_equal(class(non), "GLMModelResult")

  # Check that the result from the linear and glm function is the same
  # when using the identity link.
  est_lin <- lin$result$estimate
  names(est_lin) <- NULL
  est_non <- non$result$estimate
  names(est_non) <- NULL
  expect_equal(est_lin, est_non, tolerance=1e-10)

  # These won't be exactly equivalent, only asymptotically (?)
  var_lin <- lin$result$se
  names(var_lin) <- NULL
  var_non <- non$result$se
  names(var_non) <- NULL
  expect_equal(var_lin, var_non, tolerance=1e-2)
})

test_that("GLM full function -- linear (ANCOVA)", {
  lin <- robincar_linear(
    df=DATA,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_scheme="simple",
    adj_method="ANCOVA",
    covariate_to_include_strata=FALSE)
  non <- robincar_linear2(
    df=DATA,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_scheme="simple",
    adj_method="ANCOVA")
  expect_equal(class(non), "GLMModelResult")

  # Check that the result from the linear and glm function is the same
  # when using the identity link.
  est_lin <- lin$result$estimate
  names(est_lin) <- NULL
  est_non <- non$result$estimate
  names(est_non) <- NULL
  expect_equal(est_lin, est_non, tolerance=1e-10)

  # These won't be exactly equivalent, only asymptotically (?)
  var_lin <- lin$result$se
  names(var_lin) <- NULL
  var_non <- non$result$se
  names(var_non) <- NULL
  expect_equal(var_lin, var_non, tolerance=1e-2)
})

test_that("GLM full function -- linear (ANHECOVA)", {
  lin <- robincar_linear(
    df=DATA,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_scheme="simple",
    adj_method="ANHECOVA",
    covariate_to_include_strata=FALSE)
  non <- robincar_linear2(
    df=DATA,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_scheme="simple",
    adj_method="ANHECOVA")
  expect_equal(class(non), "GLMModelResult")

  # Check that the result from the linear and glm function is the same
  # when using the identity link.
  est_lin <- lin$result$estimate
  names(est_lin) <- NULL
  est_non <- non$result$estimate
  names(est_non) <- NULL
  expect_equal(est_lin, est_non, tolerance=1e-10)

  # These won't be exactly equivalent, only asymptotically (?)
  var_lin <- lin$result$se
  names(var_lin) <- NULL
  var_non <- non$result$se
  names(var_non) <- NULL
  expect_equal(var_lin, var_non, tolerance=1e-2)
})

# LINEAR V. GLM WITH Z ----------------------------------- #

test_that("GLM full function -- linear (ANCOVA) with Z", {
  lin <- robincar_linear(
    df=DATA,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_strata_cols=c("z1"),
    car_scheme="pocock-simon",
    adj_method="ANCOVA",
    covariate_to_include_strata=FALSE)
  non <- robincar_linear2(
    df=DATA,
    response_col="y",
    treat_col="A",
    car_scheme="pocock-simon",
    covariate_cols=c("x1", "z1"),
    car_strata_cols=c("z1"),
    adj_method="ANCOVA")
  expect_equal(class(non), "GLMModelResult")

  # Check that the result from the linear and glm function is the same
  # when using the identity link.
  est_lin <- lin$result$estimate
  names(est_lin) <- NULL
  est_non <- non$result$estimate
  names(est_non) <- NULL
  expect_equal(est_lin, est_non, tolerance=1e-10)

  # These won't be exactly equivalent, only asymptotically (?)
  var_lin <- lin$result$se
  names(var_lin) <- NULL
  var_non <- non$result$se
  names(var_non) <- NULL
  expect_equal(var_lin, var_non, tolerance=1e-2)
})

test_that("GLM full function -- linear (ANCOVA) with Z", {
  lin <- robincar_linear(
    df=DATA,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_strata_cols=c("z1"),
    car_scheme="permuted-block",
    adj_method="ANHECOVA",
    covariate_to_include_strata=FALSE)
  non <- robincar_linear2(
    df=DATA,
    response_col="y",
    treat_col="A",
    car_scheme="permuted-block",
    covariate_cols=c("x1", "z1"),
    car_strata_cols=c("z1"),
    adj_method="ANHECOVA")
  expect_equal(class(non), "GLMModelResult")

  # Check that the result from the linear and glm function is the same
  # when using the identity link.
  est_lin <- lin$result$estimate
  names(est_lin) <- NULL
  est_non <- non$result$estimate
  names(est_non) <- NULL
  expect_equal(est_lin, est_non, tolerance=1e-10)

  # These won't be exactly equivalent, only asymptotically (?)
  var_lin <- lin$result$se
  names(var_lin) <- NULL
  var_non <- non$result$se
  names(var_non) <- NULL
  expect_equal(var_lin, var_non, tolerance=1e-2)
})
