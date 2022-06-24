
DATA <- RobinCar:::data_sim
DATA$A <- as.factor(DATA$A)
DATA$y_bin <- ifelse(DATA$y > 2, 1, 0)

test_that("GLM full function -- linear", {
  lin <- robincar_linear(
    df=DATA,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_scheme="simple",
    adj_method="ANHECOVA",
    covariate_to_include_strata=FALSE,
    vcovHC="HC0")
  non <- robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_scheme="simple",
    g_family=gaussian(link="identity"),
    g_accuracy=7,
    adj_method="heterogeneous",
    covariate_to_include_strata=FALSE,
    vcovHC="HC0")
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

n <- 1000
set.seed(10)
DATA2 <- data.frame(A=rbinom(n, size=1, prob=0.5),
                 y=rbinom(n, size=1, prob=0.2),
                 x1=rnorm(n),
                 x2=rnorm(n),
                 x3=as.factor(rbinom(n, size=1, prob=0.5)),
                 z1=rbinom(n, size=1, prob=0.5),
                 z2=rbinom(n, size=1, prob=0.5))
DATA2$A <- as.factor(DATA2$A)
DATA2$x1 <- DATA2$x1 - mean(DATA2$x1)

test_that("GLM Settings", {
  non <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_scheme="simple",
    g_family=gaussian(link="identity"),
    g_accuracy=7,
    adj_method="heterogeneous",
    covariate_to_include_strata=FALSE,
    vcovHC="HC0")
  expect_equal(non$settings$adj_vars, "x")
  expect_equal(non$settings$pu_joint_z, FALSE)
  expect_equal(non$settings$adj_se_z, FALSE)

  non <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    strata_cols=c("z1"),
    car_scheme="biased-coin",
    g_family=gaussian(link="identity"),
    g_accuracy=7,
    adj_method="heterogeneous",
    covariate_to_include_strata=TRUE,
    vcovHC="HC0")
  expect_equal(non$settings$adj_vars, "joint_z_x")
  expect_equal(non$settings$pu_joint_z, FALSE)
  expect_equal(non$settings$adj_se_z, TRUE)
})

test_that("GLM full function -- binomial", {
  non <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    strata_cols=c("z1"),
    covariate_cols=c("x1"),
    car_scheme="pocock-simon",
    g_family=binomial(link="logit"),
    g_accuracy=7,
    adj_method="heterogeneous",
    covariate_to_include_strata=TRUE,
    vcovHC="HC3")
  expect_equal(class(non), "GLMModelResult")
})
