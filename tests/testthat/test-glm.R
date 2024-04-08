
library(MASS)

DATA <- RobinCar:::data_sim
DATA$A <- as.factor(DATA$A)
DATA$y_bin <- ifelse(DATA$y > 2, 1, 0)

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
  non <- robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A",
    car_scheme="simple",
    g_family=gaussian(link="identity"),
    g_accuracy=7,
    adj_method="homogeneous",
    covariate_to_include_strata=FALSE)
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
  non <- robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_scheme="simple",
    g_family=gaussian(link="identity"),
    g_accuracy=7,
    adj_method="homogeneous",
    covariate_to_include_strata=FALSE)
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
  non <- robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_scheme="simple",
    g_family=gaussian(link="identity"),
    g_accuracy=7,
    adj_method="heterogeneous",
    covariate_to_include_strata=FALSE)
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
    car_scheme="simple",
    adj_method="ANCOVA",
    covariate_to_include_strata=FALSE)
  non <- robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A",
    car_scheme="simple",
    covariate_cols=c("x1"),
    car_strata_cols=c("z1"),
    g_family=gaussian(link="identity"),
    g_accuracy=7,
    adj_method="homogeneous",
    covariate_to_include_strata=FALSE)
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
    adj_method="ANCOVA",
    covariate_to_include_strata=FALSE)
  non <- robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A",
    car_scheme="permuted-block",
    covariate_cols=c("x1"),
    car_strata_cols=c("z1"),
    g_family=gaussian(link="identity"),
    g_accuracy=7,
    adj_method="homogeneous",
    covariate_to_include_strata=FALSE)
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

# GLM TESTS ---------------------------------------- #

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

test_that("GLM full function -- NEGATIVE binomial, permuted block", {

  # Known dispersion parameter
  non <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    car_strata_cols=c("z1"),
    covariate_cols=c("x1"),
    car_scheme="permuted-block",
    g_family=negative.binomial(1),
    g_accuracy=7,
    adj_method="heterogeneous",
    covariate_to_include_strata=TRUE)
  expect_equal(class(non), "GLMModelResult")

  # Known dispersion parameter
  non <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    car_strata_cols=c("z1"),
    covariate_cols=c("x1"),
    car_scheme="permuted-block",
    g_family="nb",
    g_accuracy=7,
    adj_method="heterogeneous",
    covariate_to_include_strata=TRUE)
  expect_equal(class(non), "GLMModelResult")

  non

  # Test with formula
  non <- robincar_glm2(
    df=DATA2,
    response_col="y",
    treat_col="A",
    car_strata_cols=c("z1"),
    formula="y ~ A + z1",
    car_scheme="permuted-block",
    g_family="nb",
    g_accuracy=7)

  non

})

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
    covariate_to_include_strata=FALSE)
  expect_equal(non$settings$adj_vars, "x")
  expect_equal(non$settings$pu_joint_z, FALSE)
  expect_equal(non$settings$adj_se_z, FALSE)

  non <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_strata_cols=c("z1"),
    car_scheme="biased-coin",
    g_family=gaussian(link="identity"),
    g_accuracy=7,
    adj_method="heterogeneous",
    covariate_to_include_strata=TRUE)
  expect_equal(non$settings$adj_vars, "joint_z_x")
  expect_equal(non$settings$pu_joint_z, FALSE)
  expect_equal(non$settings$adj_se_z, TRUE)
})

test_that("GLM full function -- binomial, permuted block", {
  non <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    car_strata_cols=c("z1"),
    covariate_cols=c("x1"),
    car_scheme="permuted-block",
    g_family=binomial(link="logit"),
    g_accuracy=7,
    adj_method="heterogeneous",
    covariate_to_include_strata=TRUE)
  expect_equal(class(non), "GLMModelResult")
  expect_equal(non$result$estimate,
               c(X1=0.20774694, X2=0.15547416), tolerance=1e-5)
})


test_that("GLM full function -- binomial, pocock simon", {
  non <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    car_strata_cols=c("z1"),
    covariate_cols=c("x1"),
    car_scheme="pocock-simon",
    g_family=binomial(link="logit"),
    g_accuracy=7,
    adj_method="heterogeneous",
    covariate_to_include_strata=TRUE)
  expect_equal(class(non), "GLMModelResult")
  expect_equal(non$result$estimate,
               c(0.20774694, 0.15547416), tolerance=1e-5)
})

test_that("GLM -- no covariates", {
  non <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    car_strata_cols=c("z1"),
    covariate_to_include_strata=FALSE,
    car_scheme="biased-coin",
    g_family=binomial(link="logit"),
    g_accuracy=7,
    adj_method="homogeneous")
  expect_equal(class(non), "GLMModelResult")
  expect_equal(length(non$mod$coefficients), 2)
})

# TEST FOR THE DMAT.
test_that("GLM -- no covariates except car_strata", {
  non <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    car_scheme="biased-coin",
    car_strata_cols=c("z1"),
    covariate_to_include_strata=TRUE,
    g_family=binomial(link="logit"),
    g_accuracy=7,
    adj_method="heterogeneous")
  expect_equal(length(non$mod$coefficients), 4)
})

n <- 2000
DATA4 <- data.frame(
  y=rbinom(n=n, prob=0.5, size=1),
  TRT01P=sample(1:3, replace=TRUE, size=n),
  BWTGR1=rbinom(n=n, prob=0.1, size=1)
)
DATA4$TRT01P <- factor(DATA4$TRT01P)

test_that("GLM -- test g_family types in print function", {
  # Test with character
  res1 <- robincar_glm(
    df=DATA4,
    response_col="y",
    treat_col="TRT01P",
    car_scheme="simple",
    covariate_cols=c("BWTGR1"),
    g_family="binomial",
    adj_method="homogeneous"
  )
  print(res1)
  # Test with function
  res2 <- robincar_glm(
    df=DATA4,
    response_col="y",
    treat_col="TRT01P",
    car_scheme="simple",
    covariate_cols=c("BWTGR1"),
    g_family=binomial,
    adj_method="homogeneous"
  )
  print(res2)
  # Test with object
  res3 <- robincar_glm(
    df=DATA4,
    response_col="y",
    treat_col="TRT01P",
    car_scheme="simple",
    covariate_cols=c("BWTGR1"),
    g_family=binomial(link="logit"),
    adj_method="homogeneous"
  )
  print(res3)
  expect_equal(res1$result, res2$result)
  expect_equal(res1$result, res3$result)
})

DATA5 <- DATA4
DATA5$TRT01P <- factor(DATA5$TRT01P,
                       levels=1:3,
                       labels=c("placebo", "trt1", "trt2"))

test_that("GLM -- test that it does not matter if treatment levels
          are labeled or not.", {
  # Un-labeled treatment groups
  res1 <- robincar_glm(
    df=DATA4,
    response_col="y",
    treat_col="TRT01P",
    car_scheme="simple",
    covariate_cols=c("BWTGR1"),
    g_family=binomial(link="logit"),
    adj_method="homogeneous"
  )
  # Labeled treatment groups
  res2 <- robincar_glm(
    df=DATA5,
    response_col="y",
    treat_col="TRT01P",
    car_scheme="simple",
    covariate_cols=c("BWTGR1"),
    g_family=binomial(link="logit"),
    adj_method="homogeneous"
  )
  expect_equal(res1$result$estimate, res2$result$estimate)
  expect_equal(res1$result$se, res2$result$se)
})
