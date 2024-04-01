
n <- 1000
set.seed(10)
DATA2 <- data.frame(A1=rbinom(n, size=1, prob=0.5),
                    A2=rbinom(n, size=1, prob=0.5),
                    y=rbinom(n, size=1, prob=0.2),
                    x1=rnorm(n),
                    x2=rnorm(n),
                    x3=as.factor(rbinom(n, size=1, prob=0.5)),
                    z1=rbinom(n, size=1, prob=0.5),
                    z2=rbinom(n, size=1, prob=0.5))
DATA2$A <- as.factor(DATA2$A1)
DATA2$x1 <- DATA2$x1 - mean(DATA2$x1)
DATA2$A3LEVEL <- as.factor(DATA2$A1 + DATA2$A2)

test_that("GLM formula versus homogeneous - simple", {
  specs <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1", "z1"),
    car_scheme="simple",
    g_family=binomial(link="logit"),
    g_accuracy=7,
    adj_method="homogeneous",
    covariate_to_include_strata=TRUE)

  formula <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    formula="y ~ A + x1 + z1",
    car_scheme="simple",
    g_family=binomial(link=logit),
    g_accuracy=7
  )

  expect_equal(specs$result$estimate, formula$result$estimate)
  expect_equal(specs$varcov, formula$varcov)
})


test_that("GLM formula versus homogeneous - car", {
  for(car_scheme in c("pocock-simon", "permuted-block", "biased-coin")){
    specs <- robincar_glm(
      df=DATA2,
      response_col="y",
      treat_col="A",
      car_strata_cols=c("z1"),
      covariate_cols=c("x1"),
      car_scheme=car_scheme,
      g_family=binomial(link="logit"),
      g_accuracy=7,
      adj_method="homogeneous",
      covariate_to_include_strata=TRUE)

    formula <- robincar_glm(
      df=DATA2,
      response_col="y",
      treat_col="A",
      formula="y ~ A + x1 + z1",
      car_strata_cols=c("z1"),
      car_scheme=car_scheme,
      g_family=binomial(link=logit),
      g_accuracy=7
    )

    expect_equal(specs$result$estimate, formula$result$estimate)
    expect_equal(specs$varcov, formula$varcov)
  }
})

test_that("GLM formula versus heterogeneous - simple", {
  specs <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1", "z1"),
    car_scheme="simple",
    g_family=binomial(link="logit"),
    g_accuracy=7,
    adj_method="heterogeneous",
    covariate_to_include_strata=TRUE)

  formula <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    formula="y ~ A + A*x1 + A*z1",
    car_scheme="simple",
    g_family=binomial(link=logit),
    g_accuracy=7
  )

  expect_equal(specs$result$estimate, formula$result$estimate)
  expect_equal(specs$varcov, formula$varcov)
})

test_that("GLM formula versus homogeneous - car", {
  for(car_scheme in c("pocock-simon", "permuted-block", "biased-coin")){
    specs <- robincar_glm(
      df=DATA2,
      response_col="y",
      treat_col="A",
      car_strata_cols=c("z1"),
      covariate_cols=c("x1"),
      car_scheme=car_scheme,
      g_family=binomial(link="logit"),
      g_accuracy=7,
      adj_method="heterogeneous",
      covariate_to_include_strata=TRUE)

    formula <- robincar_glm(
      df=DATA2,
      response_col="y",
      treat_col="A",
      formula="y ~ A + A*x1 + A*z1",
      car_strata_cols=c("z1"),
      car_scheme=car_scheme,
      g_family=binomial(link=logit),
      g_accuracy=7
    )

    expect_equal(specs$result$estimate, formula$result$estimate)
    expect_equal(specs$varcov, formula$varcov)
  }
})

test_that("GLM formula versus homogeneous - simple, no Z cov", {
  specs <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    covariate_cols=c("x1"),
    car_scheme="simple",
    g_family=binomial(link="logit"),
    g_accuracy=7,
    adj_method="homogeneous")

  formula <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    formula="y ~ A + x1",
    car_scheme="simple",
    g_family=binomial(link=logit),
    g_accuracy=7
  )

  expect_equal(specs$result$estimate, formula$result$estimate)
  expect_equal(specs$varcov, formula$varcov)
})


test_that("GLM formula versus homogeneous - car, no cov Z", {
  for(car_scheme in c("pocock-simon", "permuted-block", "biased-coin")){

    # TODO: What desired behavior do we want here for pocock-simon?
    # The two calls below do not produce the same results because
    # we don't enforce that formula uses z as covariates.

    specs <- robincar_glm(
      df=DATA2,
      response_col="y",
      treat_col="A",
      car_strata_cols=c("z1"),
      covariate_cols=c("x1"),
      car_scheme=car_scheme,
      g_family=binomial(link="logit"),
      g_accuracy=7,
      covariate_to_include_strata=FALSE,
      adj_method="homogeneous")

    formula <- robincar_glm(
      df=DATA2,
      response_col="y",
      treat_col="A",
      formula="y ~ A + x1",
      car_strata_cols=c("z1"),
      car_scheme=car_scheme,
      g_family=binomial(link=logit),
      g_accuracy=7
    )

    expect_equal(specs$result$estimate, formula$result$estimate)
    expect_equal(specs$varcov, formula$varcov)
  }
})

test_that("GLM formula - POCOCK SIMON no X", {

  specs <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    car_strata_cols=c("z1"),
    car_scheme="pocock-simon",
    g_family=binomial(link="logit"),
    g_accuracy=7,
    covariate_to_include_strata=FALSE,
    adj_method="heterogeneous")

  formula <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    formula="y ~ A",
    car_strata_cols=c("z1"),
    car_scheme="pocock-simon",
    g_family=binomial(link=logit),
    g_accuracy=7
  )
  # if we calibrate formula with joint Z, then this should be equivalent
  # to the initial specs model
  calib <- robincar_calibrate(
    formula, joint=TRUE
  )

  expect_equal(specs$result$estimate, calib$result$estimate)
  expect_equal(specs$varcov, calib$varcov)

  expect_equal(specs$result$estimate, formula$result$estimate)
  # TODO: But why would these be different but estimates the same?
  # it's because of using HC3. Calibrating within joint levels of Z using our pocock-simon trick
  # to get variance does not change the "p" of the glm, so HC3 differs. HC0 does not.
  # Do we need to do anything or is this ok?
  #
  # Skipping HC0 for now.
  expect_equal(specs$varcov, formula$varcov)
})

test_that("GLM formula - POCOCK SIMON no X, 3 levels of A", {

  specs <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A3LEVEL",
    car_strata_cols=c("z1"),
    car_scheme="pocock-simon",
    g_family=binomial(link="logit"),
    g_accuracy=7,
    covariate_to_include_strata=FALSE,
    adj_method="heterogeneous")

  formula <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A3LEVEL",
    formula="y ~ A3LEVEL",
    car_strata_cols=c("z1"),
    car_scheme="pocock-simon",
    g_family=binomial(link=logit),
    g_accuracy=7
  )
  # if we calibrate formula with joint Z, then this should be equivalent
  # to the initial specs model
  calib <- robincar_calibrate(
    formula, joint=TRUE
  )

  expect_equal(specs$result$estimate, calib$result$estimate)
  expect_equal(specs$varcov, calib$varcov)

  expect_equal(specs$result$estimate, formula$result$estimate)
  # TODO: But why would these be different but estimates the same?
  # it's because of using HC3. Calibrating within joint levels of Z using our pocock-simon trick
  # to get variance does not change the "p" of the glm, so HC3 differs. HC0 does not.
  # Do we need to do anything or is this ok?
  expect_equal(specs$varcov, formula$varcov)
})

