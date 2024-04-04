
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

  formula <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    formula="y ~ A + x1 + z1",
    car_scheme="simple",
    g_family=binomial(link=logit),
    g_accuracy=7
  )

  formula2 <- robincar_glm2(
    df=DATA2,
    response_col="y",
    treat_col="A",
    formula="y ~ A + x1 + z1",
    car_scheme="simple",
    g_family=binomial(link=logit),
    g_accuracy=7
  )

  expect_equal(formula$result$estimate, formula2$result$estimate)
  expect_equal(formula$varcov, formula2$varcov)
})


test_that("GLM formula versus homogeneous - car", {
  for(car_scheme in c("pocock-simon", "permuted-block", "biased-coin")){

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

    formula2 <- robincar_glm2(
      df=DATA2,
      response_col="y",
      treat_col="A",
      formula="y ~ A + x1 + z1",
      car_strata_cols=c("z1"),
      car_scheme=car_scheme,
      g_family=binomial(link=logit),
      g_accuracy=7
    )

    expect_equal(formula$result$estimate, formula2$result$estimate)
    expect_equal(formula$varcov, formula2$varcov)
  }
})

test_that("GLM formula versus heterogeneous - simple", {

  formula <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    formula="y ~ A + A*x1 + A*z1",
    car_scheme="simple",
    g_family=binomial(link=logit),
    g_accuracy=7
  )

  formula2 <- robincar_glm2(
    df=DATA2,
    response_col="y",
    treat_col="A",
    formula="y ~ A + A*x1 + A*z1",
    car_scheme="simple",
    g_family=binomial(link=logit),
    g_accuracy=7
  )

  expect_equal(formula$result$estimate, formula2$result$estimate)
  expect_equal(formula$varcov, formula2$varcov)
})

test_that("GLM formula versus homogeneous - car", {
  for(car_scheme in c("pocock-simon", "permuted-block", "biased-coin")){

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

    formula2 <- robincar_glm2(
      df=DATA2,
      response_col="y",
      treat_col="A",
      formula="y ~ A + A*x1 + A*z1",
      car_strata_cols=c("z1"),
      car_scheme=car_scheme,
      g_family=binomial(link=logit),
      g_accuracy=7
    )

    expect_equal(formula$result$estimate, formula2$result$estimate)
    expect_equal(formula$varcov, formula2$varcov)
  }
})

test_that("GLM formula versus homogeneous - simple, no Z cov", {

  formula <- robincar_glm(
    df=DATA2,
    response_col="y",
    treat_col="A",
    formula="y ~ A + x1",
    car_scheme="simple",
    g_family=binomial(link=logit),
    g_accuracy=7
  )

  formula2 <- robincar_glm2(
    df=DATA2,
    response_col="y",
    treat_col="A",
    formula="y ~ A + x1",
    car_scheme="simple",
    g_family=binomial(link=logit),
    g_accuracy=7
  )

  expect_equal(formula$result$estimate, formula2$result$estimate)
  expect_equal(formula$varcov, formula2$varcov)
})


test_that("GLM formula versus homogeneous - car, no cov Z", {
  for(car_scheme in c("pocock-simon", "permuted-block", "biased-coin")){

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

    formula2 <- robincar_glm2(
      df=DATA2,
      response_col="y",
      treat_col="A",
      formula="y ~ A + x1",
      car_strata_cols=c("z1"),
      car_scheme=car_scheme,
      g_family=binomial(link=logit),
      g_accuracy=7
    )

    expect_equal(formula$result$estimate, formula2$result$estimate)
    expect_equal(formula$varcov, formula2$varcov)
  }
})

