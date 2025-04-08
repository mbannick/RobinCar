
library(matrixcalc)

DATA <- RobinCar:::data_sim
DATA$A <- as.factor(DATA$A)
DATA$y_bin <- ifelse(DATA$y > 2, 1, 0)

test_that("ANCOVA for GLM vartypes.", {

  res1 <- robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A",
    car_scheme="simple",
    g_family=gaussian(link="identity"),
    g_accuracy=7,
    formula="y ~ A + x1",
    variance_type=1)

  res2 <- robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A",
    car_scheme="simple",
    g_family=gaussian(link="identity"),
    g_accuracy=7,
    formula="y ~ A + x1", variance_type=2)

  res3 <- robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A",
    car_scheme="simple",
    g_family=gaussian(link="identity"),
    g_accuracy=7,
    formula="y ~ A + x1",
    variance_type=3)

  # These will be equal -- variance type should not
  # change the point estimate
  expect_equal(res1$result$estimate, res2$result$estimate)
  expect_equal(res2$result$estimate, res3$result$estimate)

  # These won't be exactly equivalent, only asymptotically
  # Just make sure things aren't really different
  expect_equal(res1$result$se, res2$result$se, tolerance=1e-2)
  expect_equal(res1$result$se, res3$result$se, tolerance=1e-2)

  # Check that variance type 3 is PD matrix -- this should
  # be true for ANCOVA
  expect_true(is.positive.definite(res3$varcov))

})
