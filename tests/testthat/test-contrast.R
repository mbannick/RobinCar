

test_that("Test contrast function ordering", {

  data <- RobinCar:::data_sim %>%
    # recode A from 0 to "a", 1 to "b", 2 to "c"
    mutate(A=factor(A, levels=0:2, labels=c("b", "c", "a")))

  my.f<-function(x){
    x[2]-x[1] # should be "c"-"b" based on the ordering
  }

  fit.anova <- robincar_linear(df = data,
                               response_col="y",
                               treat_col="A",
                               car_scheme="simple",
                               adj_method="ANOVA",
                               contrast_h=my.f)

  # Check levels are in the right order
  lev <- fit.anova$main$result$treat
  expect_equal(lev, c("b", "c", "a"))

  # Check that the estimates are in correct order
  anova <- mean(data$y[data$A == "c"]) - mean(data$y[data$A == "b"])
  expect_equal(as.vector(fit.anova$contrast$result$estimate[1]), anova)

})

test_that("Test warnings for ratio and odds ratio functions", {

  data <- RobinCar:::data_sim %>%
    # recode A from 0 to "a", 1 to "b", 2 to "c"
    mutate(A=factor(A, levels=0:2, labels=c("a", "b", "c"))) %>%
    mutate(y_bin=ifelse(y > mean(y), 1, 0)) # create binary outcome

  warning_txt <- "Using a ratio or odds ratio is not recommended.
            Consider using log_ratio or log_odds_ratio for better
            performance of the variance estimator."

  # Test for linear model with ratio
  expect_warning(
    fit <- robincar_glm(
      df=data,
      response_col="y_bin",
      treat_col="A",
      formula= y_bin ~ A,
      contrast_h="odds_ratio"
    ),
    warning_txt
  )

  # Test for linear model with ratio
  expect_warning(
    fit <- robincar_glm(
      df=data,
      response_col="y",
      treat_col="A",
      formula= y ~ A,
      contrast_h="ratio"
    ),
    warning_txt
  )

})



test_that("Test warnings for contrast functions that are not probabilities", {
  # We want these to give warnings when someone
  # specifies that they want a log odds ratio.

  data <- RobinCar:::data_sim %>%
    # recode A from 0 to "a", 1 to "b", 2 to "c"
    mutate(A=factor(A, levels=0:2, labels=c("a", "b", "c"))) %>%
    mutate(y_bin=ifelse(y > mean(y), 1, 0)) # create binary outcome

  # Test for linear model with ratio
  expect_warning(
    fit <- robincar_glm(
      df=data,
      response_col="y",
      treat_col="A",
      formula= y ~ A,
      contrast_h="log_odds_ratio"
    ),
    "Estimates are not between 0 and 1: are you sure you want an odds ratio?"
  )

})

test_that("Compare log odds ratio SEs to GLM", {
  # We want these to give warnings when someone
  # specifies that they want a log odds ratio.

  data <- RobinCar:::data_sim %>%
    # recode A from 0 to "a", 1 to "b", 2 to "c"
    mutate(A=factor(A, levels=0:2, labels=c("a", "b", "c"))) %>%
    mutate(y_bin=ifelse(y > mean(y), 1, 0)) # create binary outcome

  # Test for linear model with ratio
  fit <- robincar_glm(
    df=data,
    response_col="y_bin",
    treat_col="A",
    formula= y_bin ~ A,
    contrast_h="log_odds_ratio"
  )

  or <- with(data, glm(y_bin ~ A, family=binomial(link="logit")))

  # Check that log odds ratios match
  expect_equal(
    unname(or$coefficients[2:3]),
    unname(fit$contrast$result$estimate), tolerance=1e-6)

})

