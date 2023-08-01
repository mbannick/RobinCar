## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
rm(list=ls())
library(RobinCar)

## -----------------------------------------------------------------------------
data <- RobinCar:::data_sim
data$A <- as.factor(data$A)
 # add some noise as an extra covariate as toy example
data$extra_cov <- rnorm(nrow(data))

## -----------------------------------------------------------------------------
data$y_bin <- ifelse(data$y > 2, 1, 0)

## -----------------------------------------------------------------------------
result <- robincar_glm(
  df=data,
  response_col="y_bin",
  treat_col="A",
  strata_cols=c("z1"),
  covariate_cols=c("x3"),
  car_scheme="biased-coin",
  g_family=binomial(link="logit"),
  covariate_to_include_strata=TRUE,
  g_accuracy=7,
  adj_method="homogeneous"
)
print(result)

## -----------------------------------------------------------------------------
robincar_calibrate(
  result=result,
  joint=FALSE
)

## -----------------------------------------------------------------------------
robincar_calibrate(
  result=result,
  joint=TRUE
)

## -----------------------------------------------------------------------------
robincar_calibrate(
  result=result,
  joint=TRUE,
  add_x=c("extra_cov")
)

