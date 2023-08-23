## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
rm(list=ls())
library(RobinCar)
library(dplyr)

## -----------------------------------------------------------------------------
set.seed(10)
n <- 100

# Simulate data
data.simu0 <- data_gen2(
  n=n,
  theta=0,
  randomization="permuted_block",
  p_trt=0.5,
  case="case1"
)

# Restructure data
data.simu <- data.simu0 %>%
  tidyr::pivot_longer(
    cols=starts_with("strata"),
    names_prefix="strata",
    names_to="strt"
  ) %>%
  dplyr::filter(value == 1) %>%
  dplyr::select(-value) %>%
  dplyr::mutate(strt=forcats::as_factor(strt)) %>%
  dplyr::select(t, strt) %>%
  dplyr::left_join(data.simu0, .)

## -----------------------------------------------------------------------------
cox <- robincar_coxscore(
  df=data.simu,
  treat_col="I1",
  response_col="t",
  event_col="delta",
  strata_cols="strt",
  covariate_cols=c("model_w3"),
  car_scheme="permuted-block",
  ref_arm=0
)
print(cox)

## -----------------------------------------------------------------------------
cl <- robincar_logrank(
  df=data.simu,
  treat_col="I1",
  response_col="t",
  event_col="delta",
  strata_cols="strt",
  covariate_cols=c("model_w3"),
  car_scheme="permuted-block",
  adj_method="CL",
  ref_arm=0
)
print(cl)

## -----------------------------------------------------------------------------
csl <- robincar_logrank(
  df=data.simu,
  treat_col="I1",
  response_col="t",
  event_col="delta",
  strata_cols="strt",
  covariate_cols=c("model_w3"),
  car_scheme="permuted-block",
  adj_method="CSL",
  ref_arm=0
)
print(csl)

