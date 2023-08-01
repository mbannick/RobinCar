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

## -----------------------------------------------------------------------------
fit.anova <- robincar_linear(df = data, 
                             response_col="y",
                             treat_col="A",
                             covariate_cols=c("x1", "x3"),
                             car_scheme="simple",
                             adj_method="ANOVA",
                             vcovHC="HC0")

## -----------------------------------------------------------------------------
fit.ancova <- robincar_linear(df = data, 
                              response_col="y",
                              treat_col="A",
                              strata_cols=c("z1", "z2"),
                              covariate_cols=c("x1", "x3"),
                              car_scheme="biased-coin",
                              adj_method="ANCOVA",
                              vcovHC="HC0")

## -----------------------------------------------------------------------------
fit.anhecova <- robincar_linear(df = data, 
                                response_col="y",
                                treat_col="A",
                                strata_cols=c("z1", "z2"),
                                covariate_cols=c("x1", "x3"),
                                car_scheme="pocock-simon",
                                adj_method="ANHECOVA",
                                vcovHC="HC0",
                                contrast_h="diff")

## -----------------------------------------------------------------------------
data$y_bin <- ifelse(data$y > 2, 1, 0)

## -----------------------------------------------------------------------------
glm.homogeneous<-robincar_glm(df = data,
                              response_col="y_bin",
                              treat_col="A",
                              strata_cols=c("z1"),
                              car_scheme="pocock-simon",
                              g_family=binomial(link="logit"),
                              g_accuracy=7,
                              adj_method="heterogeneous",
                              covariate_to_include_strata=TRUE,
                              vcovHC="HC3")
glm.homogeneous

## -----------------------------------------------------------------------------
glm.heterogeneous<-robincar_glm(df = data,
                              response_col="y_bin",
                              treat_col="A",
                              strata_cols=c("z2"),
                              covariate_cols=c("x1"),
                              car_scheme="biased-coin",
                              g_family=binomial(link="logit"),
                              g_accuracy=7,
                              covariate_to_include_strata=TRUE,
                              adj_method="heterogeneous",
                              vcovHC="HC1")
glm.heterogeneous

## -----------------------------------------------------------------------------
glm.formula<-robincar_glm(df = data,
                          response_col="y_bin",
                          treat_col="A",
                          formula="y_bin ~ A + z2:x1",
                          strata_cols=c("z2"),
                          car_scheme="biased-coin",
                          g_family=binomial(link="logit"),
                          g_accuracy=7,
                          vcovHC="HC0")
glm.formula

## -----------------------------------------------------------------------------
glm.heterogeneous<-robincar_glm(df = data,
                              response_col="y_bin",
                              treat_col="A",
                              strata_cols=c("z2"),
                              covariate_cols=c("x1"),
                              car_scheme="biased-coin",
                              g_family=binomial(link="logit"),
                              g_accuracy=7,
                              covariate_to_include_strata=TRUE,
                              adj_method="heterogeneous",
                              vcovHC="HC0",
                              contrast_h="diff")
glm.heterogeneous

