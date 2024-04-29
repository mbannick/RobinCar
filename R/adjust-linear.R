# Fits a linear model with settings based on model and data.
# The linear model is not the "true" model, but it's the linear model
# that will aid in adjustment for ANCOVA and ANHECOVA, and used
# to estimate the treatment effects for ANOVA (equivalent to a sample)
# mean per group.

# This function allows us to implement the negative binomial model through MASS
#' @importFrom MASS glm.nb
fitmod <- function(family, ...){
  useMASS <- FALSE
  if(is.character(family)){
    if(family == "nb"){
      useMASS <- TRUE
    }
  }
  if(useMASS){
    mod <- MASS::glm.nb(...)
    # copy over the 'model' attribute to the 'data'
    # because it has a different name in glm.nb
    mod$data <- mod$model
  } else {
    mod <- stats::glm(..., family=family)
  }
  return(mod)
}

#' @importFrom stats gaussian
linmod <- function(model, data, family, center=TRUE){
  UseMethod("linmod", model)
}

#' @importFrom stats gaussian
#' @exportS3Method
linmod.ANOVA <- function(model, data, family=stats::gaussian, center=TRUE){
  df <- data.frame(
    treat=data$treat,
    response=data$response
  )
  mod <- fitmod(family=family, formula=response ~ 0 + treat, data=df)
  return(mod)
}

#' @importFrom stats gaussian
#' @exportS3Method
linmod.ANCOVA <- function(model, data, family=stats::gaussian, center=TRUE){
  df <- data.frame(
    treat=data$treat,
    response=data$response
  )
  dmat <- get.dmat(data, model$adj_vars)
  if(center) dmat <- .center.dmat(dmat)
  df <- cbind(df, dmat)
  if(center){
    mod <- fitmod(family=family, formula=response ~ 0 + treat + ., data=df)
  } else {
    mod <- fitmod(family=family, formula=response ~ 1 + treat + ., data=df)
  }
  return(mod)
}

#' @exportS3Method
linmod.ANHECOVA <- function(model, data, family=stats::gaussian, center=TRUE){
  df <- data.frame(
    treat=data$treat,
    response=data$response
  )
  dmat <- get.dmat(data, model$adj_vars)
  if(center) dmat <- .center.dmat(dmat)
  df <- cbind(df, dmat)

  if(center){
    mod <- fitmod(family=family, response ~ 0 + treat:., data=df)
  } else {
    mod <- fitmod(family=family, response ~ 1 + treat:., data=df)
  }
  return(mod)
}

#' @exportS3Method
linmod.CUSTOM <- function(model, data, family=stats::gaussian, center=TRUE){
  df <- data.frame(
    treat=data$treat,
    response=data$response
  )
  dmat <- get.dmat(data, model$adj_vars)
  if(center) dmat <- .center.dmat(dmat)
  df <- cbind(df, dmat)
  mod <- fitmod(family=family, stats::as.formula(data$formula), data=df)

  return(mod)
}

# Perform adjustment for linear models, including ANOVA, ANCOVA,
# and ANHECOVA, with or without covariate-adaptive randomization.
#' @exportS3Method
adjust.LinModel <- function(model, data, ...){

  # Fit a model with the settings in model
  mod <- linmod(model, data, family=gaussian)

  # Get the simple randomization variance and adjust if necessary
  asympt.variance <- vcov_car(model, data, mod)
  vcov_wt <- get.vcovHC(model$vcovHC, n=data$n, p=mod$rank)
  variance <- asympt.variance * vcov_wt / data$n

  # Extract estimates and create results data
  est <- stats::coef(mod)[1:data$k]
  result <- format_results(data$treat_levels, est, variance)

  return(
    structure(
      class="LinModelResult",
      list(result=result, varcov=variance, settings=model,
           data=data, mod=mod)
    )
  )
}
