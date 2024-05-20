# Fits a working model based on GLM and computes AIPW estimator

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

# Perform GLM adjustment, based on the classes
# of the model. Will perform adjustment based on the linear
# model type of `model` and also do G-computation or AIPW
# based on the second model type of `model`.
#' @importFrom stats setNames
#' @exportS3Method
adjust.GLMModel <- function(model, data, ...){

  # 1. Set up data frame
  df <- data$df
  columns <- colnames(df)
  columns[which(columns == data$treat_col)] <- "treat"
  columns[which(columns == data$response_col)] <- "response"
  df <- setNames(df, columns)

  # 2. Fit GLM working model
  glmod <- fitmod(
    family=model$g_family,
    stats::as.formula(data$formula),
    data=df)

  # 3. Predict potential outcome for each treatment group
  # to compute \hat{\mu}_a for a = 1, ..., k
  predict_a <- function(a){

    # Set up data frame
    df_a <- df
    df_a$treat <- rep(a, data$n)
    df_a$treat <- factor(df_a$treat, levels=data$treat_levels)

    # Create predictions
    predictions <- stats::predict(glmod, newdata=df_a, type="response")

    return(predictions)
  }

  preds <- lapply(data$treat_levels, predict_a)
  muhat <- do.call(cbind, preds)

  # 4. Save the g-computation estimate of \theta for diagnostic purposes
  g.estimate <- colMeans(muhat)

  # 5. Use the \hat{\mu} to compute the adjusted \hat{\mu} for AIPW estimator
  glmod.adjusted <- get.mutilde(model, data, muhat)
  mutilde <- glmod.adjusted$mutilde

  # 6. Compute AIPW estimator of \theta
  estimate <- colMeans(mutilde)

  # 7. Compute the asymptotic variance
  asympt.variance <- vcov_car(model, data, glmod, mutilde)
  df_adjust <- glmod.adjusted$df_adjust
  vcov_wt <- get.vcovHC(model$vcovHC, n=data$n, p=df_adjust)
  variance <- asympt.variance * vcov_wt / data$n

  # 8. Format results
  result <- format_results(data$treat_levels, estimate, variance)

  return(
    structure(
      class="GLMModelResult",
      list(result=result, varcov=variance, settings=model,
           data=data, mod=glmod, mu_a=mutilde, g.estimate=g.estimate)
    )
  )
}
