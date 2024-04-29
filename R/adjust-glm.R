predictions <- function(model, data, mod){
  UseMethod("predictions", model)
}

#' @exportS3Method
predictions.GLMModel <- function(model, data, mod, ...){
  df <- data.frame(
    treat=data$treat,
    response=data$response
  )
  dmat <- get.dmat(data, model$adj_vars)
  if(!is.null(dmat)){
    df <- cbind(df, dmat)
  }
  preds <- stats::predict(mod, newdata=df, type="response")
  return(preds)
}

get.muhat <- function(model, data, mod){
  set.treat <- function(a){
    dat <- data
    dat$treat <- rep(a, data$n)
    dat$treat <- factor(dat$treat, levels=data$treat_levels)
    return(dat)
  }
  pred.treat <- function(dat) predictions.GLMModel(model, dat, mod)

  datas <- lapply(data$treat_levels, set.treat)
  preds <- lapply(datas, pred.treat)

  pred_cols <- do.call(cbind, preds)
  return(pred_cols)
}

# Perform GLM adjustment, based on the classes
# of the model. Will perform adjustment based on the linear
# model type of `model` and also do G-computation or AIPW
# based on the second model type of `model`.
#' @exportS3Method
adjust.GLMModel <- function(model, data, ...){

  # Get the generalized linear model and estimates by AIPW
  glmod <- linmod(model, data, family=model$g_family, center=FALSE)

  # Get g-computation predictions
  muhat <- get.muhat(model, data, glmod)
  g.estimate <- colMeans(muhat)

  # Get mutilde from the GLM model, then estimate the treatment means by
  # taking the mean over all of the potential outcomes
  glmod.adjusted <- get.mutilde(model, data, muhat)
  mutilde <- glmod.adjusted$mutilde

  # Degree of freedom adjustment
  df_adjust <- glmod.adjusted$df_adjust

  estimate <- colMeans(mutilde)

  # Compute the asymptotic variance
  asympt.variance <- vcov_car(model, data, glmod, mutilde)

  vcov_wt <- get.vcovHC(model$vcovHC, n=data$n, p=df_adjust)
  variance <- asympt.variance * vcov_wt / data$n

  result <- format_results(data$treat_levels, estimate, variance)

  return(
    structure(
      class="GLMModelResult",
      list(result=result, varcov=variance, settings=model,
           data=data, mod=glmod, mu_a=mutilde, g.estimate=g.estimate)
    )
  )
}
