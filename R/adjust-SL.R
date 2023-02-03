
#' Perform GLM adjustment, based on the classes
#' of the model. Will perform adjustment based on the linear
#' model type of `model` and also do G-computation or AIPW
#' based on the second model type of `model`.
adjust.SLModel <- function(model, data){

  dmat <- get.dmat(data, model$adj_vars)

  # TODO: Need to have option for multiple treatment groups
  # TODO: This fails if there are not enough people in each
  #       treatment group... how to fix?
  aipw <- AIPW$new(
    Y=data$response,
    A=as.integer(data$treat) - 1,
    W.Q=dmat,
    W.g=rep(1, data$n),
    Q.SL.library=model$SL_libraries,
    g.SL.library="SL.mean",
    k_split=model$k_split,
    verbose=T
  )
  # TODO: Need option for whether or not to perform a stratified fit
  # maybe in the adj_method option for homogeneous or heterogeneous?
  aipw$stratified_fit()
  mutilde <- do.call(cbind, aipw$obs_est[1:data$k])
  estimate <- colMeans(mutilde)

  # This is just a placeholder for the data frame, it's not
  # actually using this GLM model fit
  df <- data.frame(response=data$response, treat=data$treat, dmat)
  glmod <- glm(response ~ ., data=df)

  # Compute the asymptotic variance -- not specific to a partition
  asympt.variance <- vcov_car(model, data, glmod, mutilde)
  # TODO: What do we put for the rank of the AIPW model?
  vcov_wt <- get.vcovHC(model$vcovHC, n=data$n, p=glmod$rank)
  variance <- asympt.variance * vcov_wt / data$n

  result <- format.results(data$treat_levels, estimate, variance)

  return(
    structure(
      class="SLModelResult",
      list(result=result, varcov=variance, settings=model,
           data=data, mod=aipw, mu_a=mutilde)
    )
  )
}
