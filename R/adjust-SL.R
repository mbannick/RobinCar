
# Perform GLM adjustment, based on the classes
# of the model. Will perform adjustment based on the linear
# model type of `model` and also do G-computation or AIPW
# based on the second model type of `model`.
#' @rdname adjust.SLModel
#' @exportS3Method
adjust.SLModel <- function(model, data, ...){

  SL.libraries <- model$SL_libraries
  for(learner in model$SL_learners){
    SL.libraries <- c(SL.libraries, learner$names)
  }

  # Placeholder for mutilde
  muhat <- matrix(data=NA, nrow=data$n, ncol=data$k)
  browser()
  if(is.null(data$joint_strata_levels) | (model$stratify_fit == F)){
    data$joint_strata_levels <- c("0")
    data$joint_strata <- rep("0", data$n)
  }

  for(sl in data$joint_strata_levels){

    indices <- which(data$joint_strata == sl)

    dmat <- get.dmat(data, model$adj_vars, indices=indices)

    # TODO: Need to have option for multiple treatment groups
    # TODO: This fails if there are not enough people in each
    #       treatment group... how to fix?
    aipw <- AIPW$new(
      Y=data$response[indices],
      A=as.integer(data$treat[indices]) - 1,
      W.Q=dmat,
      W.g=rep(1, length(indices)),
      Q.SL.library=SL.libraries,
      g.SL.library="SL.mean",
      k_split=model$k_split,
      verbose=T,
      save.sl.fit=T
    )

    # TODO: Need option for whether or not to perform a stratified fit
    # maybe in the adj_method option for homogeneous or heterogeneous?
    aipw$stratified_fit()
    muhat_ind <- do.call(cbind, aipw$obs_est[1:data$k])

    muhat[indices, ] <- muhat_ind

  }

  adjustment <- get.mutilde(model, data, muhat)
  mutilde <- adjustment$mutilde
  estimate <- colMeans(mutilde)

  # This is just a placeholder for the data frame, it's not
  # actually using this GLM model fit
  dmat <- get.dmat(data, model$adj_vars)
  df <- data.frame(response=data$response, treat=data$treat, dmat)
  glmod <- stats::glm(response ~ ., data=df)

  # Compute the asymptotic variance -- not specific to a partition
  asympt.variance <- vcov_car(model, data, glmod, mutilde)
  vcov_wt <- get.vcovHC(model$vcovHC, n=data$n, p=NULL)
  variance <- asympt.variance * vcov_wt / data$n

  result <- format_results(data$treat_levels, estimate, variance)

  return(
    structure(
      class="SLModelResult",
      list(result=result, varcov=variance, settings=model,
           data=data, mod=aipw, mu_a=mutilde)
    )
  )
}
