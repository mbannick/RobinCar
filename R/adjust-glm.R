
predictions <- function(model, data, mod){
  UseMethod("predictions", model)
}

predictions.GLMModel <- function(model, data, mod){
  df <- data.frame(
    treat=data$treat,
    response=data$response
  )
  dmat <- .get.dmat(data, model$adj_vars)
  df <- cbind(df, dmat)
  preds <- stats::predict(mod, newdata=df, type="response")
  return(preds)
}

.get.muhat <- function(model, data, mod){

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

.get.mutilde <- function(model, data, mod){

  # Get response
  y <- data$response

  # Get g-computation predictions
  muhat <- .get.muhat(model, data, mod)

  # Get treatment group indicators and outcome
  t_ids <- sapply(data$treat_levels, function(x) data$treat == x)

  if(!model$pu_joint_z){
    # Only look within treatment groups
    center_ids <- as.list(data.frame(t_ids))
    center_mus <- as.list(data.frame(muhat))
  } else {
    # Will create matrix for combination of strata and treatment groups
    sl <- data$joint_strata_levels
    tl <- data$treat_levels

    # Get indicators for joint strata group
    s_ids <- sapply(sl, function(x) data$joint_strata == x)

    # Get joint levels of treatment and strata groups
    center_ids <- as.list(data.frame(
      t_ids[, rep(1:ncol(t_ids), each=length(sl))] &
        s_ids[, rep(1:ncol(s_ids), times=length(tl))]
    ))

    # Repeat the mu columns for each strata
    center_mus <- as.list(data.frame(muhat[, rep(1:ncol(muhat), each=length(sl))]))
  }

  # Compute AIPW estimator by re-centering predictions
  # within treatment groups
  recenter <- function(u, i) u - sum(u[i])/sum(i) + sum(y[i])/sum(i)
  mutilde <- mapply(FUN=recenter,
                    u=center_mus,
                    i=center_ids)

  # Check prediction un-biasedness for the original muhat
  # g-computation just for warning/error reporting,
  # only up to level of accuracy specified by the user
  check.pu <- function(u, i) round(mean(u[i] - y[i]),
                                   digits=model$g_accuracy)

  # Calculate the group-specific residuals
  # where the groups are defined as treatment or treatment x strata above
  resid <- mapply(
    FUN=check.pu,
    u=center_mus,
    i=center_ids
  )

  # Report warning or error messages, whatever
  # is passed through the model settings, if not prediction unbiased.
  if(!all(resid == 0)){
    for(func in list(model$pu_funcs)){
      func()
    }
  }

  return(mutilde)
}

# Perform GLM adjustment, based on the classes
# of the model. Will perform adjustment based on the linear
# model type of `model` and also do G-computation or AIPW
# based on the second model type of `model`.
adjust.GLMModel <- function(model, data){

  # Get the generalized linear model and estimates by AIPW
  glmod <- linmod(model, data, family=model$g_family, center=FALSE)

  # Get mutilde from the GLM model, then estimate the treatment means by
  # taking the mean over all of the potential outcomes
  mutilde <- .get.mutilde(model, data, glmod)
  estimate <- colMeans(mutilde)

  # Compute the asymptotic variance
  asympt.variance <- vcov_car(model, data, glmod, mutilde)
  vcov_wt <- .get.vcovHC(model$vcovHC, n=data$n, p=glmod$rank)
  variance <- asympt.variance * vcov_wt / data$n

  result <- format.results(data$treat_levels, estimate, variance)

  return(list(result=result, varcov=variance, mutilde=mutilde, settings=model))
}
