# Gets the vcovHC weights for the sample size and number of parameters
get.vcovHC <- function(vcovHC, n, p){
  if(vcovHC == "HC0"){
    wgt <- 1
  } else if(vcovHC == "HC1") {
    wgt <- n / (n - p)
  } else if(vcovHC == "HC3") {
    wgt <- (n / (n - p)) ** 2
  } else {
    stop("Unrecognized vcovHC method.")
  }
  return(wgt)
}

# Get design matrix for the specified adjustment variables
# using the data stored.
get.dmat <- function(data, adj_vars){
  if(is.null(adj_vars)){
    dmat <- NULL
  } else if(adj_vars == "x"){
    dmat <- data$covariate
  } else if(adj_vars == "z"){
    dmat <- data$car_strata
  } else if(adj_vars == "joint_z"){
    dmat <- as.matrix(data$joint_strata, ncol=1)
    colnames(dmat) <- "joint_strata"
  } else if(adj_vars == "joint_z_x"){
    joint_strata <- data$joint_strata
    dmat <- cbind(data$covariate, joint_strata)
  } else if(adj_vars == "formula"){
    dmat <- data$formula_vars
  } else {
    stop(paste("Unrecognized adjustment variable type ", adj_vars))
  }
  return(dmat)
}

#' @importFrom stats model.matrix
.center.dmat <- function(dmat){
  modmat <- stats::model.matrix(~ 0 + ., data=data.frame(dmat))
  modmat <- t(t(modmat) - colMeans(modmat))
  return(modmat)
}

get.mutilde <- function(model, data, muhat){
  # Get response
  y <- data$response

  # Get treatment group indicators and outcome
  t_ids <- sapply(data$treat_levels, function(x) data$treat == x)

  if(!model$pu_joint_z){
    # Only look within treatment groups
    center_ids <- as.list(data.frame(t_ids))
    center_mus <- as.list(data.frame(muhat))
  } else {

    # Will create matrix for combination of car_strata and treatment groups
    sl <- data$joint_strata_levels
    tl <- data$treat_levels

    # Get indicators for joint car_strata group
    s_ids <- sapply(sl, function(x) data$joint_strata == x)

    # Get joint levels of treatment and car_strata groups
    center_ids <- as.list(data.frame(
      t_ids[, rep(1:ncol(t_ids), each=length(sl))] &
        s_ids[, rep(1:ncol(s_ids), times=length(tl))]
    ))

    # Repeat the mu columns for each car_strata
    center_mus <- as.list(data.frame(muhat[, rep(1:ncol(muhat), each=length(sl))]))
  }

  # Compute AIPW estimator by re-centering predictions
  # within treatment groups
  recenter <- function(u, i){
    cent <- sum(u[i])/sum(i) - sum(y[i])/sum(i)
    return(u - cent)
  }
  mutilde <- mapply(FUN=recenter,
                    u=center_mus,
                    i=center_ids)

  if(model$pu_joint_z){
    new_mutilde <- matrix(data=NA, nrow=data$n, ncol=ncol(muhat))
    for(s_col in 1:ncol(s_ids)){
      scol_seq <- seq(s_col, ncol(mutilde), by=length(sl))
      new_mutilde[s_ids[, s_col], ] <- mutilde[s_ids[, s_col], scol_seq]
    }
    mutilde <- new_mutilde
  }

  # Check prediction un-biasedness for the original muhat
  # g-computation just for warning/error reporting,
  # only up to level of accuracy specified by the user
  check.pu <- function(u, i) round(mean(u[i] - y[i]),
                                   digits=model$g_accuracy)

  # Calculate the group-specific residuals
  # where the groups are defined as treatment or treatment x car_strata above
  resid <- mapply(
    FUN=check.pu,
    u=center_mus,
    i=center_ids
  )

  dfree <- 0
  # Report warning or error messages, whatever
  # is passed through the model settings, if not prediction unbiased.
  if(!all(resid == 0)){
    if(!is.list(model$pu_funcs)){
      funcs <- list(model$pu_funcs)
    } else {
      funcs <- model$pu_funcs
    }
    for(func in funcs){
      func()
    }
    if(model$pu_joint_z){
      # Adjust the degrees of freedom for the model if prediction unbiasedness does not hold
      # and we end up having to adjust within levels of z
      # TODO: This is not going to work when there are covariates in the model...
      # Might need to do a more conservative adjustment
      # dfree <- (dfree - data$k) + data$k * length(data$joint_strata_levels)
      # dfree <- data$k * length(data$joint_strata_levels)
    }
  }

  return(list(mutilde=mutilde, df_adjust=dfree))
}
