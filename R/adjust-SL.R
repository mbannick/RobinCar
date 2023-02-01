
#' Perform GLM adjustment, based on the classes
#' of the model. Will perform adjustment based on the linear
#' model type of `model` and also do G-computation or AIPW
#' based on the second model type of `model`.
adjust.SLModel <- function(model, data){
  browser()

  estimates <- list()

  # TODO: Estimate \pi_a
  pi_a <- ...

  for(i in 1:model$k_split){ # k splits
    # Get covariates for outcome regression -- everything but the k'th partition
    dmat <- get.dmat(data[!data$k], model$adj_vars)
    # TODO: Use superlearner to estimate $\mu_k(X_i)$ for each observation in k
    # and treatment group
    sl_muk <- ...

    # TODO: The mutilde function needs to be modified so that
    # it uses \hat{\pi}_a on the whole sample. Should pass in a value for
    # \hat{\pi}_a here so that we can override the default behavior of get.mutilde
    mutilde <- get.mutilde(model, data, muhat=sl_muk, pi=pi_a)
    estimates[[i]] <- colMeans(mutilde)
  }

  # TODO: Then average the estimates of \check{\theta}_k to get \hat{\theta}

  # Compute the asymptotic variance -- not specific to a partition
  asympt.variance <- vcov_car(model, data, glmod, mutilde)
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
