# FUNCTIONS TO COMPUTE THE ASYMPTOTIC VARIANCE
# OF ESTIMATES FROM SIMPLE AND COVARIATE-ADAPTIVE RANDOMIZATION

# Gets the diagonal sandwich variance component
# for all linear models in the asymptotic variance formula.
#' @importFrom rlang .data
#' @importFrom dplyr mutate group_by summarize
vcov_sr_diag <- function(data, mod, residual=NULL){
  # Calculate the SD of the residuals from the model fit (or responses),
  # in order to compute sandwich variance -- this is
  # the asymptotic variance, not yet divided by n.
  if(is.null(residual)){
    residual <- stats::residuals(mod)
  }
  result <- dplyr::tibble(
    resid=residual,
    treat=data$treat
  ) %>%
    dplyr::group_by(.data$treat) %>%
    dplyr::summarize(se=stats::sd(.data$resid), .groups="drop") %>%
    dplyr::mutate(se=.data$se*sqrt(1/data$pie))

  return(diag(c(result$se)**2))
}

#' @importFrom rlang .data
#' @importFrom dplyr filter
get.erb <- function(model, data, mod, mu_hat=NULL){

  if(is.null(mu_hat)){
    mu_hat <- stats::predict(mod)
  }
  residual <- data$response - mu_hat

  # Calculate Omega Z under simple
  omegaz_sr <- omegaz.closure("simple")(data$pie)

  # Adjust this variance for Z

  # Calculate Omega Z under the covariate adaptive randomization scheme
  # Right now this will only be zero, but it's here for more generality
  # if we eventually include the urn design
  omegaz <- model$omegaz_func(data$pie)

  # Calculate the expectation of the residuals within each level of
  # car_strata variables Z
  dat <- dplyr::tibble(
    treat=data$treat,
    car_strata=data$joint_strata,
    resid=residual
  ) %>%
    dplyr::group_by(.data$treat, .data$car_strata) %>%
    dplyr::summarize(mean=mean(.data$resid), .groups="drop")

  # Calculate car_strata levels and proportions for
  # the outer expectation
  strata_levels <- data$joint_strata %>% levels
  strata_props <- data$joint_strata %>% table %>% proportions

  # Estimate R(B) by first getting the conditional expectation
  # vector for a particular car_strata (vector contains
  # all treatment groups), then dividing by the pi_t
  .get.cond.exp <- function(s) dat %>%
    dplyr::filter(.data$car_strata==s) %>%
    dplyr::arrange(.data$treat) %>%
    dplyr::pull(mean)

  .get.rb <- function(s) diag(.get.cond.exp(s) / c(data$pie))

  # Get the R(B) matrix for all car_strata levels
  rb_z <- lapply(strata_levels, .get.rb)

  # Compute the R(B)[Omega_{SR} - Omega_{Z_i}]R(B) | Z_i
  # for each Z_i
  rb_omega_rb_z <- lapply(rb_z, function(x) x %*% (omegaz_sr - omegaz) %*% x)

  # Compute the outer expectation
  ERB <- mapply(function(x, y) x*y, x=rb_omega_rb_z, y=strata_props, SIMPLIFY=FALSE)
  ERB <- Reduce("+", ERB)
  return(ERB)
}

# Gets AIPW asymptotic variance under simple randomization
#' @importFrom stats cov
#' @importFrom stats var
vcov_car <- function(model, data, mod, mutilde){

  # Get predictions for observed treatment group
  preds <- matrix(nrow=data$n, ncol=1)
  for(t_id in 1:length(data$treat_levels)){
    t_group <- data$treat == data$treat_levels[t_id]
    preds[t_group] <- mutilde[, t_id][t_group]
  }

  # Compute residuals for the diagonal portion
  residual <- data$response - preds

  # Get covariance between observed Y and predicted \mu counterfactuals
  get.cov.Ya <- function(a){
    t_group <- data$treat == a
    cv <- stats::cov(data$response[t_group], mutilde[t_group, ])
    return(cv)
  }
  # Covariance matrix between Y and \mu
  cov_Ymu <- t(sapply(data$treat_levels, get.cov.Ya))

  # Get specific variance type for simple randomization
  if(model$variance_type == 1){

    # We are avoiding the situation where we have issues in estimating
    # the variance when many (or all) of the Y_a are zero in one group
    var_mutilde <- stats::var(mutilde)
    diag_mutilde <- diag(diag(var_mutilde))
    diag_covYmu <- diag(diag(cov_Ymu))
    diag_pi <- diag(1/c(data$pie))

    # Formula for variance calculation, doing a decomposition of variance
    # rather than calculating variance of residual
    diagmat <- vcov_sr_diag(data, mod, residual=data$response) +
      (diag_mutilde - 2 * diag_covYmu) * diag_pi

    # Sum of terms to compute simple randomization variance
    v <- diagmat + cov_Ymu + t(cov_Ymu) - stats::var(mutilde)

  } else if(model$variance_type == 2) {

    # This is the simplest type which just calculates
    # the variance of the residuals directly
    diagmat <- vcov_sr_diag(data, mod, residual=residual)
    v <- diagmat + cov_Ymu + t(cov_Ymu) - stats::var(mutilde)

  } else if(model$variance_type == 3){

    # This variance type is based on @Eureeca's work
    # which has a SR variance matrix which is for sure positive definite
    # in some special cases
    #
    # Biggest difference is that it uses observations in group a and b
    # to calculate the covariance between \mu_a(X) and \mu_b(X), rather than the
    # whole sample.
    #
    # This is helpful when there are large chance imbalances
    # in the distribution of X across the treatment groups, since
    # the covariance of Y_a and \mu_a(X) can only be calculated within group a.

    K <- length(data$treat_levels)
    v <- matrix(data=NA, nrow=K, ncol=K)

    for(a in 1:K){
      for(b in a:K){

        # Diagonal entries in the vcov matrix
        if(a == b){

          t_group <- data$treat == data$treat_levels[[a]]

          term1 <- mean(as.numeric(t_group) * residual^2) / data$pie[[a]]^2
          term2 <- cov(residual[which(t_group)], mutilde[which(t_group), a])
          term3 <- var(mutilde[which(t_group), a])

          v[a, a] <- term1 + 2 * term2 + term3

        # Off-diagonal entries in the vcov matrix
        } else {

          a_group <- data$treat == data$treat_levels[[a]]
          b_group <- data$treat == data$treat_levels[[b]]

          term1 <- cov(data$response[which(a_group)], mutilde[which(a_group), b])
          term2 <- cov(data$response[which(b_group)], mutilde[which(b_group), a])
          term3 <- cov(mutilde[which(a_group), c(a, b)])[1, 2]
          term4 <- cov(mutilde[which(b_group), c(a, b)])[1, 2]

          v[a, b] <- v[b, a] <- term1 + term2 - (term3 + term4) / 2

        }
      }
    }

  } else {
    stop("Unrecognized variance type.")
  }

  # Adjust for Z if necessary
  if(!is.null(model$omegaz_func)) v <- v - get.erb(model, data, mod, mu_hat=preds)

  return(v)
}
