
#' @importFrom dplyr group_by mutate ungroup
get.ordered.data.est <- function(df, lin_predictors){

  data <- df %>% dplyr::mutate(lin_pred=lin_predictors)

  data <- data %>%
    dplyr::group_by(.data$car_strata) %>%
    dplyr::mutate(
      s0_seq          = exp(.data$lin_pred)
    ) %>% dplyr::mutate(
      mu_num          = .data$s0_seq * .data$Y1,
      mu_denom        = (.data$s0_seq * .data$Y1 + .data$Y0)
    ) %>% dplyr::mutate(
      mu_t            = .data$mu_num / .data$mu_denom,
      score_i         = .data$event * (.data$trt1 - .data$mu_t),
      O.hat1          = .data$event * .data$Y0 / .data$mu_denom -
        .data$s0_seq * cumsum(.data$event * .data$Y0 / .data$mu_denom**2),
      O.hat0          = .data$event * .data$s0_seq * .data$Y1 / .data$mu_denom -
        .data$s0_seq * cumsum(.data$event * .data$Y1 / .data$mu_denom**2)
    ) %>% dplyr::mutate(
      O.hat           = .data$O.hat0 * (.data$trt0 == 1) +
        .data$O.hat1 * (.data$trt1 == 1)
    ) %>% dplyr::ungroup()

  return(data)

}

# Define score of theta function
score.theta <- function(theta, df){
  ordered <- get.ordered.data.est(lin_predictors=theta, df=df)
  score <- mean(ordered$score_i)
  return(score)
}

#' @importFrom dplyr mutate group_by summarise arrange
#' @importFrom stats uniroot
#' @exportS3Method
adjust.CovHR <- function(model, data, ...){

  # Creates data
  df <- create.tte.df(model, data)

  # Risk set and failure time times # THIS IS WHERE THE TEST FAILS
  df <- process.tte.df(df, ref_arm=model$ref_arm)
  df_process <- df # save unmodified for later

  # Get un-adjusted estimate
  score.unadj <- function(theta) score.theta(theta, df=df_process)
  thetaLhat <- stats::uniroot(score.unadj, interval=model$interval)$root

  # Perform covariate adjustment
  df <- get.ordered.data.est(df=df, lin_predictors=thetaLhat)

  # Perform covariate adjustment
  df <- get.tte.adjustment(df, model, data)

  df <- df %>%
    dplyr::mutate(
      uu_cl = .data$trt1 * .data$adjust1 - .data$trt0 * .data$adjust0,
    )
  cov_adjust <- mean(df$uu_cl)

  # Get adjusted estimate
  score.adj <- function(theta) score.theta(theta, df=df_process) - cov_adjust
  thetaCLhat <- stats::uniroot(score.adj, interval=model$interval)$root

  # Use thetaLhat rather than thetaCLhat below Both result in consistent
  # variance estimators but using thetaLhat is expected behavior.
  df <- df %>% dplyr::mutate(
    ssig_l = .data$event * exp(thetaLhat) * .data$Y0 * .data$Y1 /
      (exp(thetaLhat) * .data$Y1 + .data$Y0)**2
  )

  # Summarize by car_strata (if CL, then single car_strata)
  ss <- df %>%
    dplyr::filter(!is.na(.data$uu_cl)) %>%
    dplyr::group_by(.data$car_strata) %>%
    dplyr::summarise(
      var_adj = model$p_trt * (1 - model$p_trt) * unique(.data$bsigb) * dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$car_strata)

  # Final quantities for the C(S)L statistic
  sig2_L  <- mean(df$ssig_l)
  sig2_CL <- sig2_L - sum(ss$var_adj) / data$n
  var_est <- sig2_CL / sig2_L**2

  se_theta_L <- sqrt(1 / sig2_L / data$n)
  se_theta_CL <- sqrt(var_est / data$n)

  result <- list(
    strata_sum=ss,
    theta_L=thetaLhat,
    se_theta_L=se_theta_L,
    theta_CL=thetaCLhat,
    se_theta_CL=se_theta_CL
  )

  return(
    structure(
      class=c("TTEResultEst", "TTEResult"),
      list(result=result, settings=model, data=data)
    )
  )
}
