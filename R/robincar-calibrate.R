#' Perform linear or joint calibration
#'
#' Uses linear or joint calibration to "calibrate" the estimates from a linear or GLM-type adjustment.
#' Linear calibration fits a linear model with treatment (and treatment-by-covariate interactions) and with the predicted \eqn{\hat{\bm \mu}(X_i) = (\hat{\mu}_1(X_i), \dots, \hat{\mu}_K(X_i))} as constructed covariates where \eqn{K} is the number of treatment groups;
#' joint calibration also includes \eqn{Z_i} the strata variables as covariates.
#'
#' @param result A GLMModelResult
#' @param joint If true, then performs joint calibration
#'              with the \eqn{\hat{\bm \mu}(X_i)} and strata \eqn{Z_i}
#'              to achieve universality and efficiency gain
#'              rather than just linear calibration that uses \eqn{\hat{\bm \mu}(X_i)}.
#' @param add_x Additional x to use in the calibration. Must have been in
#'              the original dataset that robincar_glm was called on.
#' @export
#'
#' @returns A result object that has the same structure as \link[RobinCar:robincar_glm]{RobinCar::robincar_glm()}, with the argument `result` included as "original" in the list.
robincar_calibrate <- function(result, joint=FALSE,
                               add_x=NULL){

  # Get the g-computation predictions
  mu_names <- paste0("mu_", result$data$treat_levels)
  mu_a <- data.frame(result$mu_a)
  colnames(mu_a) <- mu_names

  # Get new dataset to pass to new robincar_glm func
  newdat <- cbind(result$original_df, mu_a)

  # Add on more covariates
  if(!is.null(add_x)){
    for(cname in c(add_x)){
      if(cname %in% colnames(newdat)) next
      if(!cname %in% colnames(result$original_df)) .miss.covariate.calibrate()
      newdat <- cbind(newdat, result$original_df[cname])
    }
  }

  # Run robincar_glm using the heterogeneous
  # working model and identity link
  if(joint){
    # Use JOINT strata levels with formula notation
    strata <- colnames(result$data$car_strata)
    strata <- paste0(strata, collapse="*")
    covariate_cols <- c(mu_names, add_x, strata)
  } else {
    covariate_cols <- c(mu_names, add_x)
  }

  cal_result <- robincar_linear(
    df=newdat,
    response_col=result$data$response_col,
    treat_col=result$data$treat_col,
    covariate_cols=covariate_cols,
    car_strata_cols=colnames(result$data$car_strata),
    car_scheme=result$settings$car_scheme,
    adj_method="ANHECOVA"
  )
  # Designate as calibration result
  # for the print statement to work
  class(cal_result) <- c(
    "CalibrationResult",
    class(cal_result)
  )
  # Save additional settings
  cal_result[["original"]] <- result
  cal_result[["joint"]] <- joint
  cal_result[["add_x"]] <- add_x

  return(cal_result)
}
