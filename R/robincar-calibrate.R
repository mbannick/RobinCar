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
  newdat <- data.frame(
    response=result$data$response,
    treat=result$data$treat,
    mu_a
  )

  # Add on covariates and car_strata from the original
  # robincar_glm function call
  if(!is.null(result$data$covariate)){
    newdat <- cbind(newdat, result$data$covariate)
  }
  if(!is.null(add_x)){
    for(cname in c(add_x)){
      if(cname %in% colnames(newdat)) next
      if(!cname %in% colnames(result$original_df)) .miss.covariate.calibrate()
      newdat <- cbind(newdat, result$original_df[cname])
    }
  }
  if(!is.null(result$data$car_strata)){
    for(cname in colnames(result$data$car_strata)){
      if(cname %in% colnames(newdat)) next
      newdat <- cbind(newdat, result$data$car_strata[cname])
    }
  }

  # Run robincar_glm using the heterogeneous
  # working model and identity link
  cal_result <- robincar_glm(
    df=newdat,
    response_col="response",
    treat_col="treat",
    covariate_cols=c(mu_names, add_x),
    car_strata_cols=colnames(result$data$car_strata),
    covariate_to_include_strata=joint,
    car_scheme=result$settings$car_scheme,
    adj_method="heterogeneous",
    g_family=gaussian(link="identity"),
    g_accuracy=result$settings$g_accuracy
    # vcovHC=vcovHC
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
