#' Run a linear or joint calibration on
#' a GLMModelResult object.
#'
#' @param result A GLMModelResult
#' @param joint If true, then performs joint calibration
#'              with the mu and Z
#'              to achieve universality and efficiency gain
#'              rather than just linear calibration that uses mu.
#' @param add_x Additional x to use in the calibration. Must have been in
#'              the original dataset that robincar_glm was called on.
#' @param vcovHC Which type of heteroskedasticity-consistent
#'               standard error estimates to use.
#' @export
robincar_calibrate <- function(result, joint=FALSE,
                               add_x=NULL, vcovHC="HC0"){

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

  # Add on covariates and strata from the original
  # robincar_glm function call
  if(!is.null(result$data$covariate)){
    newdat <- cbind(newdat, result$data$covariate)
  }
  if(!is.null(add_x)){
    if(!add_x %in% colnames(data)) .miss.covariate.calibrate()
    newdat <- cbind(newdat, result$original_df[add_x])
  }
  if(!is.null(result$data$strata)){
    newdat <- cbind(newdat, result$data$strata)
  }

  # Run robincar_glm using the heterogeneous
  # working model and identity link
  cal_result <- robincar_glm(
    df=newdat,
    response_col="response",
    treat_col="treat",
    covariate_cols=c(mu_names, add_x),
    strata_cols=colnames(result$data$strata),
    covariate_to_include_strata=joint,
    car_scheme=result$settings$car_scheme,
    adj_method="heterogeneous",
    g_family=gaussian(link="identity"),
    g_accuracy=result$settings$g_accuracy,
    vcovHC=vcovHC
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
