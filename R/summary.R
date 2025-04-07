#' Print linear model result
#'
#' @param x A LinModelResult object
#' @param ... Additional arguments
#' @export
#'
#' @returns Prints the treatment mean estimates (and variances) based on a linear working model,
#' along with the settings used. See \link[RobinCar:robincar_linear]{RobinCar::robincar_linear()}.
print.LinModelResult <- function(x, ...){
  descript(x)
  cat("\n\n")
  cat("Estimates:\n")
  print(x$result)
  cat("\nVariance-Covariance Matrix:\n")
  print(x$varcov)
}

#' Print glm model result
#'
#' @param x A GLMModelResult object
#' @param ... Additional arguments
#' @export
#'
#' @returns Prints the treatment mean estimates (and variances) based on a GLM working model,
#' along with the settings used. See \link[RobinCar:robincar_glm]{RobinCar::robincar_glm()}.
print.GLMModelResult <- function(x, ...){
  descript(x)
  cat("\n\n")
  cat("Estimates:\n")
  print(x$result)
  cat("\nVariance-Covariance Matrix:\n")
  print(x$varcov)
}

#' Print contrast result
#'
#' @param x A ContrastResult object
#' @param ... Additional arguments
#' @export
#'
#' @returns Prints estimates (and variances) of treatment contrasts based on a linear or GLM working model,
#' along with the settings used. See \link[RobinCar:robincar_contrast]{RobinCar::robincar_contrast()}
print.ContrastResult <- function(x, ...){
  if("DIFF" %in% class(x$settings)){
    c_type <- "linear contrast"
  } else if("RATIO" %in% class(x$settings)){
    c_type <- "ratio contrast"
  } else {
    c_type <- "custom contrast function"
  }
  output <- sprintf("Treatment group contrasts using %s", c_type)
  cat(output)
  cat("\n\n")
  cat("Contrasts:\n")
  print(x$result)
  cat("\nVariance-Covariance Matrix for Contrasts:\n")
  print(x$varcov)
}

#' Print TTE result
#'
#' @param x A TTEResult object
#' @param ... Additional arguments
#'
#' @importFrom data.table data.table setorder setnames .SD :=
#' @import data.table
#' @export
#'
#' @returns Prints results of time-to-event covariate adjusted analyses including covariate-adjusted (stratified) logrank,
#' robust Cox score, and covariate-adjusted hazard ratio. Prints summary statistics about number of observations and events, possibly by strata,
#' and the test statistics and/or estimates, and p-values. See \link[RobinCar:robincar_tte]{RobinCar::robincar_tte()} and \link[RobinCar:robincar_covhr]{RobinCar::robincar_covhr()}.
print.TTEResult <- function(x, ...){

  if(!is.null(x$data$covariate)){
    covariates <- colnames(x$data$covariate)
  } else {
    covariates <- c()
  }
  if(!is.null(x$data$car_strata)){
    car_strata <- colnames(x$data$car_strata)
  } else {
    car_strata <- c()
  }
  if(x$settings$method == "CSL"){
    if((x$settings$car_strata)){

    }
  }
  if((x$settings$method == "CL") | (x$settings$method == "CSL" & !x$settings$car_strata)){
    if((x$settings$adj_cov) | x$settings$adj_strata){
      if("TTEResultEst" %in% class(x)){
        cat("Performed covariate-adjusted cox hazard ratio estimation with covariates ",
            paste0(c(covariates, car_strata), sep=", "))
      } else {
        cat("Performed covariate-adjusted logrank test with covariates ",
            paste0(c(covariates, car_strata), sep=", "))
      }
    } else {
      if("TTEResultEst" %in% class(x)){
        cat("Performed cox hazard ratio estimation.")
      } else {
        cat("Performed logrank test.")
      }
    }
  } else if(x$settings$method == "CSL"){
    if((x$settings$adj_cov)){
      if("TTEResultEst" %in% class(x)){
        cat("Performed covariate-adjusted cox hazard ratio estimation with covariates ",
            paste0(c(covariates), sep=", "),
            " and stratifying by ", paste0(c(car_strata), sep=", "))
      } else {
        cat("Performed covariate-adjusted stratified logrank test with covariates ",
            paste0(c(covariates), sep=", "),
            " and stratifying by ", paste0(c(car_strata), sep=", "))
      }
    } else {
      if("TTEResultEst" %in% class(x)){
        cat("Performed stratified cox hazard ratio estimation stratifying by ",
            paste0(c(car_strata), sep=", "))
      } else {
        cat("Performed stratified logrank test stratifying by ",
            paste0(c(car_strata), sep=", "))
      }
    }
  } else if(x$settings$method == "coxscore"){
    cat("Performed coxscore test")
    if(x$settings$adj_cov){
      cat(", \n  adjusting for covariates", paste0(covariates, collapse=", "))
    }
    if(x$settings$adj_strata){
      cat(", \n  adjusting SE for car_strata", paste0(car_strata, collapse=", "))
    }
  }
  cat("\n------------------------------------\n")

  # Make visible binding for global variables
  # Recommended by data.table developers
  # These are names of columns of the data.table created below
  # and including them as NULL allows the R CMD CHECK to pass
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/

  df <- data.table::data.table(observed=x$data$event, treat=x$data$treat)
  df[, treat := as.character(treat)]

  df[, N := 1]
  id.cols <- c("treat")
  if(x$settings$car_strata){
    df$car_strata <- x$data$joint_strata
    id.cols <- c(id.cols, "car_strata")
    data.table::setorder(df, car_strata, treat)
  } else {
    data.table::setorder(df, treat)
  }
  txtitle <- "Treatment Group"

  summ <- df[, lapply(.SD, sum), by=id.cols, .SDcols=c("N", "observed")]
  summ[, name := paste0(x$data$treat_col, " = ", treat)]

  if(x$settings$car_strata){
    summ[, strata_col := paste0("car_strata = ", car_strata)]
    summ <- summ[, c("strata_col", "name", "N", "observed"), with=F]
    data.table::setnames(summ, c("Strata", "Treatment", "N.total", "N.events"))
  } else {
    summ <- summ[, c("name", "N", "observed")]
    data.table::setnames(summ, c("Treatment", "N.total", "N.events"))
  }

  print(summ)
  cat("\nReference arm is ", x$data$treat_col, "=", x$settings$ref_arm, "\n")
  if("TTEResultEst" %in% class(x)){
    stat <- x$result$theta_CL/x$result$se_theta_CL
    cat(
      "\nTest Stat:", x$result$theta_CL/x$result$se_theta_CL,
      "\n2-side p-value:", 2*(1-stats::pnorm(abs(stat))),
      "\nHR:", exp(x$result$theta_CL),
      "\nLog HR:", (x$result$theta_CL),
      "\nLog HR SE:", x$result$se_theta_CL
    )
  } else {
    cat(
      "\nScore function:", x$result$U,
      "\nStandard error:", x$result$se,
      "\nTest statistic:", x$result$statistic,
      "\n2-side p-value:", 2*(1-stats::pnorm(abs(x$result$statistic)))
    )
  }
}

#' Print calibration result
#'
#' @param x A GLMModel result. If you'd like to calibrate a linear
#'          adjustment, use `robincar_glm` instead of `robincar_linear`.
#' @param ... Additional arguments
#' @export
#'
#' @returns Prints the treatment mean estimates (and variances) based on a calibration on top of a
#' GLM working model, along with the settings used. See \link[RobinCar:robincar_calibrate]{RobinCar::robincar_calibrate()}.
print.CalibrationResult <- function(x, ...){
  if(x$joint){
    type <- "Joint calibration"
  } else {
    type <- "Linear calibration"
  }

  if(!is.null(x$add_x)){
    add_x_text <- paste0(", additionally adjusting for covariates",
                  paste0(x$add_x, sep=", "), ",")
  } else {
    add_x_text <- ""
  }
  output <- sprintf("%s%s on the following model result:",
                    type, add_x_text)

  cat(output)
  cat("\n")
  cat("--------------------------------------------\n")
  descript(x$original)
  cat("\n\nEstimates:\n")
  print(x$result)
  cat("\nVariance-Covariance Matrix:\n")
  print(x$varcov)
}


#' Print MH result
#'
#' @param x A ContrastResult object
#' @param ... Additional arguments
#' @export
#'
#' @returns Prints estimates (and variances) of treatment contrasts based on MH risk difference or ATE
print.MHResult <- function(x, ...){
  if("MH" == x$settings$estimand){
    est_type <- "Mantel-Haenszel risk difference"
  } else {
    est_type <- "ATE"
  }
  if("GR" == x$settings$ci_type){
    c_type <- "Greenland's estimator"
  } else if("mGR" == x$settings$ci_type){
    c_type <- "modified Greenland's estimator"
  } else{
    c_type <- "Sato's estimator"
  }
  output_title <- sprintf("Treatment group contrasts based on %s", est_type)
  output_estimand <- sprintf("Estimand: %s", est_type)
  output_strat <- sprintf("Stratified by %s", paste(x$settings$strata_cols, collapse=", "))
  output_ci <- sprintf("SE calculated via %s", c_type)
  cat(output_title)
  cat("\n")
  cat(output_estimand)
  cat("\n")
  cat(output_strat)
  cat("\n")
  cat(output_ci)
  cat("\n\n")
  cat("Contrasts:\n")
  print(x$result)
}
