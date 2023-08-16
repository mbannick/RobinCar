#' Print linear model result
#'
#' @param x A LinModelResult object
#' @param ... Additional arguments
#' @export
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
#' @importFrom data.table data.table setorder setnames
#' @export
print.TTEResult <- function(x, ...){

  if(!is.null(x$data$covariate)){
    covariates <- colnames(x$data$covariate)
  } else {
    covariates <- c()
  }
  if(!is.null(x$data$strata)){
    strata <- colnames(x$data$strata)
  } else {
    strata <- c()
  }
  if(x$settings$method == "CSL"){
    if((x$settings$car_strata)){

    }
  }
  if((x$settings$method == "CL") | (x$settings$method == "CSL" & !x$settings$car_strata)){
    if((x$settings$adj_cov) | x$settings$adj_strata){
      if("TTEResultEst" %in% class(x)){
        cat("Performed covariate-adjusted cox hazard ratio estimation with covariates ",
            paste0(c(covariates, strata), sep=", "))
      } else {
        cat("Performed covariate-adjusted logrank test with covariates ",
            paste0(c(covariates, strata), sep=", "))
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
            " and stratifying by ", paste0(c(strata), sep=", "))
      } else {
        cat("Performed covariate-adjusted stratified logrank test with covariates ",
            paste0(c(covariates), sep=", "),
            " and stratifying by ", paste0(c(strata), sep=", "))
      }
    } else {
      if("TTEResultEst" %in% class(x)){
        cat("Performed stratified cox hazard ratio estimation stratifying by ",
            paste0(c(strata), sep=", "))
      } else {
        cat("Performed stratified logrank test stratifying by ",
            paste0(c(strata), sep=", "))
      }
    }
  } else if(x$settings$method == "coxscore"){
    cat("Performed coxscore test")
    if(x$settings$adj_cov){
      cat(", \n  adjusting for covariates", paste0(covariates, collapse=", "))
    }
    if(x$settings$adj_strata){
      cat(", \n  adjusting SE for strata", paste0(strata, collapse=", "))
    }
  }
  cat("\n------------------------------------\n")

  df <- data.table::data.table(observed=x$data$event, treat=x$data$treat)
  df[, treat := as.character(treat)]

  df[, N := 1]
  id.cols <- c("treat")
  data.table::setorder(df, treat)
  txtitle <- "Treatment Group"

  if(x$settings$car_strata){
    df$strata <- x$data$joint_strata
    id.cols <- c(id.cols, "strata")
    data.table::setorder(df, strata, treat)
  }
  summ <- df[, lapply(.SD, sum), by=id.cols, .SDcols=c("N", "observed")]
  summ[, name := paste0(x$data$treat_col, " = ", treat)]

  if(x$settings$car_strata){
    summ[, strata_col := paste0("strata = ", strata)]
    summ <- summ[, .(strata_col, name, N, observed)]
    data.table::setnames(summ, c("Strata", "Treatment", "N.total", "N.events"))
  } else {
    summ <- summ[, .(name, N, observed)]
    data.table::setnames(summ, c("Treatment", "N.total", "N.events"))
  }

  print(summ)
  cat("\nReference arm is ", x$data$treat_col, "=", x$settings$ref_arm, "\n")
  if("TTEResultEst" %in% class(x)){
    stat <- x$result$theta_CL/x$result$se_theta_CL
    cat(
      "\nTest Stat:", x$result$theta_CL/x$result$se_theta_CL,
      "\n2-side p-value:", 2*(1-pnorm(abs(stat))),
      "\nHR:", x$result$theta_CL,
      "\nLog HR:", log(x$result$theta_CL),
      "\nHR SE:", x$result$se_theta_CL
    )
  } else {
    cat(
      "\nScore function:", x$result$U,
      "\nStandard error:", x$result$se,
      "\nTest statistic:", x$result$statistic,
      "\n2-side p-value:", 2*(1-pnorm(abs(x$result$statistic)))
    )
  }
}

#' Print calibration result
#'
#' @param x A GLMModel result. If you'd like to calibrate a linear
#'          adjustment, use `robincar_glm` instead of `robincar_linear`.
#' @param ... Additional arguments
#' @export
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
