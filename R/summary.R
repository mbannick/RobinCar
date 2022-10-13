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
                  paste0(add_x, sep=", "), ",")
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
