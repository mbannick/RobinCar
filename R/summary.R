#' Print linear model result
#'
#' @param x A LinModelResult object
#' @param ... Additional arguments
#' @export
print.LinModelResult <- function(x, ...){
  output <- c()
  output <- c(
    output,
    sprintf("Treatment group mean estimates using fit from an %s model",
            class(x$settings)[2])
  )
  if(class(x$settings)[2] != "ANOVA"){
    k <- length(x$data$treat_levels)
    mat <- .get.dmat(x$data, x$settings$adj_vars)
    name <- colnames(mat)
    cov_name <- name[name != "joint_strata"]
    if("joint_strata" %in% name){
      cov_name <- c(
        cov_name,
        sprintf("joint levels of %s",
                paste0(names(x$data$strata), collapse=", "))
      )
    }
    output <- c(
      output,
      sprintf("\n  using adjustment variables: %s",
              paste0(cov_name, collapse=", "))
    )
    if(class(x$settings)[2] == "ANHECOVA"){
      output <- c(
        output,
        "\n  and their interactions with treatment."
      )
    } else {
      output <- c(output, ".")
    }
  }
  output <- c(
    output,
    sprintf("\n\nUsed %s-type of heteroskedasticity-consistent variance estimates ",
            x$settings$vcovHC)
  )
  if(!is.null(x$settings$omegaz_func)){
    if(all(x$settings$omegaz_func(x$data$pie) == 0)){
      strata <- colnames(x$data$strata)
      output <- c(
        output,
        sprintf("\n  and adjusted variance-covariance matrix for randomization strata consistent with the
              %s design: %s.",
              x$settings$car_scheme,
              paste0(strata, collapse=", "))
      )
    }
  }
  for(o in output){
    cat(o)
  }
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

