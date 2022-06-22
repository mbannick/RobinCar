
#' @export
print.LinModelResult <- function(res){
  output <- c()
  output <- c(
    output,
    sprintf("Treatment group mean estimates using fit from an %s model",
            class(res$settings)[2])
  )
  if(class(res$settings)[2] != "ANOVA"){
    k <- length(res$data$treat_levels)
    mat <- .get.dmat(res$data, res$settings$adj_vars)
    name <- colnames(mat)
    cov_name <- name[name != "joint_strata"]
    if("joint_strata" %in% name){
      cov_name <- c(
        cov_name,
        sprintf("joint levels of %s",
                paste0(names(res$data$strata), collapse=", "))
      )
    }
    output <- c(
      output,
      sprintf("\n  using adjustment variables: %s",
              paste0(cov_name, collapse=", "))
    )
    if(class(res$settings)[2] == "ANHECOVA"){
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
            res$settings$vcovHC)
  )
  if(!is.null(res$settings$omegaz_func)){
    if(all(res$settings$omegaz_func(res$data$pie) == 0)){
      strata <- colnames(res$data$strata)
      output <- c(
        output,
        sprintf("\n  and adjusted variance-covariance matrix for randomization strata consistent with the
              %s design: %s.",
              res$settings$car_scheme,
              paste0(strata, collapse=", "))
      )
    }
  }
  for(o in output){
    cat(o)
  }
  cat("\n\n")
  cat("Estimates:\n")
  print(res$result)
  cat("\nVariance-Covariance Matrix:\n")
  print(res$varcov)
}

#' @export
print.ContrastResult <- function(res){
  if("DIFF" %in% class(res$settings)){
    c_type <- "linear contrast"
  } else if("RATIO" %in% class(res$settings)){
    c_type <- "ratio contrast"
  } else {
    c_type <- "custom contrast function"
  }
  output <- sprintf("Treatment group contrasts using %s", c_type)
  cat(output)
  cat("\n\n")
  cat("Contrasts:\n")
  print(res$result)
  cat("\nVariance-Covariance Matrix for Contrasts:\n")
  print(res$varcov)
}

