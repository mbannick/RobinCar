
#' @exportS3Method
descript.LinModelResult <- function(x, ...){
  output <- c()
  output <- c(
    output,
    sprintf("Treatment group mean estimates using fit from an %s model",
            class(x$settings)[2])
  )
  if(class(x$settings)[2] != "ANOVA"){
    k <- length(x$data$treat_levels)
    mat <- get.dmat(x$data, x$settings$adj_vars)
    name <- colnames(mat)
    cov_name <- name[name != "joint_strata"]
    if("joint_strata" %in% name){
      cov_name <- c(
        cov_name,
        sprintf("joint levels of %s",
                paste0(names(x$data$car_strata), collapse=", "))
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
      car_strata <- colnames(x$data$car_strata)
      output <- c(
        output,
        sprintf("\nand adjusted variance-covariance matrix for randomization car_strata %s \nconsistent with the
              %s design.",
              paste0(car_strata, collapse=", "),
              x$settings$car_scheme)
      )
    }
  }
  for(o in output){
    cat(o)
  }
}
