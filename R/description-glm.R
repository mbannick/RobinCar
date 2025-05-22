
#' @exportS3Method
descript.GLMModelResult <- function(x, ...){
  output <- c()
  if("AIPW" %in% class(x$settings)){
    etype <- "aipw-type"
  } else {
    etype <- "g-computation-type"
  }

  if(is.character(x$settings$g_family)){
    if(x$settings$g_family == "nb"){
      family <- list(
        family="negative binomial with unknown dispersion",
        link="log"
      )
    } else {
      family <- get(x$settings$g_family)()
    }
  } else if(is.function(x$settings$g_family)){
    family <- (x$settings$g_family)()
  } else {
    family <- x$settings$g_family
  }

  if(!is.null(x$data$formula)){

    if(is.character(x$settings$g_family)){
      if(x$settings$g_family == "nb"){
        form <- x$data$formula
      } else {
        form <- x$mod$formula
      }
    } else {
      form <- x$mod$formula
    }

    output <- c(
      output,
      sprintf("Treatment group mean estimates from a GLM working model of family %s and link %s using formula: \n",
              family$family, family$link),
      form
    )
  } else {
    output <- c(
      output,
      sprintf("Treatment group mean estimates using a %s estimator",
              etype),
      sprintf("\nfrom a GLM working model of family %s and link %s",
              family$family, family$link)
    )
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
    if(!is.null(x$settings$adj_vars)){
      output <- c(
        output,
        sprintf("\nusing adjustment variables: %s",
                paste0(cov_name, collapse=", "))
      )
      if("ANHECOVA" %in% class(x$settings)){
        output <- c(
          output,
          "\nand their interactions with treatment."
        )
      } else if("ANCOVA" %in% class(x$settings)){
        output <- c(output, ".")
      } else {
        stop("Error with the type of model.")
      }
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
        sprintf("\nand adjusted variance-covariance matrix for randomization car_strata %s \nconsistent with the %s design.",
              paste0(car_strata, collapse=", "),
              x$settings$car_scheme)
      )
    }
  }
  for(o in output){
    if("formula" %in% class(o)){
      print(o, showEnv=FALSE)
    } else {
      cat(o)
    }
  }
}
