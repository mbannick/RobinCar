#' Generate description for glm model result
#'
#' @param x A GLMModelResult object
#' @param ... Additional arguments
descript.GLMModelResult <- function(x, ...){
  output <- c()
  if("AIPW" %in% class(x$settings)){
    etype <- "aipw-type"
  } else {
    etype <- "g-computation-type"
  }
  if(is.character(x$settings$g_family)){
    family <- get(x$settings$g_family)()
  } else if(is.function(x$settings$g_family)){
    family <- (x$settings$g_family)()
  } else {
    family <- x$settings$g_family
  }
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
              paste0(names(x$data$strata), collapse=", "))
    )
  }
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
        sprintf("\nand adjusted variance-covariance matrix for randomization strata %s \nconsistent with the %s design.",
              paste0(strata, collapse=", "),
              x$settings$car_scheme)
      )
    }
  }
  for(o in output){
    cat(o)
  }
}
