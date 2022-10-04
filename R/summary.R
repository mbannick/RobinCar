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

#' Print glm model result
#'
#' @param x A GLMModelResult object
#' @param ... Additional arguments
#' @export
print.GLMModelResult <- function(x, ...){
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
    sprintf("\n  from a GLM working model of family %s and link %s",
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
    sprintf("\n  using adjustment variables: %s",
            paste0(cov_name, collapse=", "))
  )
  if(class(x$settings)[2] == "heterogeneous"){
    output <- c(
      output,
      "\n  and their interactions with treatment."
    )
  } else {
    output <- c(output, ".")
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

#' Print TTE result
#'
#' @param x A TTEResult object
#' @param ... Additional arguments
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
      cat("Performed covariate-adjusted logrank test with covariates ",
          paste0(c(covariates, strata), sep=", "))
    } else {
      cat("Performed logrank test.")
    }
  } else if(x$settings$method == "CSL"){
    if((x$settings$adj_cov)){
      cat("Performed covariate-adjusted stratified logrank test with covariates ",
          paste0(c(covariates), sep=", "),
          " and stratifying by ", paste0(c(strata), sep=", "))
    } else {
      cat("Performed stratified logrank test stratifying by ",
          paste0(c(strata), sep=", "))
    }
  } else if(x$settings$method == "coxscore"){
    cat("Performed coxscore test ")
    if(x$settings$adj_cov){
      cat("adjusting for covariates ", paste0(covariates, sep=", "))
    }
    if(x$settings$adj_strata){
      cat("adjusting SE for strata ", paste0(strata, sep=", "))
    }
  }
  cat("\n----------------------------\n\n")

  df <- data.table(observed=x$data$event, treat=x$data$treat)
  df[, treat := as.character(treat)]

  df[, N := 1]
  id.cols <- c("treat")
  setorder(df, treat)
  txtitle <- "Treatment Group"

  if(x$settings$car_strata){
    df$strata <- x$data$joint_strata
    id.cols <- c(id.cols, "strata")
    setorder(df, strata, treat)
  }
  summ <- df[, lapply(.SD, sum), by=id.cols, .SDcols=c("N", "observed")]
  summ[, name := paste0(x$data$treat_col, " = ", treat)]

  if(x$settings$car_strata){
    summ[, strata_col := paste0("strata = ", strata)]
    summ <- summ[, .(strata_col, name, N, observed)]
    setnames(summ, c("Strata", "Treatment", "N.total", "N.events"))
  } else {
    summ <- summ[, .(name, N, observed)]
    setnames(summ, c("Treatment", "N.total", "N.events"))
  }

  print(summ)
  cat("\nReference arm is ", x$data$treat_col, "=", x$settings$ref_arm, "\n")
  cat("\nScore function:", x$result$U,
      "\nStandard error:", x$result$se,
      "\nTest statistic:", x$result$statistic,
      "\n2-side p-value:", 2*(1-pnorm(abs(x$result$statistic))))
}



