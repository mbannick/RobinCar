#' Covariate adjustment using linear working model
#'
#' Estimate treatment-group-specific response means and (optionally)
#' treatment group contrasts using a linear working model for continuous outcomes.
#'
#' @param df A data.frame with the required columns
#' @param treat_col Name of column in df with treatment variable
#' @param response_col Name of the column in df with response variable
#' @param strata_cols Names of columns in df with strata variables
#' @param covariate_cols Names of columns in df with covariate variables
#' @param car_scheme Name of the type of covariate-adaptive randomization scheme. One of: "simple", "pocock-simon", "biased-coin", "permuted-block".
#' @param adj_method Name of linear adjustment method to use. One of: "ANOVA", "ANCOVA", "ANHECOVA".
#' @param vcovHC Type of heteroskedasticity-consistent variance estimates. One of: "HC0", "HC1", "HC3".
#' @param covariate_to_include_strata Whether to include strata variables in covariate adjustment. Defaults to F for ANOVA and ANCOVA; defaults to T for ANHECOVA. User may override by passing in this argument.
#' @param contrast An optional function to specify a desired contrast
#'
#' @import dplyr
#' @export
robincar_linear <- function(df,
                            treat_col, response_col, strata_cols=NULL, covariate_cols=NULL,
                            car_scheme="simple", adj_method="ANOVA", vcovHC="HC0",
                            covariate_to_include_strata=NULL,
                            conf_level=0.95,
                            contrast_h=NULL, contrast_dh=NULL){

  .check.car_scheme(car_scheme)
  .check.adj_method.linear(adj_method)
  .check.vcovHC(vcovHC)

  # Create data object and validate
  data <- .make.data(
    df=df, classname="RoboDataLinear",
    treat_col=treat_col,
    response_col=response_col,
    strata_cols=strata_cols,
    covariate_cols=covariate_cols
  )
  validate(data)

  # Create model object
  model <- .make.model(
    data=data,
    adj_method=adj_method,
    car_scheme=car_scheme,
    vcovHC=vcovHC,
    covariate_to_include_strata=covariate_to_include_strata
  )

  # Perform linear adjustment
  result <- adjust(model, data)

  # Create transformation object
  if(!is.null(contrast_h)){
    c_result <- robincar_contrast(
      result=result,
      contrast_h=contrast_h,
      contrast_dh=contrast_dh
    )
    result <- list(
      main=result, contrast=c_result
    )
  }

  print(result)
  return(result)
}
