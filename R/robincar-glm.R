#' Covariate adjustment using generalized linear working model
#'
#' WARNING: this function is still under development.
#' Estimate treatment-group-specific response means and (optionally)
#' treatment group contrasts using a generalized linear working model.
#'
#'
#' @param df A data.frame with the required columns
#' @param treat_col Name of column in df with treatment variable
#' @param response_col Name of the column in df with response variable
#' @param strata_cols Names of columns in df with strata variables
#' @param covariate_cols Names of columns in df with covariate variables
#' @param car_scheme Name of the type of covariate-adaptive randomization scheme. One of: "simple", "pocock-simon", "biased-coin", "permuted-block".
#' @param adj_method Name of adjustment method to use, one of "heterogeneous" (interaction model) or "homogeneous"
#' @param vcovHC Type of heteroskedasticity-consistent variance estimates. One of: "HC0", "HC1", "HC3".
#' @param covariate_to_include_strata Whether to include strata variables in covariate adjustment. Defaults to F for ANOVA and ANCOVA; defaults to T for ANHECOVA. User may override by passing in this argument.
#' @param contrast_h An optional function to specify a desired contrast
#' @param contrast_dh An optional jacobian function for the contrast (otherwise use numerical derivative)
#' @param g_family Family that would be supplied to glm(...), e.g., binomial. If no link specified, will use default link, like behavior in glm.
#' @param g_accuracy Level of accuracy to check prediction un-biasedness.
#' @param formula An optional formula to use for adjustment specified using as.formula("..."). This overrides strata_cols and covariate_cols.
#'
#' @import dplyr
#' @import magrittr
#' @export
robincar_glm <- function(df,
                         treat_col, response_col, strata_cols=NULL, covariate_cols=NULL,
                         car_scheme="simple", adj_method="heterogeneous", vcovHC="HC0",
                         covariate_to_include_strata=NULL,
                         g_family=stats::gaussian, g_accuracy=7, formula=NULL,
                         contrast_h=NULL, contrast_dh=NULL){

  .check.car_scheme(car_scheme)
  .check.adj_method.glm(adj_method)
  .check.vcovHC(vcovHC)

  data <- .make.data(
    df=df, classname="RoboDataGLM",
    treat_col=treat_col,
    response_col=response_col,
    strata_cols=strata_cols,
    covariate_cols=covariate_cols,
    formula=formula
  )
  validate(data)

  # Create model object
  model <- .make.model(
    data=data,
    adj_method=adj_method,
    car_scheme=car_scheme,
    vcovHC=vcovHC,
    covariate_to_include_strata=covariate_to_include_strata,
    g_family=g_family,
    g_accuracy=g_accuracy
  )

  # Perform adjustment
  result <- adjust(model, data)

  # Create transformation object
  if(!is.null(contrast_h)){
    c_result <- robincar_contrast(
      result=result,
      contrast_h=contrast_h,
      contrast_dh=contrast_dh
    )
    result <- list(main=result, contrast=c_result)
  }

  return(result)
}
