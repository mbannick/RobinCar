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
#' @param car_strata_cols Names of columns in df with car_strata variables
#' @param covariate_cols Names of columns in df with covariate variables
#' @param car_scheme Name of the type of covariate-adaptive randomization scheme. One of: "simple", "pocock-simon", "biased-coin", "permuted-block".
#' @param adj_method Name of adjustment method to use, one of "heterogeneous" (interaction model) or "homogeneous"
# @param vcovHC Type of heteroskedasticity-consistent variance estimates. One of: "HC0", "HC1", "HC3".
#' @param covariate_to_include_strata Whether to include car_strata variables in covariate adjustment. Defaults to F for homogeneous; defaults to T for heterogeneous. User may override by passing in this argument.
#' @param contrast_h An optional function to specify a desired contrast
#' @param contrast_dh An optional jacobian function for the contrast (otherwise use numerical derivative)
#' @param g_family Family that would be supplied to glm(...), e.g., binomial. If no link specified, will use default link, like behavior in glm.
#'                 If you wish to use a negative binomial working model with an unknown dispersion parameter, then use `g_family="nb"`.
#' @param g_accuracy Level of accuracy to check prediction un-biasedness.
#' @param formula An optional formula to use for adjustment specified using as.formula("..."). This overrides car_strata_cols and covariate_cols.
#'
#' @export
robincar_glm <- function(df,
                         treat_col, response_col, car_strata_cols=NULL, covariate_cols=NULL,
                         car_scheme="simple", adj_method="heterogeneous", # vcovHC="HC0",
                         covariate_to_include_strata=NULL,
                         g_family=stats::gaussian, g_accuracy=7, formula=NULL,
                         contrast_h=NULL, contrast_dh=NULL){

  .check.car_scheme(car_scheme)
  .check.adj_method.glm(adj_method)
  # .check.vcovHC(vcovHC)
  # for now, until we implement a better degrees of freedom correction
  # only allowing HC0.
  vcovHC <- "HC0"

  data <- .make.data(
    df=df, classname="RoboDataGLM",
    treat_col=treat_col,
    response_col=response_col,
    car_strata_cols=car_strata_cols,
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
  # This is to save the original dataset for calibration add_x later on
  result$original_df <- df

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
