#' Covariate adjustment using linear working model
#'
#' Estimate treatment-group-specific response means and (optionally)
#' treatment group contrasts using a linear working model for continuous outcomes.
#'
#' * Adjustment method "ANOVA" fits a linear model with formula `Y ~ A` where
#' `A` is the treatment group indicator and `Y` is the response.
#' * "ANCOVA" fits a linear model with `Y ~ A + X` where `X` are the variables
#' specified in the `covariate_cols` argument.
#' * "ANHECOVA" fits a linear model with `Y ~ A*X`, the main effects and treatment-by-covariate
#' interactions.
#'
#' @param df A data.frame with the required columns
#' @param treat_col Name of column in df with treatment variable
#' @param response_col Name of the column in df with response variable
#' @param car_strata_cols Names of columns in df with car_strata variables
#' @param covariate_cols Names of columns in df with covariate variables. **If you want to include the strata variables as covariates also, add them here.**
#' @param car_scheme Name of the type of covariate-adaptive randomization scheme. One of: "simple", "pocock-simon", "biased-coin", "permuted-block".
#' @param adj_method Name of linear adjustment method to use. One of: "ANOVA", "ANCOVA", "ANHECOVA".
#' @param contrast_h An optional function to specify a desired contrast
#' @param contrast_dh An optional jacobian function for the contrast (otherwise use numerical derivative)
#' @export
#'
#' @returns See value of \link[RobinCar:robincar_glm]{RobinCar::robincar_glm()}, this function is a wrapper
#' using a linear link function.
robincar_linear <- function(df,
                             treat_col, response_col, car_strata_cols=NULL, covariate_cols=NULL,
                             car_scheme="simple", adj_method="ANOVA",
                             contrast_h=NULL, contrast_dh=NULL){

  # Define the formula based on adj_method:
  formula <- .create.formula(
    adj_method=adj_method,
    treat_col=treat_col,
    response_col=response_col,
    covariate_cols=covariate_cols)

  obj <- robincar_glm(
    df=df,
    treat_col=treat_col,
    response_col=response_col,
    car_strata_cols=car_strata_cols,
    car_scheme=car_scheme,
    formula=formula,
    contrast_h=contrast_h,
    contrast_dh=contrast_dh,
    g_family=stats::gaussian,
    g_accuracy=7
  )

  return(obj)
}
