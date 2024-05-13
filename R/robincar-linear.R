#' Covariate adjustment using linear working model
#'
#' Estimate treatment-group-specific response means and (optionally)
#' treatment group contrasts using a linear working model for continuous outcomes.
#'
#' @param df A data.frame with the required columns
#' @param treat_col Name of column in df with treatment variable
#' @param response_col Name of the column in df with response variable
#' @param car_strata_cols Names of columns in df with car_strata variables
#' @param covariate_cols Names of columns in df with covariate variables
#' @param car_scheme Name of the type of covariate-adaptive randomization scheme. One of: "simple", "pocock-simon", "biased-coin", "permuted-block".
#' @param adj_method Name of linear adjustment method to use. One of: "ANOVA", "ANCOVA", "ANHECOVA".
#' @param covariate_to_include_strata Whether to include car_strata variables in covariate adjustment. Defaults to F for ANOVA and ANCOVA; defaults to T for ANHECOVA. User may override by passing in this argument.
#' @param contrast_h An optional function to specify a desired contrast
#' @param contrast_dh An optional jacobian function for the contrast (otherwise use numerical derivative)
#' @export
#'
#' @returns See value of \link[RobinCar:robincar_glm]{RobinCar::robincar_glm()}, this function is a wrapper using an identity link for the glm.
robincar_linear <- function(df,
                            treat_col, response_col, car_strata_cols=NULL, covariate_cols=NULL,
                            car_scheme="simple", adj_method="ANOVA",
                            covariate_to_include_strata=NULL,
                            contrast_h=NULL, contrast_dh=NULL){

  .check.adj_method.linear(adj_method)
  if(adj_method == "ANOVA"){
    covariate_cols <- NULL
    covariate_to_include_strata <- FALSE
    adj_method <- "homogeneous"
  } else if(adj_method == "ANCOVA"){
    adj_method <- "homogeneous"
  } else {
    adj_method <- "heterogeneous"
  }

  result <- robincar_glm(
    df=df,
    treat_col=treat_col,
    response_col=response_col,
    covariate_cols=covariate_cols,
    car_strata_cols=car_strata_cols,
    car_scheme=car_scheme,
    covariate_to_include_strata=covariate_to_include_strata,
    contrast_h=contrast_h,
    contrast_dh=contrast_dh,
    g_family=stats::gaussian,
    g_accuracy=7,
    adj_method=adj_method
  )

  return(result)
}

#' Covariate adjustment using linear working model, with simplified interface.
#'
#' Estimate treatment-group-specific response means and (optionally)
#' treatment group contrasts using a linear working model for continuous outcomes.
#'
#' Note the difference between \link[RobinCar:robincar_linear]{RobinCar::robincar_linear()} and \link[RobinCar:robincar_linear2]{RobinCar::robincar_linear2()} with regards to their `covariate_cols` argument.
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
#' @returns See value of \link[RobinCar:robincar_linear]{RobinCar::robincar_linear()}, this function is a wrapper.
robincar_linear2 <- function(df,
                             treat_col, response_col, car_strata_cols=NULL, covariate_cols=NULL,
                             car_scheme="simple", adj_method="ANOVA",
                             contrast_h=NULL, contrast_dh=NULL){

  obj <- robincar_linear(
    df=df,
    treat_col=treat_col,
    response_col=response_col,
    car_strata_cols=car_strata_cols,
    covariate_cols=covariate_cols,
    car_scheme=car_scheme,
    adj_method=adj_method,
    covariate_to_include_strata=FALSE,
    contrast_h=contrast_h,
    contrast_dh=contrast_dh
  )

  return(obj)
}
