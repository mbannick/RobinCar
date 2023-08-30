#' Covariate-adjusted estimators for time to event data
#'
#' Estimate a covariate-adjusted hazard ratio (`adj_method="CL"`),
#' or a covariate-adjusted stratified hazard ratio (`adj_method="CSL"`).
#'
#'
#' @param df A data.frame with the required columns
#' @param treat_col Name of column in df with treatment variable
#' @param response_col Name of the column in df with response variable
#' @param event_col Name of column in df with event indicator
#'                  (0/FALSE=no event, 1/TRUE=event)
#' @param strata_cols Names of columns in df with strata variables
#' @param covariate_cols Names of columns in df with covariate variables
#' @param car_scheme Name of the type of covariate-adaptive randomization scheme. One of: "simple", "pocock-simon", "biased-coin", "permuted-block".
#' @param ref_arm Reference arm of the treatment group, defaults to NULL,
#'                which results in using the first element of `unique(data[, treat_col])`.
#' @param p_trt Treatment allocation ratio for the reference arm.
#' @param adj_method Adjustment method (one of "CL", "CSL")
#' @param interval Interval for uniroot function
#'
#' @export
robincar_covhr <- function(df,
                           treat_col, response_col, event_col,
                           strata_cols=NULL, covariate_cols=NULL,
                           p_trt=0.5, ref_arm=NULL,
                           car_scheme="simple",
                           adj_method="CL",
                           interval=c(-10, 10)){

  .check.car_scheme(car_scheme)

  data <- .make.data(
    df=df,
    classname="RoboDataTTE",
    treat_col=treat_col,
    response_col=response_col,
    event_col=event_col,
    strata_cols=strata_cols,
    covariate_cols=covariate_cols
  )
  validate(data, ref_arm)

  # Create model object
  model <- .make.model(
    data=data,
    adj_method=adj_method,
    car_scheme=car_scheme,
    p_trt=p_trt,
    ref_arm=ref_arm,
    interval=interval
  )

  # Append the CovHR classification
  # so that we perform estimation rather than testing
  class(model) <- c("CovHR", class(model))

  # Perform adjustment
  result <- adjust(model, data)
  result$original_df <- df

  return(result)

}
