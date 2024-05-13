#' Covariate adjustment for time to event data
#'
#' Perform a covariate-adjusted logrank test (`adj_method="CL"`),
#' covariate-adjusted stratified logrank test (`adj_method="CSL"`),
#' or a covariate-adjusted robust Cox score test (`adj_method="coxscore"`).
#'
#' `robincar_coxscore` and `robincar_logrank` are wrapper functions around
#' `robincar_tte`.
#'
#' @param df A data.frame with the required columns
#' @param treat_col Name of column in df with treatment variable
#' @param response_col Name of the column in df with response variable
#' @param event_col Name of column in df with event indicator
#'                  (0/FALSE=no event, 1/TRUE=event)
#' @param car_strata_cols Names of columns in df with car_strata variables
#' @param covariate_cols Names of columns in df with covariate variables
#' @param car_scheme Name of the type of covariate-adaptive randomization scheme. One of: "simple", "pocock-simon", "biased-coin", "permuted-block".
#' @param ref_arm Reference arm of the treatment group, defaults to NULL,
#'                which results in using the first element of `unique(data[, treat_col])`.
#' @param p_trt Treatment allocation ratio for the reference arm.
#' @param adj_method Adjustment method (one of "CL", "CSL", or "coxscore")
#' @param sparse_remove Remove sparse car_strata from calculation
#'
#' @export
#'
#' @returns For adjustment method "CL" or "CSL", see value of \link[RobinCar:robincar_linear]{RobinCar::robincar_logrank()}; for adjustment method "coxscore" see value of \link[RobinCar:robincar_coxscore]{RobinCar::robincar_coxscore()}.
robincar_tte <- function(df,
                         treat_col, response_col, event_col,
                         adj_method,
                         car_strata_cols=NULL, covariate_cols=NULL,
                         p_trt=0.5, ref_arm=NULL, sparse_remove=TRUE,
                         car_scheme="simple"){

  .check.car_scheme(car_scheme)

  data <- .make.data(
    df=df,
    classname="RoboDataTTE",
    treat_col=treat_col,
    response_col=response_col,
    event_col=event_col,
    car_strata_cols=car_strata_cols,
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
    sparse_remove=sparse_remove
  )

  # Perform adjustment
  result <- adjust(model, data)
  result$original_df <- df

  return(result)

}
