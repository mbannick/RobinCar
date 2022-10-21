#' Robust cox score adjustment
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
#'
#' @export
robincar_coxscore <- function(...){

  result <- robincar_tte(
    adj_method="coxscore", ...
  )

  return(result)

}
