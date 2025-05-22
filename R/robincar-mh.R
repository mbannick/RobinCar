#' Estimate Mantel-Haenszel Risk Difference
#'
#' @description
#' This function estimates Mantel-Haenszel risk difference and average treatment effect.
#'
#' @param df A data.frame with the required columns
#' @param treat_col Name of column in df with treatment variable. Must be binary
#' @param response_col Name of the column in df with response variable
#' @param strata_cols Names of columns in df with strata variables
#' @param estimand A character string specifying the estimand. One of "MH" or "ATE" (default). See Details
#' @param ci_type A character string specifying the type of confidence interval. One of "GR", "mGR" (default), "Sato"

#' @details
#' The estimand of interest can be either Mantel-Haenszel risk difference or Average Treatment Effect (ATE). 
#' The latter is the default option of `estimand`. When `estimand="ATE"`, `ci_type` is limited to the modified Greenland variance estimator (mGR).
#'  Otherwise, Greenland's variance estimator (GR) and Sato's variance estimator are optional.
#' 
#' @examples
#' df <- RobinCar:::data_sim
#' df$y_bin = ifelse(df$y>2.5, 1, 0)
#' robincar_mh(df = df[df$A!=2,],
#'             treat_col = "A",
#'             response_col = "y_bin",
#'             strata_cols = c("z1", "z2"),
#'             estimand = "MH",
#'             ci_type = "mGR")
#' 
#' @export
robincar_mh <- function(df, treat_col, response_col, strata_cols, estimand="ATE", ci_type='mGR') {
  
  .check.estimand_ci(estimand, ci_type)
  
  data <- .make.data(
    df = df,
    classname = "RoboDataMH",
    treat_col = treat_col,
    response_col = response_col,
    car_strata_cols = strata_cols,
    formula = NULL
  )
  validate(data)
  
  model <- .make.model(
    data = data,
    estimand = estimand,
    ci_type = ci_type,
    strata_cols = strata_cols
  )
  
  
  result <- adjust(model, data)
  
  return(result)
}
