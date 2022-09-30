#' Covariate adjustment using generalized linear working model
#'
#' WARNING: this function is still under development.
#' Estimate treatment-group-specific response means and (optionally)
#' treatment group contrasts using a generalized linear working model.
#'
#' @param inheritParams robincar_coxscore
#' @param adj_method Adjustment method (one of "CL", "CSL", or "coxscore")
#'
#' @import dplyr
#' @import magrittr
#' @export
robincar_tte <- function(df,
                         treat_col, response_col, event_col,
                         adj_method,
                         strata_cols=NULL, covariate_cols=NULL,
                         p_trt=0.5, ref_arm=NULL, sparse_remove=TRUE,
                         car_scheme="simple"){

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
    sparse_remove=sparse_remove
  )

  # Perform adjustment
  result <- adjust(model, data)

  return(result)

}
