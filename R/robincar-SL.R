#' Covariate adjustment using working models from the super learner libraries
#' through the AIPW package with cross-fitting.
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
#' @param vcovHC Type of heteroskedasticity-consistent variance estimates. One of: "HC0", "HC1", "HC3".
#' @param covariate_to_include_strata Whether to include strata variables in covariate adjustment. Defaults to F for ANOVA and ANCOVA; defaults to T for ANHECOVA. User may override by passing in this argument.
#' @param contrast_h An optional function to specify a desired contrast
#' @param contrast_dh An optional jacobian function for the contrast (otherwise use numerical derivative)
#' @param SL_libraries Vector of super-learner libraries to use for the covariate adjustment (see SuperLearner::listWrappers())
#' @param k_split Number of splits to use in cross-fitting
#' @param centering WARNING: Advanced use only. By default is NULL, so will perform expected behavior
#'                           in documentation.
#'                  If centering = 'none', then this will force it to be a g-computation estimator,
#'                     even if prediction unbiasedness does not hold.
#'                  If centering = 'tx', then this will center within treatment groups, which is the
#'                     AIPW estimator (equivalent to g-computation if prediction unbiasedness holds).
#'                  TODO: If centering = 'tx-strata', then this will center within treatment group and
#'                     within joint strata, which allows for the AIPW estimator to be used with Pocock-Simon.
#' @export
#' @importFrom SuperLearner listWrappers
#' @import SuperLearner
robincar_SL <- function(df,
                         treat_col, response_col, strata_cols=NULL, covariate_cols=NULL,
                         car_scheme="simple", vcovHC="HC0",
                         covariate_to_include_strata=NULL,
                         SL_libraries=c(), k_split=2,
                         contrast_h=NULL, contrast_dh=NULL){

  .check.car_scheme(car_scheme)
  .check.vcovHC(vcovHC)
  .check.sl.libraries(SL_libraries)

  # Add a k split index to .make.data, so that we can index
  # the dataset on the kth index set for cross-fitting
  data <- .make.data(
    df=df, classname="RoboDataSL",
    treat_col=treat_col,
    response_col=response_col,
    strata_cols=strata_cols,
    covariate_cols=covariate_cols
  )
  validate(data)

  # Create model object
  model <- .make.model(
    data=data,
    car_scheme=car_scheme,
    vcovHC=vcovHC,
    covariate_to_include_strata=covariate_to_include_strata,
    SL_libraries=SL_libraries,
    k_split=k_split
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

robincar_SL_median <- function(n_times, seed, ...){

  # Set the seed before running any of the robincar function
  # that will do sample splitting.
  set.seed(seed)
  res <- list()
  for(i in n_times){
    res[[i]] <- robincar_SL(...)
  }

  # TODO: Extract each of the point estimates and variance-covariance
  # matrices to get the median-adjusted results.

}
