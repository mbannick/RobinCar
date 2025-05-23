#' Covariate adjustment using generalized linear working model
#'
#' Estimate treatment-group-specific response means and (optionally)
#' treatment group contrasts using a generalized linear working model.
#'
#' The output is the AIPW estimator given by (for each treatment group \eqn{a}):
#'
#' \deqn{\frac{1}{n} \sum_{i=1}^{n} \hat{\mu}_a(X_i) + \frac{1}{n_a} \sum_{i:A_i=a} \{Y_i - \hat{\mu}(X_i)\}}
#'
#' where \eqn{Y_i} is the outcome, \eqn{A_i} is the treatment assignment, \eqn{X_i} are the covariates, \eqn{n_a = \sum_{i=1}^{n} A_i=a},
#' and \eqn{\hat{\mu}_a} is the estimated conditional mean function based on the GLM working model.
#' This working model has treatment \eqn{a}-specific coefficients if `adj_method` is "heterogeneous". Otherwise, they are shared across the treatment arms.
#' Alternatively, if `formula` is used, the working model can be specified according to the user.
#'
#' Importantly, the estimated variance accounts for misspecification of the working model, and for covariate-adaptive randomization.
#' The variance estimator is given by
#' \deqn{\hat{V} = \hat{V}_{\rm SR} - \hat{V}_{\Omega}}
#' where \eqn{\hat{V}_{\rm SR}} is the contribution to the variance under simple randomization, and \eqn{\hat{V}_{\Omega}} is a term
#' that only appears when a covariate-adaptive randomization scheme is used.
#' The \eqn{\hat{V}_{\Omega}} is the second line of \eqn{\hat{V}} in \insertCite{bannickGeneralFormCovariate2025}{RobinCar}.
#'
#' There are three different estimators available for \eqn{\hat{V}_{\rm SR}}, which the user
#' can choose with the argument \code{variance_type}. We describe these here.
#'
#' The three variance types are given as follows:
#'
#' \itemize{
#'  \item{
#'    Type 1 (default):
#'    \deqn{\mathrm{diag}\left[\hat{\pi}_a^{-1} (\mathrm{Var}_a(Y_i) - 2\hat{Q}_{a,a} + \hat{\Sigma}_{a,a} ), a = 1, \dots, K \right] + \hat{Q} + \hat{Q}^T - \hat{\Sigma}}
#'  }
#'  \item{
#'    Type 2:
#'    \deqn{\mathrm{diag}\left[\hat{\pi}_a^{-1} (\mathrm{Var}_a(Y_i - \hat{\mu}_a(X_i)) - 2\hat{Q}_{a,a} + \hat{\Sigma}_{a,a} ), a = 1, \dots, K \right] + \hat{Q} + \hat{Q}^T - \hat{\Sigma}}
#'  }
#'  \item{
#'    Type 3:
#'    \deqn{\mathrm{diag}\left[\hat{\pi}_a^{-1} E_a([Y_i - \hat{\mu}_a(X_i)]^2), a = 1, \dots, K \right] + \hat{A}}
#'  }
#' }
#'
#' where \eqn{\hat{\pi}_a} is the treatment proportion for group a,
#' \eqn{\hat{Q}_{a,b} = \mathrm{Cov}_a(Y_i, \hat{\mu}_b(X_i))},
#' \eqn{\hat{\Sigma}_{a,b} = \mathrm{Cov}(\hat{\mu}_a(X_i), \hat{\mu}_b(X_i))},
#' and the matrix \eqn{\hat{A}} has diagonal entries for \eqn{(a, a)} given by
#' \deqn{2\mathrm{Cov}_a(Y_i - \hat{\mu}_a(X_i), \hat{\mu}_a(X_i)) + \mathrm{Var}_a(\hat{\mu}_a(X_i))}
#' and off-diagonal entries for \eqn{(a, b)} given by
#' \deqn{\mathrm{Cov}_a(Y_i, \hat{\mu}_b(X_i)) + \mathrm{Cov}_b(Y_i, \hat{\mu}_a(X_i)) - (1/2) \left[\mathrm{Cov}_a(\hat{\mu}_a(X_i), \hat{\mu}_b(X_i)) + \mathrm{Cov}_b(\hat{\mu}_a(X_i), \hat{\mu}_b(X_i)) \right].}
#' We use \eqn{E_a}, \eqn{\mathrm{Var}_a}, and \eqn{\mathrm{Cov}_a} to refer to the empirical expectation, variance, and
#'  covariance among observations in group a only, and \eqn{\mathrm{Cov}} is the covariance within
#' the entire sample.
#'
#' Please see the Supplemental Material Sect. H of \insertCite{bannickGeneralFormCovariate2025}{RobinCar} for a discussion
#' of the merits of each type of variance estimator. Briefly, we recommend
#' variance types 1 generally, and variance type 3 if it is anticipated
#' that the distribution of \eqn{X} varies substantially over treatment groups.
#'
#' @param df A data.frame with the required columns
#' @param treat_col Name of column in df with treatment variable
#' @param response_col Name of the column in df with response variable
#' @param formula The formula to use for adjustment specified using as.formula("..."). This overrides car_strata_cols and covariate_cols.
#' @param car_strata_cols Names of columns in df with car_strata variables
#' @param car_scheme Name of the type of covariate-adaptive randomization scheme. One of: "simple", "pocock-simon", "biased-coin", "permuted-block".
#' @param contrast_h An optional function to specify a desired contrast
#' @param contrast_dh An optional jacobian function for the contrast (otherwise use numerical derivative)
#' @param g_family Family that would be supplied to glm(...), e.g., binomial. If no link specified, will use default link, like behavior in glm.
#'                 If you wish to use a negative binomial working model with an unknown dispersion parameter, then use `g_family="nb"`.
#' @param g_accuracy Level of accuracy to check prediction un-biasedness.
#' @param variance_type The type of variance estimator to use, type 1, 2, or 3. All three are asymptotically equivalent. See details for more.
#' @export
#'
#' @importFrom Rdpack reprompt
#'
#' @returns If `contrast_h` argument is used, outputs a `main` and a `contrast` object. The `main` object has the following structure:
#'
#'  \item{result}{A \link[dplyr:tibble]{dplyr::tibble()} with the treatment label, treatment mean estimate using AIPW, estimated SE, and p-value based on a z-test with estimate and SE.}
#'  \item{varcov}{The variance-covariance matrix for the treatment mean estimates.}
#'  \item{settings}{List of model settings used in covariate adjustment.}
#'  \item{original_df}{The original dataset provided by the user.}
#'  \item{mod}{The fit from the \link[stats:glm]{glm()} working model used for covariate adjustment.}
#'  \item{mu_a}{Predicted potential outcomes for each treatment category (columns) and individual (rows). These are the \eqn{\hat{\mu}_a}}.
#'  \item{g.estimate}{The G-computation estimate based only on \eqn{\frac{1}{n} \sum_{i=1}^{n} \hat{\mu}_a(X_i)}. This is equivalent to the AIPW estimate when a canonical link function is used.}
#'  \item{data}{Attributes about the dataset.}
#'
#'  The `contrast` object has a structure that is documented in \link[RobinCar:robincar_contrast]{RobinCar::robincar_contrast()}.
#'
#' @references
#'   \insertAllCited{}
robincar_glm <- function(df,
                         treat_col, response_col,
                         formula=NULL, car_strata_cols=NULL,
                         car_scheme="simple",
                         g_family=stats::gaussian, g_accuracy=7,
                         contrast_h=NULL, contrast_dh=NULL,
                         variance_type=1){

  .check.car_scheme(car_scheme, car_strata_cols)
  .check.variance_type(variance_type)

  vcovHC <- "HC0"

  data <- .make.data(
    df=df, classname="RoboDataGLM",
    treat_col=treat_col,
    response_col=response_col,
    car_strata_cols=car_strata_cols,
    formula=formula
  )
  validate(data)

  # Create model object
  model <- .make.model(
    data=data,
    car_scheme=car_scheme,
    vcovHC=vcovHC,
    g_family=g_family,
    g_accuracy=g_accuracy,
    variance_type=variance_type
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
