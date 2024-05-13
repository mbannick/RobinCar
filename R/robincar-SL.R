#' BETA: Covariate adjustment using working models from the super learner libraries
#' through the AIPW package with cross-fitting.
#'
#' Estimate treatment-group-specific response means and (optionally)
#' treatment group contrasts using a generalized linear working model.
#'
#' *WARNING: This function is still under development and has not been extensively tested.*
#' This function currently only works for two treatment groups.
#' Before using this function, you must load the \link{SuperLearner} library with
#' `library(SuperLearner)`, otherwise the function call will fail.
#'
#' @param df A data.frame with the required columns
#' @param treat_col Name of column in df with treatment variable
#' @param response_col Name of the column in df with response variable
#' @param car_strata_cols Names of columns in df with car_strata variables
#' @param covariate_cols Names of columns in df with covariate variables
#' @param car_scheme Name of the type of covariate-adaptive randomization scheme. One of: "simple", "pocock-simon", "biased-coin", "permuted-block".
#' @param covariate_to_include_strata Whether to include car_strata variables in covariate adjustment. Defaults to F for ANOVA and ANCOVA; defaults to T for ANHECOVA. User may override by passing in this argument.
#' @param contrast_h An optional function to specify a desired contrast
#' @param contrast_dh An optional jacobian function for the contrast (otherwise use numerical derivative)
#' @param SL_libraries Vector of super-learner libraries to use for the covariate adjustment (see \link[SuperLearner:listWrappers]{SuperLearner::listWrappers})
#' @param SL_learners Optional list of super-learner "learners" to use for the covariate adjustment (see \link[SuperLearner:create.Learner]{SuperLearner::create.Learner())}
#' @param k_split Number of splits to use in cross-fitting
#' @param g_accuracy Level of accuracy to check prediction un-biasedness (in digits).
#'
#'
#' @export
#' @importFrom SuperLearner listWrappers
#' @import SuperLearner
#' @import AIPW
#'
#' @examples
#'
#' library(SuperLearner)
#' library(ranger)
#' n <- 1000
#' set.seed(10)
#' DATA2 <- data.frame(A=rbinom(n, size=1, prob=0.5),
#'                     y=rbinom(n, size=1, prob=0.2),
#'                     x1=rnorm(n),
#'                     x2=rnorm(n),
#'                     x3=as.factor(rbinom(n, size=1, prob=0.5)),
#'                     z1=rbinom(n, size=1, prob=0.5),
#'                     z2=rbinom(n, size=1, prob=0.5))
#' DATA2[, "y"] <- NA
#' As <- DATA2$A == 1
#' DATA2[DATA2$A == 1, "y"] <- rbinom(
#'   sum(As),
#'   size=1,
#'   prob=exp(DATA2[As,]$x1)/(1+exp(DATA2[As,]$x1)))
#' DATA2[DATA2$A == 0, "y"] <- rbinom(
#'   n-sum(As),
#'   size=1,
#'   prob=exp(1 +
#'     5*DATA2[!As,]$x1 + DATA2[!As,]$x2)/
#'     (1+exp(1 + 5*DATA2[!As,]$x1 + DATA2[!As,]$x2)))
#' DATA2$A <- as.factor(DATA2$A)
#'
#' sl.mod <- robincar_SL(
#'   df=DATA2,
#'   response_col="y",
#'   treat_col="A",
#'   car_strata_cols=c("z1"),
#'   covariate_cols=c("x1"),
#'   SL_libraries=c("SL.ranger"),
#'   car_scheme="permuted-block",
#'   covariate_to_include_strata=TRUE
#' )
#'
#' sl.mod$result
#'
#' @returns See value of \link[RobinCar:robincar_glm]{RobinCar::robincar_glm}, but the working model for \eqn{\hat{\mu}(X_i)} is based on the \link{AIPW} package that uses specified SuperLearner libraries and cross-fitting.
#' Also, `mod` attribute is an object of class \link[AIPW:AIPW]{AIPW::AIPW}.
robincar_SL <- function(df,
                        treat_col, response_col, car_strata_cols=NULL, covariate_cols=NULL,
                        car_scheme="simple",
                        covariate_to_include_strata=NULL,
                        SL_libraries=c(), SL_learners=c(),
                        k_split=2,
                        g_accuracy=7,
                        contrast_h=NULL, contrast_dh=NULL){

  .check.car_scheme(car_scheme)
  # .check.sl.libraries(SL_libraries)

  if(length(SL_libraries) + length(SL_learners) == 0){
    stop("Must provide at least one SL_libraries or SL_learners.")
  }

  # Add a k split index to .make.data, so that we can index
  # the dataset on the kth index set for cross-fitting
  data <- .make.data(
    df=df, classname="RoboDataSL",
    treat_col=treat_col,
    response_col=response_col,
    car_strata_cols=car_strata_cols,
    covariate_cols=covariate_cols
  )
  validate(data)

  # Create model object
  model <- .make.model(
    data=data,
    car_scheme=car_scheme,
    vcovHC="HC0",
    covariate_to_include_strata=covariate_to_include_strata,
    SL_libraries=SL_libraries,
    SL_learners=SL_learners,
    k_split=k_split,
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

#' BETA: Covariate adjustment using working models from the super learner libraries
#' through the AIPW package with cross-fitting, with median adjustment.
#'
#' Estimate treatment-group-specific response means and (optionally)
#' treatment group contrasts using a generalized linear working model.
#' Perform median adjustment to limit randomness induced from cross-fitting.
#'
#' *WARNING: This function is still under development and has not been extensively tested.*
#' This function currently only works for two treatment groups.
#' Before using this function, you must load the SuperLearner library with
#' `library(SuperLearner)`, otherwise the function call will fail.
#'
#' @param n_times Number of times to run the robincar_SL function
#' @param seed Seed to set before running the set of functions
#' @inheritParams robincar_SL
#' @export
#'
#' @returns See value of \link[RobinCar:robincar_SL]{RobinCar::robincar_SL}.
#' Attributes `mods` and `mu_as` are lists of `mod` and `mu_a` attributes, respectively, for each replicate of `robincar_SL` used in the median.
robincar_SL_median <- function(n_times, seed, df,
                               treat_col, response_col,
                               car_strata_cols=NULL, covariate_cols=NULL,
                               car_scheme="simple",
                               covariate_to_include_strata=NULL,
                               SL_libraries=c(), SL_learners=c(),
                               k_split=2,
                               g_accuracy=7,
                               contrast_h=NULL, contrast_dh=NULL){

  if(n_times %% 2 == 0) stop("Must have an odd number of n_times.")

  # Set the seed before running any of the robincar function
  # that will do sample splitting.
  set.seed(seed)
  res <- list()
  for(i in 1:n_times){
    cat(".")
    res[[i]] <- robincar_SL(
      df=df,
      treat_col=treat_col,
      response_col=response_col,
      car_strata_cols=car_strata_cols,
      covariate_cols=covariate_cols,
      car_scheme=car_scheme,
      covariate_to_include_strata=covariate_to_include_strata,
      SL_libraries=SL_libraries,
      SL_learners=SL_learners,
      k_split=k_split,
      g_accuracy=g_accuracy,
      contrast_h=contrast_h,
      contrast_dh=contrast_dh)
  }

  estimates <- sapply(res, function(x) x$result$estimate)
  varcovs <- lapply(res, function(x) x$varcov)

  # Get median of estimates
  estimate <- apply(estimates, 1, stats::median)

  # Create varcov_tilde which adds on the residual
  varcov_tilde <- list()
  for(i in 1:n_times){
    c_est <- estimates[, i] - estimate # centered estimates
    varcov_tilde[[i]] <- varcovs[[i]] + c_est %*% t(c_est)
  }
  varcov_opnorm <- sapply(varcov_tilde, function(x) norm(x, type="2"))
  # TODO: What if there are multiple matrices with the same operator norm?
  idx <- which(varcov_opnorm == stats::median(varcov_opnorm))
  variance <- varcov_tilde[[idx]]

  data <- res[[1]]$data

  result <- format_results(data$treat_levels, estimate, variance)

  mods <- lapply(res, function(x) x$mod)
  mu_as <- lapply(res, function(x) x$mu_a)

  return(
    structure(
      class="SLModelResult",
      list(result=result, varcov=variance,
           settings=res[[1]]$settings,
           data=res[[1]]$data, mod=mods, mu_a=mu_as,
           original_df=res[[1]]$original_df)
    )
  )
}
