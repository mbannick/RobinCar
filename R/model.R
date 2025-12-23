.make.model <- function(data, ...){
  UseMethod(".make.model", data)
}

# Makes a model class for the specified glm adjustment method
# with settings for covariate randomization
# scheme and vcovHC type.
.make.model.RoboDataGLM <- function(data, car_scheme, vcovHC,
                                    g_family, g_accuracy, variance_type, ...) {

  # Get logic for car scheme + adjustment methods
  z_exists <- !is.null(data$car_strata)
  logic <- glmlogic(car_scheme=car_scheme, z_exists=z_exists)

  model <- structure(
    list(
      vcovHC=vcovHC,
      g_family=g_family,
      g_accuracy=g_accuracy,
      car_scheme=car_scheme,
      pu_joint_z=logic$pu_joint_z,
      pu_funcs=logic$pu_funcs,
      omegaz_func=logic$omegaz_func,
      variance_type=variance_type
    ),
    class=c(logic$method)
  )

  return(model)
}

# Makes a model class for the specified glm adjustment method
# with settings for covariate randomization
# scheme and vcovHC type.
.make.model.RoboDataSL <- function(data, car_scheme, vcovHC,
                                   covariate_to_include_strata,
                                   SL_libraries, SL_learners,
                                   k_split,
                                   g_accuracy, variance_type, ...) {

  x_exists <- !is.null(data$covariate)
  z_exists <- !is.null(data$car_strata)

  # Get logic for adjustment methods
  logic <- SLlogic(car_scheme=car_scheme,
                    x_exists=x_exists, z_exists=z_exists,
                    cov_strata=covariate_to_include_strata)
  model <- structure(
    list(
      vcovHC=vcovHC,
      SL_libraries=SL_libraries,
      SL_learners=SL_learners,
      k_split=k_split,
      car_scheme=car_scheme,
      adj_se_z=logic$adj_se_z,
      adj_vars=logic$adj_vars,
      pu_joint_z=logic$pu_joint_z,
      pu_funcs=logic$pu_funcs,
      omegaz_func=logic$omegaz_func,
      g_accuracy=g_accuracy,
      variance_type=variance_type
    ),
    class=c(logic$method)
  )

  return(model)
}


# Makes a model class for the specified TTE adjustment method
# with settings for covariate randomization
# scheme and vcovHC type.
.make.model.RoboDataTTE <- function(data, adj_method, car_scheme,
                                    p_trt, ref_arm, return_influence, ...) {

  x_exists <- !is.null(data$covariate)
  z_exists <- !is.null(data$car_strata)

  # Get logic for adjustment methods
  logic <- ttelogic(adj_method=adj_method, car_scheme=car_scheme,
                    x_exists=x_exists, z_exists=z_exists)

  if(adj_method %in% c("CL", "CSL")){
    classtype <- "LogRank"
  } else if(adj_method %in% c("coxscore")) {
    classtype <- "CoxScore"
  } else {
    stop("Unrecognized adjustment method.")
  }

  model <- structure(
    list(
      method=logic$method,
      adj_cov=logic$adj_cov,
      adj_strata=logic$adj_strata,
      car_strata=logic$car_strata,
      car_scheme=car_scheme,
      p_trt=p_trt,
      ref_arm=ref_arm,
      return_influence=return_influence,
      ...
    ),
    class=c(classtype, logic$method)
  )

  if ("return_influence" %in% names(model)){
    if (!is.logical(model$return_influence)) stop("return_influence must be either TRUE or FALSE.")
    if (model$return_influence) {
      if (!"id_col" %in% names(data)) stop("id_col must not be NULL if return_influence is TRUE")
      if (!is.null(data$covariate)){
        if ("id" %in% colnames(data$covariate)){
          stop("if return_influence is TRUE, no covariate can be called 'id'.")
        }
      }
    }
  }
  
  return(model)
}


# Makes a model class for the MH method
# with settings for estimand and CI type
.make.model.RoboDataMH <- function(data, estimand, ci_type, strata_cols, ...) {
  
  model <- structure(
    list(
      estimand = estimand,
      ci_type = ci_type,
      strata_cols = strata_cols
    ),
    class=c("MHModel")
  )
  
  return(model)
}
