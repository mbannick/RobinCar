# These are the covariate randomization schemes
# that are applicable; leaves out simple randomization
# and pocock-simon minimization because those are dealt with separately
APPLICABLE.SCHEMES <- c("permuted-block",
                        "biased-coin",
                        "urn")

glmlogic <- function(adj_method, x_exists, z_exists, car_scheme, cov_strata,
                     formula=NULL){

  # Method to use for the GLM G-Computation adjustment,
  # and whether we should do AIPW on top of it.
  if(adj_method == "heterogeneous"){
    method <- "ANHECOVA"
  } else {
    method <- "ANCOVA"
  }
  adj_vars <- NULL

  aipw <- TRUE
  adj_se_z <- FALSE

  # Whether to check prediction un-biasedness in all joint levels
  # of Z, or just overall, and whether to report an error or warning
  pu_joint_z <- FALSE
  # Biased action will be one of "", "warning", "error"
  pu_funcs <- function(){}

  omegaz_func <- NULL

  if(!is.null(formula)){
    method <- "CUSTOM"
    adj_vars <- "formula"
    if(x_exists) .form.warn()
    if(z_exists & car_scheme == "pocock-simon"){
      pu_funcs <- c(.pu.z.calibrate, .pu.z.warn)
      pu_joint_z <- TRUE
    } else {
      pu_funcs <- .pu.warn
      pu_joint_z <- FALSE
    }
    if(z_exists & (car_scheme != "simple")){
      omegaz_func <- omegaz.closure(car_scheme)
    }
  } else {
    if(car_scheme == "simple"){
      if(z_exists) .z.exist.warn.simple()
      if(x_exists){
        adj_vars <- "x"
        pu_funcs <- .pu.warn
      } else {
        # .x.miss.warn()
        method <- "ANOVA"
        aipw <- FALSE
      }
    }
    if(car_scheme == "pocock-simon"){
      if(!z_exists){
        .z.miss.err()
      } else {
        pu_joint_z <- TRUE
        omegaz_func <- omegaz.closure(car_scheme)
        if(cov_strata){
          if(x_exists){
            adj_vars <- "joint_z_x"
            pu_funcs <- .pu.z.warn
          } else {
            adj_vars <- "joint_z"
            pu_funcs <- c(.pu.z.calibrate, .pu.z.warn)
          }
        } else {
          if(x_exists){
            adj_vars <- "x"
            pu_funcs <- .pu.z.warn
          } else {
            method <- "ANOVA"
            pu_funcs <- .pu.z.warn
          }
        }
      }
    }
    if(car_scheme %in% APPLICABLE.SCHEMES){

      if(z_exists){
        omegaz_func <- omegaz.closure(car_scheme)

        adj_se_z <- TRUE
        if(x_exists){
          if(cov_strata){
            adj_vars <- "joint_z_x"
            pu_funcs <- .pu.warn
          } else {
            adj_vars <- "x"
          }
        } else {
          if(cov_strata){
            adj_vars <- "joint_z"
            pu_funcs <- .pu.warn
          } else {
            method <- "ANOVA"
            aipw <- FALSE
          }
        }
      } else {
        .z.miss.err()
      }
    }
  }
  if(aipw){
    method=c("AIPW", "GCOMP", "GLMModel", method)
  } else {
    method=c("GCOMP", "GLMModel", method)
  }
  return(list(
    method=method,
    adj_vars=adj_vars,
    adj_se_z=adj_se_z,
    pu_joint_z=pu_joint_z,
    pu_funcs=pu_funcs,
    omegaz_func=omegaz_func
  ))
}
