SLlogic <- function(x_exists, z_exists, car_scheme, cov_strata){

  adj_vars <- NULL

  aipw <- TRUE
  adj_se_z <- FALSE

  # Whether to check prediction un-biasedness in all joint levels
  # of Z, or just overall, and whether to report an error or warning
  pu_joint_z <- FALSE
  # Biased action will be one of "", "warning", "error"
  pu_funcs <- function(){}

  omegaz_func <- NULL

  if(car_scheme == "simple"){
    if(z_exists) .z.exist.warn.simple()
    if(x_exists){
      adj_vars <- "x"
      pu_funcs <- .pu.warn
    } else {
      .x.miss.err()
    }
  }
  if(car_scheme == "pocock-simon"){
    if(!z_exists){
      .z.miss.err()
    } else {
      omegaz_func <- omegaz.closure(car_scheme)
      if(!cov_strata){
        .z.include.warn()
      }
      pu_joint_z <- TRUE
      if(x_exists){
        adj_vars <- "joint_z_x"
        pu_funcs <- .pu.z.warn
      } else {
        adj_vars <- "joint_z"
        pu_funcs <- c(.pu.z.calibrate, .pu.z.warn)
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
        }
      }
    } else {
      .z.miss.err()
    }
  }
  method=c("SLModel")
  return(list(
    method=method,
    adj_vars=adj_vars,
    adj_se_z=adj_se_z,
    pu_joint_z=pu_joint_z,
    pu_funcs=pu_funcs,
    omegaz_func=omegaz_func
  ))
}
