# These are the covariate randomization schemes
# that are applicable; leaves out simple randomization
# and pocock-simon minimization because those are dealt with separately
APPLICABLE.SCHEMES <- c("permuted-block",
                        "biased-coin",
                        "urn")

glmlogic <- function(z_exists, car_scheme){

  # PREDICTION UNBIASEDNESS -------------------------
  # Whether to check prediction un-biasedness in all joint levels
  # of Z, or just overall, and whether to report an error or warning
  pu_joint_z <- FALSE

  # Biased action will be one of "", "warning", "error"
  pu_funcs <- function(){}

  if(z_exists & car_scheme == "pocock-simon"){

    pu_funcs <- c(.pu.z.calibrate, .pu.z.warn)
    pu_joint_z <- TRUE

  } else {

    pu_funcs <- .pu.warn
    pu_joint_z <- FALSE

  }

  # CAR OMEGA FUNCTION -----------------------------
  # Get the Omega Z function
  omegaz_func <- NULL

  if(z_exists & (car_scheme != "simple")){
    omegaz_func <- omegaz.closure(car_scheme)
  }

  return(list(
    method=c("GLMModel"),
    pu_joint_z=pu_joint_z,
    pu_funcs=pu_funcs,
    omegaz_func=omegaz_func
  ))
}
