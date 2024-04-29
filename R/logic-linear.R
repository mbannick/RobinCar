# These are the covariate randomization schemes
# that are applicable; leaves out simple randomization
# and pocock-simon minimization because those are dealt with separately
APPLICABLE.SCHEMES <- c("permuted-block",
                        "biased-coin",
                        "urn")

GENERIC.SETTINGS <- function(method, car_scheme){
  return(list(
    method=method,
    car_scheme=car_scheme,
    omegaz_func=NULL,
    adj_vars=NULL
  ))
}

lmlogic <- function(meth, x_exists, z_exists, car_scheme, cov_strata){
  UseMethod("lmlogic", meth)
}

#' @exportS3Method
lmlogic.ANOVA <- function(meth, x_exists, z_exists, car_scheme, cov_strata){

  s <- GENERIC.SETTINGS("ANOVA", car_scheme)

  if(car_scheme == "simple"){
    if(x_exists) .x.exist.warn()
    if(z_exists) .z.exist.warn()
  }
  if(car_scheme == "pocock-simon"){
    .car.min.err()
  }
  if(car_scheme %in% APPLICABLE.SCHEMES){
    if(z_exists){
      if(x_exists) .x.exist.warn()
      s$omegaz_func <- omegaz.closure(car_scheme)
    } else {
      .z.miss.err()
    }
  }
  return(s)
}

#' @exportS3Method
lmlogic.ANCOVA <- function(meth, x_exists, z_exists, car_scheme, cov_strata){

  s <- GENERIC.SETTINGS("ANCOVA", car_scheme)

  if(car_scheme == "simple"){
    if(z_exists) .z.exist.warn.simple()
    if(x_exists){
      s$adj_vars <- "x"
    } else {
      .x.miss.warn()
      s$method <- "ANOVA"
    }
  }
  if(car_scheme == "pocock-simon"){
    .car.min.err()
  }
  if(car_scheme %in% APPLICABLE.SCHEMES){
    if(z_exists){
      s$omegaz_func <- omegaz.closure(car_scheme)
      if(cov_strata){
        if(x_exists){
          s$adj_vars <- "joint_z_x"
        } else {
          s$adj_vars <- "z"
        }
      } else {
        if(x_exists){
          s$adj_vars <- "x"
        } else {
          s$method <- "ANOVA"
        }
      }
    } else {
      .z.miss.err()
    }
  }
  return(s)
}

#' @exportS3Method
lmlogic.ANHECOVA <- function(meth, x_exists, z_exists, car_scheme, cov_strata){

  s <- GENERIC.SETTINGS("ANHECOVA", car_scheme)

  if(car_scheme == "simple"){
    if(z_exists) .z.exist.warn.simple()
    if(x_exists){
      s$adj_vars <- "x"
    } else {
      .x.miss.warn()
      s$method <- "ANOVA"
    }
  }
  if(car_scheme == "pocock-simon"){
    if(z_exists){
      if(x_exists){
        s$adj_vars <- "joint_z_x"
      } else {
        s$adj_vars <- "joint_z"
      }
    } else {
      .z.miss.err()
    }
  }
  if(car_scheme %in% APPLICABLE.SCHEMES){
    if(z_exists){
      s$omegaz_func <- omegaz.closure(car_scheme)
      if(cov_strata){
        if(x_exists){
          s$adj_vars <- "joint_z_x"
        } else {
          s$adj_vars <- "z"
        }
      } else {
        if(x_exists){
          s$adj_vars <- "x"
        } else {
          s$method <- "ANOVA"
        }
      }
    } else {
      .z.miss.err()
    }
  }
  return(s)
}
