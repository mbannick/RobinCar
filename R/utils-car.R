# Compute the Omega(Z) matrix
# based on the given covariate adaptive randomization scheme.
#
# This function returns a function to obtain the correct matrix
# for a specified vector of randomization probabilities pi_t.
omegaz.closure <- function(car_scheme){

  omegaz.func <- function(pi_t){
    # TODO: Check that minimization works out to 0
    if(car_scheme %in% c("simple", "pocock-simon")){
      pi_t <- c(pi_t)
      omegaz <- diag(pi_t) - pi_t %*% t(pi_t)
    } else {
      omegaz <- matrix(0, nrow=length(pi_t), ncol=length(pi_t))
    }
    return(omegaz)
  }
  return(omegaz.func)
}

# Calculate nu_d for the randomization design
nu.d <- function(car_scheme, p_trt=0.5){

  if(length(car_scheme) > 1 | length(car_scheme) == 0){
    nu_d <- NA
  }else{
    if(is.na(car_scheme)){
      nu_d <- NA
    }else{
      if(car_scheme %in% c("permuted-block", "biased-coin")){
        nu_d <- 0
      }else if(car_scheme == "simple"){
        nu_d <- p_trt * (1 - p_trt)
      }
      else{
        nu_d <- NA
      }
    }
  }
  return(nu_d)
}
