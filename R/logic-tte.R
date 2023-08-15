ALLOWED_COX_CARSCHEME <- c("simple", "permuted-block", "biased-coin")

# @title Logic for performing robincar_coxscore or robincar_logrank for time to event endpoints
# @return warning or error messages based on the user input,
#         whether or not to adjust for covariates, and whether
#         to account for covariate-adaptive randomization
#
# logic_CL <- ttelogic(x_exists=TRUE, z_exists=TRUE,
#                               car_scheme="permuted-block",
#                               adj_method="CL")
#
# logic_CSL <- ttelogic(x_exists=TRUE, z_exists=TRUE,
#                                car_scheme="permuted-block",
#                                adj_method="CSL")
#
# logic_coxscore <- ttelogic(x_exists=TRUE, z_exists=TRUE,
#                                     car_scheme="permuted-block",
#                                     adj_method="coxscore")
ttelogic <- function(x_exists, z_exists,
                     car_scheme, adj_method){

  if(adj_method == "coxscore" & !(car_scheme %in% ALLOWED_COX_CARSCHEME)){
    .cox.carscheme.error()
  }

  if(adj_method %in% c("CL", "CSL", "coxscore")){

    if(x_exists){
      adj_cov <- TRUE
    } else {
      adj_cov <- FALSE
    }
    if(car_scheme == "simple"){
      adj_strata <- FALSE
      if(z_exists){
        .z.exist.warn.simple()
        if(adj_method == "CSL") .csl.nostrata.warn()
      }
    } else {
      if(!z_exists) .z.miss.err()
      adj_strata <- TRUE
    }
  } else{
    .tte.nomethod.error()
  }

  if(adj_method == "CSL" & car_scheme != "simple"){
    car_strata <- TRUE
  } else {
    car_strata <- FALSE
  }

  return(list(
    method=adj_method,
    adj_cov=adj_cov,
    adj_strata=adj_strata,
    car_strata=car_strata
  ))
}
