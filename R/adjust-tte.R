
# Creates a time-to-event dataset
# to be directly used in adjustment.
create.tte.df <- function(model, data){

  df <- data.frame(
    treat=data$treat,
    response=data$response,
    event=data$event,
    nu_d=nu_d(model$car_scheme)
  )
  # Adjust for x-covariates
  if(model$adj_cov){
    df <- cbind(df, data$covariate)
  }
  # Adjust for z-covariates
  if(model$adj_strata){
    df$carcov_z <- data$joint_strata
  } else {
    df$carcov_z <- 0
  }
  # Covariate-adaptive randomization
  if(model$car_strata){
    df$strata <- data$joint_strata
  } else {
    df$strata <- 0
  }

  return(df)

}

# Pre-process the time to event dataset by calculating
# the risk set size, and fixing ties in the failure times.
process.tte.df <- function(df, ref_arm=NULL){

  # Get the treatment column and set the reference group
  trts <- unique(data$treat)
  if(!is.null(ref_arm)){
    trts <- c(ref_arm, trts[trts!=ref_arm])
  }

  # Calculate risk set sizes at each failure time
  data <- data %>%
    dplyr::arrange(car_strata, response) %>%
    group_by(car_strata) %>%
    mutate(Y=n():1,
           trt0=as.integer(treat == trts[1]),
           trt1=as.integer(treat == trts[2]),
           Y0=cumsum(trt0[n():1])[n():1],
           Y1=cumsum(trt1[n():1])[n():1])

  # Fix ties
  ties <- which(diff(data[, response]) == 0) + 1
  if(length(ties) > 0){
    for(i in ties){
      data[i, c("Y0", "Y1")] <- data[i-1, c("Y0", "Y1")]
    }
  }
  data <- data %>% ungroup()

  return(data)
}
