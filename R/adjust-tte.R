
# Creates a time-to-event dataset
# to be directly used in adjustment.
create.tte.df <- function(model, data){

  df <- data.frame(
    treat=data$treat,
    response=data$response,
    event=data$event,
    nu_d=nu.d(model$car_scheme)
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

fix.ties <- function(df){
  ties <- which(base::diff(df$response) == 0) + 1
  for(i in ties){
    df[i, c("Y0", "Y1")] <- df[i-1, c("Y0", "Y1")]
  }
  return(df)
}

# Pre-process the time to event dataset by calculating
# the risk set size, and fixing ties in the failure times.
process.tte.df <- function(df, ref_arm=NULL){

  # Get the treatment column and set the reference group
  trts <- levels(df$treat)
  if(!is.null(ref_arm)){
    index <- which(trts == ref_arm)
    trts <- c(ref_arm, trts[-index])
  }

  # Calculate risk set sizes at each failure time
  df <- df %>%
    dplyr::arrange(.data$strata, .data$response) %>%
    group_by(.data$strata) %>%
    mutate(Y=n():1,
           trt0=as.integer(.data$treat == trts[1]),
           trt1=as.integer(.data$treat == trts[2]),
           Y0=cumsum(.data$trt0[n():1])[n():1],
           Y1=cumsum(.data$trt1[n():1])[n():1])

  # Fix ties
  df <- fix.ties(df)
  df <- df %>% ungroup()

  return(df)
}

#' @import dplyr
get.ordered.data <- function(df, ref_arm){

  df <- df %>%
    group_by(.data$strata) %>%
    mutate(
      mu_t            = .data$Y1 / .data$Y,

      O.hat           = .data$event * (.data$trt1 - .data$mu_t) -
        .data$trt1 * cumsum(.data$event / .data$Y) +
        cumsum(.data$event * .data$Y1 / .data$Y^2),

      s0_seq          = exp(.data$lin_pred),
      s1_seq          = .data$s0_seq * .data$trt1,

      s0              = cumsum(.data$s0_seq[n():1])[n():1]/n(),
      s1              = cumsum(.data$s1_seq[n():1])[n():1]/n(),

      mu_t            = .data$s1 / .data$s0,

      mean_at_risk    = .data$event / (n()*.data$s0),
      cumsum_at_risk  = cumsum(.data$mean_at_risk),

      mean_at_risk2   = .data$event / (n()*.data$s0) * .data$mu_t,
      cumsum_at_risk2 = cumsum(.data$mean_at_risk2),

      O_i             = .data$event * (.data$trt1 - .data$mu_t) -
        .data$s1_seq * .data$cumsum_at_risk +
        .data$s0_seq * .data$cumsum_at_risk2,

      u_i             = .data$event * (.data$trt1 - .data$mu_t)
    ) %>% mutate(
      O.hat           = -.data$O.hat * (.data$trt0 == 1) +
        .data$O.hat * (.data$trt1 == 1)
    ) %>% ungroup()

  return(df)
}
