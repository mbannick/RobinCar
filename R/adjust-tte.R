
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
    df$car_strata <- data$joint_strata
  } else {
    df$car_strata <- 0
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
#' @importFrom dplyr n
process.tte.df <- function(df, ref_arm=NULL){

  # Get the treatment column and set the reference group
  trts <- levels(df$treat)
  if(!is.null(ref_arm)){
    index <- which(trts == ref_arm)
    trts <- c(ref_arm, trts[-index])
  }

  # Calculate risk set sizes at each failure time
  df <- df %>%
    dplyr::arrange(.data$car_strata, .data$response) %>%
    group_by(.data$car_strata) %>%
    mutate(Y=dplyr::n():1,
           trt0=as.integer(.data$treat == trts[1]),
           trt1=as.integer(.data$treat == trts[2]),
           Y0=cumsum(.data$trt0[dplyr::n():1])[dplyr::n():1],
           Y1=cumsum(.data$trt1[dplyr::n():1])[dplyr::n():1])

  # Fix ties
  df <- fix.ties(df)
  df <- df %>% ungroup()

  return(df)
}

#' @importFrom dplyr group_by mutate ungroup n
get.ordered.data <- function(df, ref_arm){

  df <- df %>%
    group_by(.data$car_strata) %>%
    mutate(
      mu_t            = .data$Y1 / .data$Y,

      O.hat           = .data$event * (.data$trt1 - .data$mu_t) -
        .data$trt1 * cumsum(.data$event / .data$Y) +
        cumsum(.data$event * .data$Y1 / .data$Y^2),

      s0_seq          = exp(.data$lin_pred),
      s1_seq          = .data$s0_seq * .data$trt1,

      s0              = cumsum(.data$s0_seq[dplyr::n():1])[dplyr::n():1]/dplyr::n(),
      s1              = cumsum(.data$s1_seq[dplyr::n():1])[dplyr::n():1]/dplyr::n(),

      mu_t            = .data$s1 / .data$s0,

      mean_at_risk    = .data$event / (dplyr::n()*.data$s0),
      cumsum_at_risk  = cumsum(.data$mean_at_risk),

      mean_at_risk2   = .data$event / (dplyr::n()*.data$s0) * .data$mu_t,
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

#' @importFrom dplyr mutate group_by ungroup everything
#' @importFrom stats model.matrix
get_design_matrix <- function(df, covnames){

  formula <- stats::as.formula(
    paste0("~ 1 + ", paste0(covnames, collapse="+"))
  )
  # Get design matrix based on formula
  # remove intercept column because will use centered covariates
  mat <- stats::model.matrix(formula, df)[,-1]
  mat <- data.frame(mat)
  colnames(mat) <- paste0("xmat_", colnames(mat))

  # Center each column of the design matrix
  # This shouldn't do anything for the covariates, because
  # they should already be centered in the .make.data function,
  # but this will center the car_strata variables also, if they are included.
  mat <- data.frame(mat) %>%
    dplyr::mutate(car_strata=df$car_strata) %>%
    dplyr::group_by(.data$car_strata, ) %>%
    dplyr::mutate(dplyr::across(
      .cols=everything(),
      .fns=list(center=~scale(., center=TRUE, scale=FALSE)),
      .names="{col}_{fn}"
    )) %>%
    dplyr::ungroup() %>%
    subset(select=-car_strata)
  return(mat)
}

#' @importFrom dplyr mutate
#' @importFrom tidyr replace_na
check.collinearity <- function(df, covnames, stratified){

  if(!stratified){
    if(anyNA(df$estimate)){
      .lin.dep.error()
    }
  } else {
    names <- unique(df$term)
    x_covnames <- names[grepl("^xmat_", names)]
    x_covnames <- x_covnames[!grepl("carcov_z", x_covnames)]
    if(all(is.na(df %>% dplyr::filter(grepl("xmat_", .data$term)) %>% .$estimate))){
      .lin.dep.strat.error()
    }
  }

  df <- df %>% dplyr::mutate(
    estimate=tidyr::replace_na(.data$estimate, 0)
  )

  return(df)
}

#' @import broom
#' @importFrom dplyr group_by group_modify ungroup select
#' @importFrom tidyr pivot_wider
regress_to_Ohat <- function(df, stratified){

  # Get design matrix covariate variables, using uncentered variables
  covnames <- colnames(df)[grepl("^xmat_", colnames(df))]
  covnames <- covnames[!grepl("_center", covnames)]

  # Perform a linear regression by treatment group
  res <- df %>%
    dplyr::group_by(.data$trt1) %>%
    dplyr::group_modify(~broom::tidy(
      lm(O.hat ~ 1 + ., data=.x %>%
           select("O.hat", covnames))
    ))
  res <- check.collinearity(res, covnames, stratified)

  # Extract the coefficients
  betas <- res %>%
    dplyr::ungroup() %>%
    tidyr::expand_grid(
      # car_strata-specific betas (same)
      car_strata=unique(df$car_strata), .) %>%
    dplyr::select(.data$car_strata,
                  .data$trt1,
                  .data$term,
                  .data$estimate) %>%
    tidyr::pivot_wider(
      names_from=c(.data$trt1, .data$term),
      names_prefix="trt1_",
      values_from=.data$estimate
    )

  return(betas)
}

#' @importFrom dplyr left_join mutate group_by group_modify
calculate.adjustment <- function(df, betas, covnames, stratified){

  dat <- dplyr::left_join(df, betas, by="car_strata")

  if(stratified){
    covnames_carcov <- grepl("carcov_z", covnames)
    covnames_use <- c(covnames)[!covnames_carcov]
  } else {
    covnames_use <- covnames
  }

  # If there are some covariates:
  if(length(covnames_use > 0)){
    covnames_x <- paste0("xmat_", covnames_use)
    covnames_xc <- paste0("xmat_", covnames_use, "_center")
    beta1_names <- paste0("trt1_1_xmat_", covnames_use)
    beta0_names <- paste0("trt1_0_xmat_", covnames_use)

    dat <- dat %>%
      dplyr::group_by(.data$car_strata) %>%
      dplyr::group_modify(~ {
        matcx   <- dplyr::select(.x, all_of(covnames_xc)) %>% as.matrix(ncol=p)
        matx    <- dplyr::select(.x, all_of(covnames_x)) %>% as.matrix(ncol=p)
        b1      <- dplyr::select(.x, all_of(beta1_names)) %>% head(n=1) %>% as.matrix(ncol=p)
        b0      <- dplyr::select(.x, all_of(beta0_names)) %>% head(n=1) %>% as.matrix(ncol=p)
        adjust1 <- matcx %*% t(b1)
        adjust0 <- matcx %*% t(b0)
        bsigb   <- as.numeric((b1+b0) %*% var(matx) %*% t(b1 + b0))
        out     <- data.frame(adjust1, adjust0, bsigb) %>% dplyr::as_tibble()
        cbind(.x,out)
      })
  } else {
    dat <- dat %>% dplyr::mutate(
      adjust1 = 0,
      adjust0 = 0,
      bsigb   = 0
    )
  }
  return(dat)
}

get.tte.adjustment <- function(df, model, data){
  if(model$adj_strata | model$adj_cov){
    # Add covariate names in order
    # with the car_strata covariates first
    # so that if there is collinearity, they
    # are chosen first over the x covariates.
    covnames <- c()
    if(model$adj_strata){
      covnames <- c(covnames, "carcov_z")
    }
    if(model$adj_cov){
      covnames <- c(covnames, names(data$covariate))
    }
    # Get design matrix and centered versions of variables
    xmat <- get_design_matrix(df, covnames)
    df <- cbind(df, xmat)

    # Regress to Ohat -- get betas, then calculate adjustment using betas
    betas <- regress_to_Ohat(df, stratified=(model$method == "CSL"))

    # Get new covariate names based on the design matrix (helpful for factors)
    new_covnames <- colnames(xmat)
    new_covnames <- gsub("xmat_", "", new_covnames)
    new_covnames <- gsub("_center", "", new_covnames)
    new_covnames <- unique(new_covnames)

    # Calculate adjustment using betas
    df <- calculate.adjustment(df, betas,
                               covnames=new_covnames,
                               stratified=(model$method == "CSL"))
  } else {
    df <- df %>%
      dplyr::mutate(
        adjust1 = 0,
        adjust0 = 0,
        bsigb   = 0
      )
  }
  return(df)
}
