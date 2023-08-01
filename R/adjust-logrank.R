
#' @importFrom dplyr mutate group_by ungroup everything
#' @importFrom stats model.matrix
get.design.matrix <- function(df, covnames){

  formula <- as.formula(
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
  # but this will center the strata variables also, if they are included.
  mat <- data.frame(mat) %>%
    dplyr::mutate(strata=df$strata) %>%
    dplyr::group_by(strata, ) %>%
    dplyr::mutate(across(
      .cols=everything(),
      .fns=list(center=~scale(., center=TRUE, scale=FALSE)),
      .names="{col}_{fn}"
    )) %>%
    dplyr::ungroup() %>%
    subset(select=-strata)
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
    if(all(is.na(df %>% dplyr::filter(grepl("xmat_",term)) %>% .$estimate))){
      .lin.dep.strat.error()
    }
  }

  df <- df %>% dplyr::mutate(
    estimate=tidyr::replace_na(estimate, 0)
  )

  return(df)
}

#' @import broom
#' @importFrom dplyr group_by group_modify ungroup select
#' @importFrom tidyr pivot_wider
regress.to.Ohat <- function(df, stratified){

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
      # strata-specific betas (same)
      strata=unique(df$strata), .) %>%
    dplyr::select(strata,
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

  dat <- dplyr::left_join(df, betas, by="strata")

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
      dplyr::group_by(.data$strata) %>%
      dplyr::group_modify(~ {
        matcx   <- dplyr::select(.x, all_of(covnames_xc)) %>% as.matrix(ncol=p)
        matx    <- dplyr::select(.x, all_of(covnames_x)) %>% as.matrix(ncol=p)
        b1      <- dplyr::select(.x, all_of(beta1_names)) %>% head(n=1) %>% as.matrix(ncol=p)
        b0      <- dplyr::select(.x, all_of(beta0_names)) %>% head(n=1) %>% as.matrix(ncol=p)
        adjust1 <- matcx %*% t(b1)
        adjust0 <- matcx %*% t(b0)
        bsigb   <- as.numeric((b1+b0) %*% var(matx) %*% t(b1 + b0))
        out     <- data.frame(adjust1, adjust0, bsigb) %>% as_tibble()
        cbind(.x,out)
      })
  } else {
    dat <- dat %>% mutate(
      adjust1 = 0,
      adjust0 = 0,
      bsigb   = 0
    )
  }
  return(dat)
}

#' @importFrom dplyr mutate group_by summarise arrange
adjust.LogRank <- function(model, data){

  # Creates data
  df <- create.tte.df(model, data)

  # Risk set and failure time times
  df <- process.tte.df(df, ref_arm=model$ref_arm)
  df <- df %>% mutate(lin_pred=0)

  # Get ordered data
  df <- get.ordered.data(df, ref_arm=model$ref_arm)

  if(model$adj_strata | model$adj_cov){
    # Add covariate names in order
    # with the strata covariates first
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
    xmat <- get.design.matrix(df, covnames)
    df <- cbind(df, xmat)

    # Regress to Ohat -- get betas, then calculate adjustment using betas
    betas <- regress.to.Ohat(df, stratified=(model$method == "CSL"))

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

  df <- df %>%
    dplyr::mutate(
      uu_cl  = .data$trt1 * (.data$O.hat - .data$adjust1) -
               .data$trt0 * (.data$O.hat - .data$adjust0),
      ssig_l = .data$event * .data$Y0 * .data$Y1 / .data$Y^2
    )

  # Summarize by strata (if CL, then single strata)
  ss <- df %>%
    dplyr::filter(!is.na(.data$uu_cl)) %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(
      U_SL_z  = sum(.data$uu_cl),
      var_adj = model$p_trt * (1 - model$p_trt) * unique(.data$bsigb) * n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(strata)

  # Final quantities for the C(S)L statistic
  U_CSL = mean(df$uu_cl)
  var_CSL = mean(df$ssig_l) - sum(ss$var_adj) / data$n
  se = sqrt(var_CSL / data$n)
  statistic = U_CSL / se

  result <- list(
    strata_sum=ss,
    U=U_CSL,
    se=se,
    statistic=statistic
  )

  return(
    structure(
      class="TTEResult",
      list(result=result, settings=model, data=data)
    )
  )
}
