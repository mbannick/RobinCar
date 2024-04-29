#' @importFrom survival coxph Surv
#' @importFrom stats model.matrix
get.linear.predictor <- function(df, covnames){

  coxdf <- df[, c("response", "event", covnames)]

  # Fit a cox PH model using only the covariates to adjust for
  fit <- survival::coxph(survival::Surv(response, event) ~ ., data=coxdf)
  beta <- fit$coefficients

  # Get the design matrix with only the covariates
  design <- stats::model.matrix(~ 0 + ., data=coxdf[, -c(1, 2)])
  lin_preds <- as.numeric(design %*% beta)

  return(lin_preds)
}

#' @importFrom stats var
#' @importFrom dplyr n
get.car_strata.sum <- function(df, n, p_trt, sparse_remove=FALSE){

  ss <- df %>%
    dplyr::group_by(.data$carcov_z, .drop=TRUE) %>%
    dplyr::summarise(
      cond_var0       = stats::var(.data$O_i[.data$trt0 == 1]),
      cond_var1       = stats::var(.data$O_i[.data$trt1 == 1]),
      cond_mean0      = mean(-.data$O_i[.data$trt0 == 1]),
      cond_mean1      = mean(.data$O_i[.data$trt1 == 1]),
      prob_z          = dplyr::n() / n,
      nu_d            = unique(.data$nu_d),
      .groups = "drop"
    ) %>%
    dplyr::mutate(na_ind=is.na(.data$cond_var0 + .data$cond_var1))

  if(sparse_remove){
    if(any(ss$na_ind)){
      .sparse.car_strata.warn()
      ss <- ss %>% dplyr::filter(!.data$na_ind)
    } else {
      ss <- ss %>% tidyr::replace_na(
        list(cond_var0 = 0, cond_var1 = 0)
      )
    }
  }

  ss <- ss %>%
    dplyr::mutate(
      z_var = .data$prob_z * (p_trt * .data$cond_var0 +
                        (1 - p_trt) * .data$cond_var1 +
                        .data$nu_d * (.data$cond_mean0 + .data$cond_mean1)^2)
    )

  return(ss)
}

#' @exportS3Method
adjust.CoxScore <- function(model, data, ...){

  # Creates data
  df <- create.tte.df(model, data)
  # Risk set and failure time times
  df <- process.tte.df(df, ref_arm=model$ref_arm)

  # If there are covariates to adjust for,
  # fit a Cox model to get linear predictor
  if(model$adj_cov){
    lin_preds <- get.linear.predictor(df,
                                      covnames=colnames(data$covariate))
    df <- df %>% dplyr::mutate(lin_pred=lin_preds)
  } else {
    df <- df %>% dplyr::mutate(lin_pred=0)
  }

  df <- get.ordered.data(df, ref_arm=model$ref_arm)
  strata_sum <- get.car_strata.sum(df,
                               n=data$n,
                               p_trt=model$p_trt,
                               sparse_remove=model$sparse_remove)

  # If there's stratification this calculates the
  # Robust Stratified Score Test (EQ #21)
  numerator <- sum(df$u_i)
  denominator <- sqrt(sum(strata_sum$z_var) * data$n)
  statistic <- numerator / denominator

  result <- list(
    U=numerator,
    se=denominator,
    statistic=statistic
  )

  return(
    structure(
      class="TTEResult",
      list(result=result, settings=model, data=data)
    )
  )
}
