get.linear.predictor <- function(df, covnames){

  coxdf <- df[, c("response", "event", covnames)]

  # Fit a cox PH model using only the covariates to adjust for
  fit <- survival::coxph(survival::Surv(response, event) ~ ., data=coxdf)
  beta <- fit$coefficients

  # Get the design matrix with only the covariates
  design <- model.matrix(~ 0 + ., data=coxdf[, -c(1, 2)])
  lin_preds <- as.numeric(design %*% beta)

  return(lin_preds)
}

get.ordered.data <- function(df, ref_arm){

  df <- df %>%
    group_by(car_strata) %>%
    mutate(
      mu_t            = Y1 / Y,

      O.hat           = event * (trt1 - mu_t) -
                        trt1 * cumsum(event / Y) +
                        cumsum(event * Y1 / Y^2),

      s0_seq          = exp(lin_pred),
      s1_seq          = s0_seq * trt1,

      s0              = cumsum(s0_seq[n():1])[n():1]/n(),
      s1              = cumsum(s1_seq[n():1])[n():1]/n(),

      mu_t            = s1 / s0,

      mean_at_risk    = event / (n()*s0),
      cumsum_at_risk  = cumsum(mean_at_risk),

      mean_at_risk2   = event / (n()*s0) * mu_t,
      cumsum_at_risk2 = cumsum(mean_at_risk2),

      O_i             = event * (trt1 - mu_t) -
                        s1_seq * cumsum_at_risk +
                        s0_seq * cumsum_at_risk2,

      u_i             = event * (trt1 - mu_t)
    ) %>% mutate(
      O.hat           = -O.hat * (trt0 == 1) + O.hat * (trt1 == 1)
    ) %>% ungroup()

  return(df)
}

get.strata.sum <- function(df, sparse_remove=FALSE){

  ss <- df %>%
    group_by(carcov_z, .drop=TRUE) %>%
    summarise(
      cond_var0       = var(O_i[trt0 == 1]),
      cond_var1       = var(O_i[trt1 == 1]),
      cond_mean0      = mean(-O_i[trt0 == 1]),
      cond_mean1      = mean(O_i[trt1 == 1]),
      prob_z          = n() / n,
      nu_d            = unique(nu_d),
      .groups = "drop"
    ) %>%
    mutate(na_ind=is.na(cond_var0 + cond_var1))

  if(sparse_remove){
    if(any(ss$na_ind)){
      .sparse.strata.warn()
      ss <- ss %>% filter(!na_ind)
    } else {
      ss <- ss %>% tidyr::replace_na(
        list(cond_var0 = 0, cond_var1 = 0)
      )
    }
  }

  ss <- ss %>%
    mutate(
      z_var = prob_z * (p_trt * cond_var0 +
                        (1 - p_trt) * cond_var1 +
                        nu_d * (cond_mean0 + cond_mean1)^2)
    )

  return(ss)
}

adjust.CoxScore <- function(model, data){

  # Creates data
  df <- create.tte.df(model, data)
  # Risk set and failure time times
  df <- process.tte.df(model, data)

  # If there are covariates to adjust for,
  # fit a Cox model to get linear predictor
  if(model$adj_cov){
    lin_preds <- get.linear.predictor(df,
                                      covnames=colnames(data$covariate))
    df <- df %>% mutate(lin_pred=lin_preds)
  } else {
    df <- df %>% mutate(lin_pred=0)
  }

  df <- get.ordered.data(df, ref_arm=model$ref_arm)
  strata_sum <- get.strata.sum(df, sparse_remove=model$sparse_remove)

  # If there's stratification this calculates the
  # Robust Stratified Score Test (EQ #21)
  numerator <- sum(df$u_i)
  denominator <- sqrt(sum(strata_sum$z_var) * data$n)
  statistic <- numerator / denominator

  result <- list(
    numerator=numerator,
    denominator=denominator,
    statistic=statistic
  )

  return(
    structure(
      class="CoxScoreResult",
      list(result=result, settings=model, data=data)
    )
  )
}
