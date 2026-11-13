# Function to format results.
format_results <- function(labels, estimates, varcov, label_name="treat", ...){

  # Extract estimates and create results data
  result <- dplyr::tibble(
    label=labels,
    estimate=c(estimates),
    se=diag(varcov**0.5)
  )
  colnames(result) <- c(label_name, "estimate", "se")

  # Compute p-values based on the correct variance estimates
  result <- result %>% dplyr::mutate(
    `pval (2-sided)`=2*stats::pnorm(abs(.data$estimate/.data$se), lower.tail=F)
  )
  return(result)
}



## function to apply finite sample correction
apply_correction <- function(result){
  
  if (!(result$settings$method %in% c("CL", "CSL"))) return(result)  
  
  n <- result$data$n
  n_0 <- n * result$data$pie[1,1]
  n_1 <- n * result$data$pie[2,1]
  
  ### work out how many terms come from the covariate columns...
  form <- paste0(result$data$response_col, " ~ ", paste(c("1", result$data$covariate_cols), collapse = "+"))
  fit_lin <- lm(as.formula(form), data = result$original_df)
  n_terms_cov <- length(coef(fit_lin)) - 1
  
  ## add on the covariate strata terms (if needed)
  n_terms_strata <- length(result$data$joint_strata_levels) - 1
  
  n_terms <- ifelse(result$settings$method == "CL" && result$settings$car_scheme != "simple",
                    n_terms_cov + n_terms_strata,
                    n_terms_cov)
  
  ## formula from Krotka & Magirr
  correction <- (1 + n_terms * (((n_0 / n)^2) / (n_1 - n_terms - 2) + ((n_1 / n)^2) / (n_0 - n_terms - 2))) *
    (n / (n - ((n_0 / n_1) + (n_1 / n_0) - 1) * n_terms - ((n_0 / n_1) + (n_1 / n_0))))
  
  correction <- unname(correction)
  
  ## what to do depends on estimation vs testing
  if (attr(result$settings, "class")[1] == "CovHR"){
    
    result$result$se_theta_CL <- result$result$se_theta_CL * sqrt(correction)
    
  }
  else {
    
    result$result$se <- result$result$se * sqrt(correction)
    result$result$statistic <- result$result$U / result$result$se
    
  }
  
  return(result)
  
}


# ---------------------------------------------------------------------------
# HC3 variance correction
# ---------------------------------------------------------------------------

# Compute the three additive HC3 correction terms for the variance of the
# covariate-adjusted log-rank statistic or log-HR estimator.
#
# Terms:
#   B1 -- variance from estimating A_n' delta_hat
#   B2 -- variance from estimating Sigma_X in the denominator
#   B3 -- covariance between U_oracle and the estimation error (currently not used)
#
# Returns zeros for all three terms when there are no covariates.
#
# @param df    Data frame produced inside adjust.LogRank / adjust.CovHR,
#              after get.tte.adjustment() has been called.  Must contain
#              O.hat, trt0, trt1, car_strata, and the centered covariate
#              columns (xmat_*_center).
# @param p_trt Randomization probability for the treatment arm.
# @param n     Total sample size.
#
# @return List with numeric elements B1, B2, B3.
compute_hc3_correction <- function(df, p_trt, n) {

  # Centered covariate columns produced by get_design_matrix().
  # Stratum-indicator columns (carcov_z) are excluded because they are
  # absorbed into the stratification and must not enter the HC3 sandwich.
  covnames_c  <- grep("^xmat_(?!.*carcov_z).*_center$", colnames(df),
                      perl = TRUE, value = TRUE)
  covnames_uc <- gsub("_center$", "", covnames_c)

  if (length(covnames_c) == 0L) return(list(B1 = 0, B2 = 0, B3 = 0))

  strata <- df$car_strata

  # A_n = (1/n) * sum_{i : trt=1} (X_i - Xbar_stratum).
  # covnames_c are the stratum-centred columns produced by get_design_matrix(),
  # so summing them within trt=1 gives A_n directly.
  An <- colSums(as.matrix(df[df$trt1 == 1L, covnames_c, drop = FALSE])) / n

  # SigmaX: stratum-size-weighted mean of within-stratum sample covariances,
  # computed on the stratum-centred columns.
  SigmaX <- Reduce("+", lapply(unique(strata), function(z) {
    Xz <- as.matrix(df[strata == z, covnames_c, drop = FALSE])
    nz <- nrow(Xz)
    (nz / n) * var(as.matrix(df[strata == z, covnames_c, drop = FALSE]))
  }))

  # Per-arm HC3 sandwich: S_j^{-1} * meat_j * S_j^{-1},
  # where meat_j = sum_i Xtilde_i Xtilde_i' (e_i / (1 - h_i))^2.
  #
  # Covariates must be centered within arm and within stratum (not across arms),
  arm_sandwich <- function(arm_df) {
    # Re-centre within arm, within stratum
    Xraw <- as.matrix(arm_df[, covnames_uc, drop = FALSE])
    arm_strata <- arm_df$car_strata
    Xc <- Xraw
    for (z in unique(arm_strata)) {
      idx <- arm_strata == z
      Xc[idx, ] <- scale(Xraw[idx, , drop = FALSE], center = TRUE, scale = FALSE)
    }
    O    <- arm_df$O.hat
    S    <- crossprod(Xc)
    Sinv <- solve(S)
    gam  <- as.numeric(Sinv %*% crossprod(Xc, O))
    res  <- O - as.numeric(Xc %*% gam)
    lev  <- rowSums((Xc %*% Sinv) * Xc)      # hat-matrix diagonal
    meat <- crossprod(Xc * (res / (1 - lev))) # HC3 meat
    list(Sinv = Sinv, Xc = Xc, res = res, lev = lev,
         sandwich = Sinv %*% meat %*% Sinv)
  }

  arm0  <- arm_sandwich(df[df$trt0 == 1L, , drop = FALSE])
  arm1  <- arm_sandwich(df[df$trt1 == 1L, , drop = FALSE])
  V_HC3 <- arm0$sandwich + arm1$sandwich

  B1 <- as.numeric(t(An) %*% V_HC3 %*% An)
  B2 <- p_trt * (1 - p_trt) * sum(diag(SigmaX %*% V_HC3)) / n

  # B3 = (-2/n) * An' * (S1inv sum_i v_i1 Xtilde_i1
  #                    - S0inv sum_i v_i0 Xtilde_i0),
  # where v_ij = e_ij^2 / (1 - h_ij)^2.
  v1      <- arm1$res^2 / (1 - arm1$lev)^2
  v0      <- arm0$res^2 / (1 - arm0$lev)^2
  B3_vec  <- arm1$Sinv %*% colSums(arm1$Xc * v1) -
             arm0$Sinv %*% colSums(arm0$Xc * v0)
  B3 <- as.numeric(-2 / n * t(An) %*% B3_vec)

  list(B1 = B1, B2 = B2, B3 = B3)
}


# Apply HC3 variance correction to a TTEResult or TTEResultEst.
#
# Overwrites se (log-rank) or se_theta_CL (CovHR) in place, matching the
# convention used by apply_correction() for the normal approximation.
#
# @param result              Object returned by adjust.LogRank or adjust.CovHR,
#                            with original_df already attached.
# @param df                  Processed data frame from inside the adjust.*
#                            method (contains O.hat and xmat_* columns).
# @param p_trt               Randomisation probability for the treatment arm.
# @param n                   Total sample size.
# @param sig2_L              (CovHR only) unadjusted information quantity.
# @param sig2_CL             (CovHR only) adjusted information quantity.
apply_hc3_correction <- function(result,
                                 df, p_trt, n,
                                 sig2_L = NULL, sig2_CL = NULL) {

  hc3    <- compute_hc3_correction(df, p_trt, n)
  extra  <- hc3$B1 + hc3$B2

  is_covhr <- inherits(result, "TTEResultEst")

  if (is_covhr) {
    result$result$se_theta_CL <- sqrt(sig2_CL / sig2_L**2 / n + extra / sig2_L**2)
  } else {
    var_CSL <- result$result$se^2
    result$result$se        <- sqrt(var_CSL + extra)
    result$result$statistic <- result$result$U / result$result$se
  }

  return(result)
}
