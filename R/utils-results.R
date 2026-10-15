
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