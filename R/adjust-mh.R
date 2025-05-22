#' @importFrom stats weighted.mean
estimate.mh <- function(df.mh) {
  delta <- df.mh$y1.k - df.mh$y0.k
  weight <- df.mh$n.1k * df.mh$n.0k / (df.mh$n.1k + df.mh$n.0k)
  delta_all <- weighted.mean(delta, weight)
  
  return(delta_all)
}

var.mh.GR <- function(n11k, n10k, n.1k, n.0k) {
  weight <- n.1k * n.0k / (n.1k + n.0k)
  var.k <- (n11k * (n.1k - n11k) * n.0k^3 + n10k * (n.0k - n10k) * n.1k^3) / n.1k / n.0k / (n.1k + n.0k)^2
  var.est.GR <- sum(var.k[weight != 0]) / (sum(weight)^2)
  return(var.est.GR)
}



var.mh.mGR <- function(n11k, n10k, n.1k, n.0k){
  help_f <- Vectorize(function(x) {
    if (x > 1) {
      return(x / (x - 1))
    } else {
      1
    }
  })
  
  n01k <- n.1k - n11k
  n00k <- n.0k - n10k
  weight <- n.1k * n.0k / (n.1k + n.0k)
  var.k <- weight^2 * (n11k * n01k / n.1k^3 * help_f(n.1k) + n10k * n00k /
                         n.0k^3 * help_f(n.0k))
  var.est.mGR <- sum(var.k[weight != 0]) / (sum(weight)^2)
  return(var.est.mGR)
}


var.mh.Sato <- function(n11k, n10k, n.1k, n.0k){
  weight <- n.1k * n.0k / (n.1k + n.0k)
  Pk <- (n.1k^2 * n10k - n.0k^2 * n11k + n.1k * n.0k * (n.0k - n.1k) / 2) / (n.1k + n.0k)^2
  Qk <- (n11k * (n.0k - n10k) + n10k * (n.1k - n11k)) / 2 / (n.1k  + n.0k)
  delta <- n11k / n.1k - n10k / n.0k
  delta_all <- weighted.mean(delta, weight)
  var.est.Sato <- (delta_all * sum(Pk) + sum(Qk)) / (sum(weight)^2)
  return(var.est.Sato)
}

var.ate.nu <- function(df, ATE.est){
  n.1k <- df$n.1k
  n.0k <- df$n.0k
  n11k <- df$n11k
  n10k <- df$n10k
  y1.k <- df$y1.k
  y0.k <- df$y0.k
  y1.var.k <- df$y1.var.k
  y0.var.k <- df$y0.var.k
  
  K <- length(n11k)
  n.k <- n.1k + n.0k
  n <- sum(n.k)
  pi1 <- sum(n.1k) / n
  pi0 <- sum(n.0k) / n
  rho.hat <- n.k/n
  delta.k <- n11k / n.1k - n10k / n.0k
  weights <- n.1k*n.0k / (n.1k + n.0k)
  
  p1k <- n11k / n.1k
  p0k <- n10k / n.0k
  
  p1k.sq <- p1k^2 - ifelse(n.1k>1,y1.var.k,0) / n.1k
  p0k.sq <- p0k^2 - ifelse(n.0k>1,y0.var.k,0) / n.0k
  
  delta.k.sq <- p1k.sq - 2*p1k * p0k + p0k.sq
  tmp1 <- rho.hat * (delta.k.sq - ATE.est^2)
  tmp2 <- (delta.k.sq - 2 * delta.k * ATE.est + ATE.est^2) * pi1 * pi0 *
    (n.k - 1) / n.k * (n.k - 1 - (4 * n.k - 6) * pi1 * pi0) / n
  var.est.nu <- sum(tmp1[weights!=0])*pi1^2*pi0^2 + sum(tmp2[weights!=0])
  
  var.est.nu <- var.est.nu / n / (sum(weights)/n)^2
  
  return(var.est.nu)
}

est.var.mh <- function(df.mh, estimand, ci_type){
  mh_est <- estimate.mh(df.mh)
  var_est <- with(df.mh,
                  switch(ci_type,
                         "GR"    = var.mh.GR(n11k, n10k, n.1k, n.0k),
                         "mGR"   = var.mh.mGR(n11k, n10k, n.1k, n.0k) + (estimand=="ATE") * var.ate.nu(df.mh, mh_est),
                         "Sato"  = var.mh.Sato(n11k, n10k, n.1k, n.0k),
                         stop("Invalid ci_type. Choose from 'GR', 'mGR', or 'Sato'."))
  )
  return(list(estimate = mh_est, var_est = var_est))
}

# Estimate the MH risk difference or ATE
# contrast is incorporated in adjust.MHModel as we don't have estimate for each treatment.
#' @importFrom utils combn
#' @exportS3Method
adjust.MHModel <- function(model, data, ...){
  
  estimand <- model$estimand
  ci_type <- model$ci_type
  
  treat_col <- data$treat_col
  response_col <- data$response_col
  strata <- data$joint_strata
  
  treat_levels <- data$treat_levels
  n_treat <- length(treat_levels)
  
  strata_list <- split(data$df, strata)
  
  n.k <- matrix(0, nrow = length(strata_list), ncol = n_treat)  # Total in each group
  n1.k <- matrix(0, nrow = length(strata_list), ncol = n_treat) # Events in each group
  y.var.k <- matrix(0, nrow = length(strata_list), ncol = n_treat) # Events in each group
  
  for (i in seq_along(strata_list)) {
    df_k <- strata_list[[i]]
    
    for (j in 1:n_treat) {
      t_level <- treat_levels[j]
      n.k[i, j] <- sum(df_k[[treat_col]] == t_level)   # Total subjects in treatment j
      n1.k[i, j] <- sum(df_k[[response_col]][df_k[[treat_col]] == t_level])  # Events in treatment j
      y.var.k[i, j] <- var(df_k[[response_col]][df_k[[treat_col]] == t_level])
    }
  }
  
  risk_diffs <- matrix(0, nrow = n_treat, ncol = n_treat)
  var_mh <- matrix(0, nrow = n_treat, ncol = n_treat)
  
  # Reserved for multiple treatment levels
  for (j in 1:(n_treat-1)) {
    for (k in (j+1):n_treat) {
        df_mh <- data.frame(
          n.1k = n.k[, j],
          n.0k = n.k[, k],
          n11k = n1.k[, j],
          n10k = n1.k[, k],
          y1.k = n1.k[, j] / n.k[, j],
          y0.k = n1.k[, k] / n.k[, k],
          y1.var.k = y.var.k[, j],
          y0.var.k = y.var.k[, k]
        )
        
        # Compute MH estimate and variance
        mh_result <- est.var.mh(df_mh, estimand = estimand, ci_type = ci_type)
        
        risk_diffs[j, k] <- - mh_result$estimate
        var_mh[j, k] <- mh_result$var_est
      
    }
  }
  
  name <- function(tl) combn(tl, 2, function(x) paste0("treat ", x[2], " - ", x[1]))
  
  lab <- name(treat_levels)
    
  c_est <- risk_diffs[upper.tri(risk_diffs)]
  
  vcv_ncol <- (n_treat^2 - n_treat)/2
  vcv <- diag(x = var_mh[upper.tri(var_mh)], nrow = vcv_ncol, ncol = vcv_ncol)
  result <- format_results(lab, c_est, vcv, label_name = "contrast")
  
  return(structure(class = "MHResult", list(result = result, settings = model)))
}
