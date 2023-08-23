
#' @importFrom dplyr mutate group_by summarise arrange
adjust.LogRank <- function(model, data){

  # Creates data
  df <- create.tte.df(model, data)

  # Risk set and failure time times
  df <- process.tte.df(df, ref_arm=model$ref_arm)
  df <- df %>% mutate(lin_pred=0)

  # Get ordered data
  df <- get.ordered.data(df, ref_arm=model$ref_arm)
  df <- get.tte.adjustment(df, model, data)

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
  U_CSL     <- mean(df$uu_cl)
  var_CSL   <- mean(df$ssig_l) - sum(ss$var_adj) / data$n
  se        <- sqrt(var_CSL / data$n)
  statistic <- U_CSL / se

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
