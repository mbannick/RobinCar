
#' @importFrom dplyr mutate group_by summarise arrange n
#' @exportS3Method
adjust.LogRank <- function(model, data, ...){

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
      ssig_l = .data$event * .data$Y0 * .data$Y1 * (.data$Y - .data$n.events) / .data$Y^2 / (.data$Y-1)
    )

  # Summarize by car_strata (if CL, then single car_strata)
  ss <- df %>%
    dplyr::filter(!is.na(.data$uu_cl)) %>%
    dplyr::group_by(.data$car_strata) %>%
    dplyr::summarise(
      U_SL_z  = sum(.data$uu_cl),
      var_adj = model$p_trt * (1 - model$p_trt) * unique(.data$bsigb) * dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$car_strata)

  # Final quantities for the C(S)L statistic
  U_CSL     <- sum(df$uu_cl) / data$n
  var_CSL   <- (sum(df$ssig_l, na.rm=TRUE) - sum(ss$var_adj)) / data$n^2
  se        <- sqrt(var_CSL)
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
