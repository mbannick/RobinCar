
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
