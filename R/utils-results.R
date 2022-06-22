
# Function to format results.

format.results <- function(labels, estimates, varcov, label_name="treat"){

  # Extract estimates and create results data
  result <- tibble(
    label=labels,
    estimate=c(estimates),
    se=diag(varcov**0.5)
  )
  colnames(result) <- c(label_name, "estimate", "se")

  # Compute p-values based on the correct variance estimates
  result <- result %>% mutate(
    `pval (2-sided)`=2*pnorm(abs(estimate/se), lower.tail=F)
  )
  return(result)
}
