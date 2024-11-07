library(MASS)
library(dplyr)

DATA <- RobinCar:::data_sim
DATA$A <- as.factor(DATA$A)
DATA$y_bin <- ifelse(DATA$y > 2, 1, 0)

# Ensure that the naming/capitalization
# of the treatment variable does not impact results.

DATA <- DATA %>% mutate(
  A2=paste0("TREAT", A),
  A3=ifelse(A2 == "TREAT1", "treat1", A2)
)

# Doing DATA %>% group_by(A3) %>% summarize(n()) will
# order the levels by c("TREAT0", "TREAT2", "treat1")
#
# whereas factor(DATA$A3) will order the levels
# by c("TREAT0", "TREAT1", "TREAT2").

test_that("Test ordering of treatment factors in GLM adjustment", {

  # The locale is actually the issue -- testthat
  # automatically sets locale C but the issue only
  # appears under locale != C. So here we set it to
  # en_US.UTF-8.

  original_locale <- Sys.getlocale("LC_COLLATE")
  Sys.setlocale("LC_COLLATE", "en_US.UTF-8")

  # English language names
  langENG <- robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A2",
    car_scheme="simple",
    formula="y ~ A2 + x1")

  # C language names
  langC <- robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A3",
    car_scheme="simple",
    formula="y ~ A3 + x1")

  # Check that both of the orderings worked,
  # regardless of the final order the output is in
  # -----------------------------------------

  # Benchmark ordering
  namesENG <- c("TREAT0", "TREAT1", "TREAT2")
  namesC <- c("TREAT0", "treat1", "TREAT2")

  # Current ordering
  currentENG <- langENG$result$treat
  currentC <- langC$result$treat

  # Get re-ordering indexes
  reorderENG <- unname(sapply(currentENG, function(x) which(x == namesENG)))
  reorderC <- unname(sapply(currentC, function(x) which(x == namesC)))

  # Re-order estimate output and make sure they're the same
  estENG <- unname(langENG$result$estimate[reorderENG])
  estC <- unname(langC$result$estimate[reorderC])
  expect_equal(estENG, estC)

  # Re-order variance output and make sure they're the same
  varENG <- unname(langENG$varcov[reorderENG, reorderENG])
  varC <- unname(langC$varcov[reorderC, reorderC])
  expect_equal(varENG, varC)

  # Reset the locale
  Sys.setlocale("LC_COLLATE", original_locale)

})

