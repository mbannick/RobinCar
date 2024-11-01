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

test_that("GLM full function -- linear (ANOVA)", {

  # The locale is actually the issue -- testthat
  # automatically sets locale C but the issue only
  # appears under locale != C. So here we set it to
  # en_US.UTF-8.

  original_locale <- Sys.getlocale("LC_COLLATE")
  Sys.setlocale("LC_COLLATE", "en_US.UTF-8")

  # English language ordering
  langENG <- robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A2",
    car_scheme="simple",
    formula="y ~ A2 + x1")

  # C language ordering
  langC <- robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A3",
    car_scheme="simple",
    formula="y ~ A3 + x1")

  expect_equal(langENG$result$estimate, langC$result$estimate)
  expect_equal(c(langENG$varcov), c(langC$varcov))

  # Reset the locale
  Sys.setlocale("LC_COLLATE", original_locale)

})

