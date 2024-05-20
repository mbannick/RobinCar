

test_that("Test contrast function ordering", {

  data <- RobinCar:::data_sim %>%
    # recode A from 0 to "a", 1 to "b", 2 to "c"
    mutate(A=factor(A, levels=0:2, labels=c("b", "c", "a")))

  my.f<-function(x){
    x[2]-x[1] # should be "c"-"b" based on the ordering
  }

  fit.anova <- robincar_linear(df = data,
                               response_col="y",
                               treat_col="A",
                               car_scheme="simple",
                               adj_method="ANOVA",
                               contrast_h=my.f)

  # Check levels are in the right order
  lev <- fit.anova$main$result$treat
  expect_equal(lev, c("b", "c", "a"))

  # Check that the estimates are in correct order
  anova <- mean(data$y[data$A == "c"]) - mean(data$y[data$A == "b"])
  expect_equal(as.vector(fit.anova$contrast$result$estimate[1]), anova)

})
