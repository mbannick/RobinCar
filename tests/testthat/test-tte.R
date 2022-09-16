library(survival)

test_that("No X no Z case under SR, CL = logrank test", {
  set.seed(0)
  n <- 100
  data.simu0 <- data_gen(
    n=n,
    theta=0,
    randomization="permuted_block",
    p_trt=0.5,
    case="case2"
  )

  data.simu0 <- data.simu0 %>% mutate(
      strata1=sample(letters[1:3],n,replace=TRUE),
      strata2=sample(LETTERS[4:5],n,replace=TRUE)
  )

  out <- robincar_logrank(df=data.simu0,
                          treat_col="I1",
                          p_trt=0.5,
                          ref_arm=0,
                          response_col="t",
                          event_col="delta",
                          # covariate_cols=c("model_z1"),
                          car_scheme=c("simple"),
                          adj_method=c("CL")
  )

  surv_out <- coxph(Surv(t, delta) ~ I1, data = data.simu0)
  surv_out.score <-
    sqrt(surv_out$score) * ifelse(surv_out$coefficients > 0, 1, -1)

  expect_equal(out$result$statistic, unname(surv_out.score))

})
