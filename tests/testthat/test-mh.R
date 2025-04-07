DATA <- RobinCar:::data_sim
DATA$y_bin <- ifelse(DATA$y > 2, 1, 0)
DATA_01 <- DATA[DATA$A!=2, ]
DATA_02 <- DATA[DATA$A!=1, ]
DATA_12 <- DATA[DATA$A!=0, ]

test_that("MH errors", {
  expect_error(robincar_mh(df = DATA,
                           treat_col = "A",
                           response_col = "y_bin",
                           strata_cols = c("z1", "z2"),
                           estimand = "MH",
                           ci_type = "mGR"))
  
  expect_error(robincar_mh(df = DATA_01,
                           treat_col = "A",
                           response_col = "y",
                           strata_cols = c("z1", "z2"),
                           estimand = "MH",
                           ci_type = "mGR"))
  expect_error(robincar_mh(df = DATA_12,
                           treat_col = "A",
                           response_col = "y_bin",
                           strata_cols = c("z1", "z2"),
                           estimand = "ATE",
                           ci_type = "GR"))
  
})

test_that("MH works", {
  fit.MH <- robincar_mh(df = DATA_01,
                       treat_col = "A",
                       response_col = "y_bin",
                       strata_cols = c("z1", "z2"),
                       estimand = "MH",
                       ci_type = "mGR")
  expect_equal(fit.MH$result$estimate, 0.2183309, tolerance = 1e-5)
  expect_equal(fit.MH$result$se, 0.03593319, tolerance = 1e-5)
  
  fit.ATE <- robincar_mh(df = DATA_01,
                         treat_col = "A",
                         response_col = "y_bin",
                         strata_cols = c("z1"),
                         estimand = "ATE",
                         ci_type = "mGR")
  expect_equal(fit.ATE$result$estimate, 0.2094867, tolerance = 1e-5)
  expect_equal(fit.ATE$result$se, 0.04296656, tolerance = 1e-5)
  
})