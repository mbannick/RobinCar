DATA <- RobinCar:::data_sim
DATA$y_bin <- ifelse(DATA$y > 2, 1, 0)
DATA_01 <- DATA[DATA$A!=2, ]
DATA_02 <- DATA[DATA$A!=1, ]
DATA_12 <- DATA[DATA$A!=0, ]

test_that("MH errors", {
  expect_error(robincar_mh(df = DATA,
                           treat_col = "A",
                           response_col = "y",
                           strata_cols = c("z1", "z2"),
                           estimand = "MH",
                           ci_type = "mGR"),
               "The treatment variable must be binary.")
  
  expect_error(robincar_mh(df = DATA_01,
                           treat_col = "A",
                           response_col = "y",
                           strata_cols = c("z1", "z2"),
                           estimand = "MH",
                           ci_type = "mGR"))
  expect_error(robincar_mh(df = DATA_12,
                           treat_col = "A",
                           response_col = "y",
                           strata_cols = c("z1", "z2"),
                           estimand = "ATE",
                           ci_type = "GR"))
  
})

test_that("MH works", {
  expect_no_message(robincar_mh(df = DATA_01,
                                treat_col = "A",
                                response_col = "y_bin",
                                strata_cols = c("z1", "z2"),
                                estimand = "MH",
                                ci_type = "mGR"))
  expect_no_message(robincar_mh(df = DATA_02,
                                treat_col = "A",
                                response_col = "y_bin",
                                strata_cols = c("z1"),
                                estimand = "MH",
                                ci_type = "Sato"))
  
  expect_no_message(robincar_mh(df = DATA_12,
                                treat_col = "A",
                                response_col = "y_bin",
                                strata_cols = c("z1"),
                                estimand = "ATE",
                                ci_type = "mGR"))
  expect_no_message(robincar_mh(df = DATA_12,
                                treat_col = "A",
                                response_col = "y",
                                strata_cols = c("z2"),
                                estimand = "ATE",
                                ci_type = "mGR"))
  
})