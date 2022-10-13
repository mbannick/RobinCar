test_that("logic for CL", {

  # 1. Z yes X yes
  # SR
  expect_warning(
    zyxy_sr <- ttelogic(
      x_exists=TRUE,
      z_exists=TRUE,
      car_scheme="simple",
      adj_method="CL"
    )
  )
  expect_true(zyxy_sr$adj_cov)
  expect_false(zyxy_sr$adj_strata)
  expect_false(zyxy_sr$car_strata)

  # CAR
  zyxy_car <- ttelogic(
    x_exists=TRUE,
    z_exists=TRUE,
    car_scheme=c("permuted-block"),
    adj_method=c("CL")
  )
  expect_true(zyxy_car$adj_cov)
  expect_true(zyxy_car$adj_strata)
  expect_false(zyxy_car$car_strata)


  # 2. Z yes X no
  # SR
  expect_warning(
    zyxn_sr <- ttelogic(
      x_exists=FALSE,
      z_exists=TRUE,
      car_scheme=c("simple"),
      adj_method=c("CL")
    )
  )
  expect_false(zyxn_sr$adj_cov)
  expect_false(zyxn_sr$adj_strata)
  expect_false(zyxn_sr$car_strata)

  # CAR
  zyxn_car <- ttelogic(
    x_exists=FALSE,
    z_exists=TRUE,
    car_scheme=c("biased-coin"),
    adj_method=c("CL")
  )
  expect_false(zyxn_car$adj_cov)
  expect_true(zyxn_car$adj_strata)
  expect_false(zyxn_car$car_strata)


  # 3. Z no X yes
  # SR
  znxy_sr = ttelogic(
    x_exists=TRUE,
    z_exists=FALSE,
    car_scheme=c("simple"),
    adj_method=c("CL")
  )
  expect_true(znxy_sr$adj_cov)
  expect_false(znxy_sr$adj_strata)
  expect_false(znxy_sr$car_strata)

  # CAR
  expect_error(
    znxy_car <- ttelogic(
      x_exists=TRUE,
      z_exists=FALSE,
      car_scheme=c("permuted-block"),
      adj_method=c("CL")
    )
  )

  # 4. Z no X no
  # SR
  znxn_sr <- ttelogic(
    x_exists=FALSE,
    z_exists=FALSE,
    car_scheme=c("simple"),
    adj_method=c("CL")
  )
  expect_false(znxn_sr$adj_cov)
  expect_false(znxn_sr$adj_strata)
  expect_false(znxn_sr$car_strata)

  # CAR
  expect_error(znxn_car <- ttelogic(
    x_exists=FALSE,
    z_exists=FALSE,
    car_scheme=c("permuted-block"),
    adj_method=c("CL")
  ))

})

test_that("logic for CSL", {

  # 1. Z yes X yes
  # SR
  expect_warning(
    zyxy_sr <- ttelogic(
      x_exists=TRUE,
      z_exists=TRUE,
      car_scheme="simple",
      adj_method="CSL"
    )
  )
  expect_true(zyxy_sr$adj_cov)
  expect_false(zyxy_sr$adj_strata)
  expect_false(zyxy_sr$car_strata)

  # CAR
  zyxy_car <- ttelogic(
    x_exists=TRUE,
    z_exists=TRUE,
    car_scheme=c("permuted-block"),
    adj_method=c("CSL")
  )
  expect_true(zyxy_car$adj_cov)
  expect_true(zyxy_car$adj_strata)
  expect_true(zyxy_car$car_strata)

  # 2. Z yes X no
  # SR
  expect_warning(
    zyxn_sr <- ttelogic(
      x_exists=FALSE,
      z_exists=TRUE,
      car_scheme=c("simple"),
      adj_method=c("CSL")
    )
  )
  expect_false(zyxn_sr$adj_cov)
  expect_false(zyxn_sr$adj_strata)
  expect_false(zyxn_sr$car_strata)

  # CAR
  zyxn_car <- ttelogic(
    x_exists=FALSE,
    z_exists=TRUE,
    car_scheme=c("biased-coin"),
    adj_method=c("CSL")
  )
  expect_false(zyxn_car$adj_cov)
  expect_true(zyxn_car$adj_strata)
  expect_true(zyxn_car$car_strata)


  # 3. Z no X yes
  # SR
  znxy_sr = ttelogic(
    x_exists=TRUE,
    z_exists=FALSE,
    car_scheme=c("simple"),
    adj_method=c("CSL")
  )
  expect_true(znxy_sr$adj_cov)
  expect_false(znxy_sr$adj_strata)
  expect_false(znxy_sr$car_strata)

  # CAR
  expect_error(
    znxy_car <- ttelogic(
      x_exists=TRUE,
      z_exists=FALSE,
      car_scheme=c("permuted-block"),
      adj_method=c("CSL")
    )
  )

  # 4. Z no X no
  # SR
  znxn_sr <- ttelogic(
    x_exists=FALSE,
    z_exists=FALSE,
    car_scheme=c("simple"),
    adj_method=c("CSL")
  )
  expect_false(znxn_sr$adj_cov)
  expect_false(znxn_sr$adj_strata)
  expect_false(znxn_sr$car_strata)

  # CAR
  expect_error(znxn_car <- ttelogic(
    x_exists=FALSE,
    z_exists=FALSE,
    car_scheme=c("permuted-block"),
    adj_method=c("CSL")
  ))

})

test_that("logic for CoxScore",{

  # 1. Z yes X yes
  # SR
  expect_warning(
    zyxy_sr <- ttelogic(
      x_exists=TRUE,
      z_exists=TRUE,
      car_scheme=c("simple"),
      adj_method=c("coxscore")
    )
  )
  expect_true(zyxy_sr$adj_cov)
  expect_false(zyxy_sr$adj_strata)
  expect_false(zyxy_sr$car_strata)

  # CAR
  zyxy_car <- ttelogic(
    x_exists=TRUE,
    z_exists=TRUE,
    car_scheme=c("permuted-block"),
    adj_method=c("coxscore")
  )
  expect_true(zyxy_car$adj_cov)
  expect_true(zyxy_car$adj_strata)
  expect_false(zyxy_car$car_strata)

  # 2. Z yes X no
  # SR
  expect_warning(
    zyxn_sr <- ttelogic(
      x_exists=FALSE,
      z_exists=TRUE,
      car_scheme=c("simple"),
      adj_method=c("coxscore")
    )
  )
  expect_false(zyxn_sr$adj_cov)
  expect_false(zyxn_sr$adj_strata)
  expect_false(zyxn_sr$car_strata)

  # CAR
  zyxn_car <- ttelogic(
    x_exists=FALSE,
    z_exists=TRUE,
    car_scheme = c("biased-coin"),
    adj_method = c("coxscore")
  )
  expect_false(zyxn_car$adj_cov)
  expect_true(zyxn_car$adj_strata)
  expect_false(zyxn_car$car_strata)

  # 3. Z no X yes
  # SR
  znxy_sr <- ttelogic(
    x_exists=TRUE,
    z_exists=FALSE,
    car_scheme=c("simple"),
    adj_method=c("coxscore")
  )
  expect_true(znxy_sr$adj_cov)
  expect_false(znxy_sr$adj_strata)
  expect_false(znxy_sr$car_strata)

  # CAR
  expect_error(
    znxy_car <- ttelogic(
      x_exists=TRUE,
      z_exists=FALSE,
      car_scheme=c("permuted-block"),
      adj_method=c("coxscore")
    )
  )

  # 4. Z no X no
  # SR
  znxn_sr <- ttelogic(
    x_exists=FALSE,
    z_exists=FALSE,
    car_scheme = c("simple"),
    adj_method = c("coxscore")
  )
  expect_false(znxn_sr$adj_cov)
  expect_false(znxn_sr$adj_strata)
  expect_false(znxn_sr$car_strata)

  # CAR
  expect_error(znxn_car <- ttelogic(
    x_exists=FALSE,
    z_exists=FALSE,
    car_scheme = c("permuted-block"),
    adj_method = c("coxscore")
  ))

})

