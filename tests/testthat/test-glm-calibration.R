DATA <- RobinCar:::data_sim
DATA$A <- as.factor(DATA$A)
DATA$y <- ifelse(DATA$y > 2, 1, 0)

test_that("test calibration", {
  this <- suppressWarnings(robincar_glm(
    df=DATA,
    response_col="y",
    treat_col="A",
    car_strata_cols=c("z1"),
    car_scheme="biased-coin",
    g_family=binomial(link="logit"),
    g_accuracy=7,
    formula="y ~ A + x3"))

  that <- Robin_g(
    robin.data=DATA,
    car.scheme="stratified_biased_coin",
    car.z=c("z1"),
    robin.x=c("x3"),
    robin.x_to_include_z=FALSE,
    robin.g_family=binomial(link="logit"),
    robin.g_accuracy=7,
    robin.formula=NULL,
    robin.adj_method="homogeneous"
  )

  for(joint in c(FALSE, TRUE)){
    for(add_x in c("x3", NULL)){

      calib_this <- suppressWarnings(robincar_calibrate(
        result=this,
        joint=joint,
        add_x=add_x
      ))
      calib_that <- Robin_calibrate(
        robin.g_object=that$robin.g_object,
        robin.add_x=add_x,
        robin.add_x_to_include_z=joint,
        robin.vcovHC="HC0"
      )
      calib_this_est <- calib_this$result$estimate
      calib_this_se <- calib_this$result$se

      calib_this_vc <- calib_this$varcov
      dimnames(calib_this_vc) <- NULL
      calib_that_vc <- calib_that$varcov
      dimnames(calib_that_vc) <- NULL

      names(calib_this_est) <- NULL
      names(calib_this_se) <- NULL
      expect_equal(calib_this_est, calib_that$estimation$estimate)
      # NEW: variance will no longer be the same now that we have done
      # the finite sample correction for the variance estimation
      # when there are lots of 0's in the response vector
      # expect_equal(calib_this_se, calib_that$estimation$se)
      # expect_equal(calib_this_vc, calib_that_vc)
    }
  }
})
