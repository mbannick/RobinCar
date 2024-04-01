# Compare to legacy code from Yanyao Yi

test_that("GLM legacy", {

  n <- 10000
  set.seed(10)
  DATA2 <- data.frame(A=rbinom(n, size=1, prob=0.5),
                      y=rbinom(n, size=1, prob=0.2),
                      x1=rnorm(n),
                      x2=rnorm(n),
                      x3=as.factor(rbinom(n, size=1, prob=0.5)),
                      z1=rbinom(n, size=1, prob=0.5),
                      z2=rbinom(n, size=1, prob=0.5))
  DATA2$A <- as.factor(DATA2$A)
  DATA2$x1 <- DATA2$x1 - mean(DATA2$x1)

  for(cov in list(c("x1"), NULL)){
    for(ctis in list(TRUE, FALSE)){
      for(zs in list(c("z1", "z2"), c("z1"), NULL)){
        for(meth in c("homogeneous", "heterogeneous")){
          for(scheme in c("SR",
                          "PB",
                          # "minimization",
                          # TODO: Fix these tests to be in line with new Z logic for pocock simon!
                          "stratified_biased_coin")){

            if(scheme == "SR"){
              scheme0 <- "simple"
            } else if(scheme == "PB"){
              scheme0 <- "permuted-block"
            } else if(scheme == "minimization"){
              scheme0 <- "pocock-simon"
            } else if(scheme == "stratified_biased_coin"){
              scheme0 <- "biased-coin"
            }

            if(scheme != "SR" & is.null(zs)) next
            # We can't compare the behavior in this case
            # because old code automatically includes Z when heterogeneous + pocock simon
            if(scheme == "minimization" & !ctis & meth == "heterogeneous") next

            runthis <- function(){
              return(robincar_glm(
                df=DATA2,
                response_col="y",
                treat_col="A",
                car_strata_cols=zs,
                covariate_cols=cov,
                covariate_to_include_strata=ctis,
                car_scheme=scheme0,
                g_family=binomial(link="logit"),
                g_accuracy=7,
                adj_method=meth
              ))
            }

            runthat <- function(){
              return(Robin_g(
                robin.data=DATA2,
                car.scheme=scheme,
                car.z=zs,
                robin.x=cov,
                robin.x_to_include_z=ctis,
                robin.g_family=binomial(link="logit"),
                robin.g_accuracy=7,
                robin.formula=NULL,
                robin.vcovHC="HC0",
                robin.adj_method=meth
              ))
            }
            if(scheme == "minimization" & meth == "homogeneous"){
              expect_warning(runthis())
              expect_error(runthat())
            } else {
              this <- runthis()
              that <- runthat()
              this_est <- this$result$estimate
              this_se <- this$result$se
              names(this_est) <- NULL
              names(this_se) <- NULL
              expect_equal(this_est, that$estimation$estimate)
              # NEW: variance will no longer be the same now that we have done
              # the finite sample correction for the variance estimation
              # when there are lots of 0's in the response vector
              # expect_equal(this_se, that$estimation$se)
            }
          }
        }
      }
    }
  }
})

