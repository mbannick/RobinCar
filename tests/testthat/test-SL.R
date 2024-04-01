library(SuperLearner)

n <- 1000
set.seed(10)
DATA2 <- data.frame(A=rbinom(n, size=1, prob=0.5),
                    y=rbinom(n, size=1, prob=0.2),
                    x1=rnorm(n),
                    x2=rnorm(n),
                    x3=as.factor(rbinom(n, size=1, prob=0.5)),
                    z1=rbinom(n, size=1, prob=0.5),
                    z2=rbinom(n, size=1, prob=0.5))
DATA2[, "y"] <- NA
As <- DATA2$A == 1
DATA2[DATA2$A == 1, "y"] <- rbinom(sum(As), size=1, prob=exp(DATA2[As,]$x1)/(1+exp(DATA2[As,]$x1)))
DATA2[DATA2$A == 0, "y"] <- rbinom(n-sum(As), size=1, prob=exp(1 + 5*DATA2[!As,]$x1 + DATA2[!As,]$x2)/(1+exp(1 + 5*DATA2[!As,]$x1 + DATA2[!As,]$x2)))
DATA2$A <- as.factor(DATA2$A)

test_that("simple super learner with random forest.", {

  sl.mod <- robincar_SL(
    df=DATA2,
    response_col="y",
    treat_col="A",
    car_strata_cols=c("z1"),
    covariate_cols=c("x1"),
    SL_libraries=c("SL.ranger"),
    car_scheme="permuted-block",
    covariate_to_include_strata=TRUE
  )
})

test_that("simple super learner with median adjustment, random forest.", {

  sl.mod <- robincar_SL_median(
    seed=1,
    n_times=3,
    df=DATA2,
    response_col="y",
    treat_col="A",
    car_strata_cols=c("z1"),
    covariate_cols=c("x1"),
    SL_libraries=c("SL.ranger"),
    car_scheme="permuted-block",
    covariate_to_include_strata=TRUE
  )

})

test_that("super learner with custom learner", {

  skip("This works, but not in the testthat environment.")

  create_rf <- create.Learner("SL.ranger",
                              tune = list(num.trees=c(25, 50),
                                          max.depth=c(0.1, 0.5)))

  sl.mod <- robincar_SL(
    df=DATA2,
    response_col="y",
    treat_col="A",
    car_strata_cols=c("z1"),
    covariate_cols=c("x1"),
    SL_learners=list(create_rf),
    car_scheme="permuted-block",
    covariate_to_include_strata=TRUE
  )

})
