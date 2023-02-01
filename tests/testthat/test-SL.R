n <- 1000
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

glm.mod <- robincar_glm(
  df=DATA2,
  response_col="y",
  treat_col="A",
  strata_cols=c("z1"),
  covariate_cols=c("x1"),
  car_scheme="permuted-block",
  g_family=binomial(link="logit"),
  g_accuracy=7,
  adj_method="heterogeneous",
  covariate_to_include_strata=TRUE,
  vcovHC="HC3")

sl.mod <- robincar_SL(
  df=DATA2,
  response_col="y",
  treat_col="A",
  strata_cols=c("z1"),
  covariate_cols=c("x1"),
  SL_libraries=c("SL.ranger"),
  car_scheme="permuted-block",
  covariate_to_include_strata=TRUE,
  vcovHC="HC3"
)
