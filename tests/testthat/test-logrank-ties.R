library(survival)
library(dplyr)

ting <- function(data){

  data.sort<-function(data.simu){
    n<-dim(data.simu)[1]
    data.simu.order<- data.simu[order(data.simu$t),]
    data.rev<- data.frame(data.simu.order,
                          Y=n:1,
                          Y1=cumsum(data.simu.order$I1[n:1])[n:1],
                          Y0=cumsum(data.simu.order$I0[n:1])[n:1]
    )

    data.rev <- data.rev %>% # handle ties, added by Ting
      group_by(t) %>%
      mutate(
        Y = first(Y),
        Y1 = first(Y1),
        Y0 = first(Y0),
        n.events= sum(delta)
      ) %>%
      ungroup()
    return(data.rev)
  }

  wlogrank<-function(data.simu){
    n<-dim(data.simu)[1]
    data.rev<-data.sort(data.simu)
    S.alt<-sum(data.rev$delta*(data.rev$I1*data.rev$Y0-data.rev$I0*data.rev$Y1)/data.rev$Y)
    var.alt<-sum(data.rev$delta*data.rev$Y0*data.rev$Y1*(data.rev$Y-data.rev$n.events)/data.rev$Y^2/(data.rev$Y-1),na.rm=TRUE)
    res<-S.alt/sqrt(var.alt)
    return(list(S.alt=S.alt, res=res, var.alt=var.alt))
  }

  # covariate adjusted logrank proposed in Ye, Yi, Shao (2022)
  covariate_adjusted_logrank <- function(data.simu){
    p_trt<-mean(data.simu$I1)
    n<-dim(data.simu)[1]
    data.simu$fz<-factor(1) # for illustration purpose, just one stratum
    data.rev<-data.sort(data.simu)

    data.rev <- data.rev %>% mutate(cum.delta.Y0.over.Ysq=cumsum(Y0/Y^2*delta), # added by Ting to handle ties
                                    cum.delta.Y1.over.Ysq=cumsum(Y1/Y^2*delta))
    data.rev <- data.rev %>%
      group_by(t) %>%
      mutate(
        cum.delta.Y0.over.Ysq=last(cum.delta.Y0.over.Ysq),
        cum.delta.Y1.over.Ysq=last(cum.delta.Y1.over.Ysq)
      ) %>%
      ungroup()

    mu_t<-data.rev$Y1/data.rev$Y
    data.rev$O.hat1<-data.rev$delta*(data.rev$Y0/data.rev$Y)-
      data.rev$cum.delta.Y0.over.Ysq
    data.rev$O.hat0<-data.rev$delta*(data.rev$Y1/data.rev$Y)-
      data.rev$cum.delta.Y1.over.Ysq
    data.rev$O.hat<-data.rev$I1*data.rev$O.hat1+data.rev$I0*data.rev$O.hat0

    sum(data.rev$I1*data.rev$O.hat -
          data.rev$I0*data.rev$O.hat)/n

    x.centered<-x.mat<-matrix(0,nrow=n,ncol=1)
    beta1.Ohat<-0
    beta0.Ohat<-0

    ind.na<-which(is.na(beta1.Ohat))
    if(length(ind.na)==0){
      U_CL<-sum(data.rev$I1*(data.rev$O.hat-t((x.centered)%*% beta1.Ohat )) -
                  data.rev$I0*(data.rev$O.hat-t((x.centered)%*% beta0.Ohat)))/n
    }else{
      warning("Removing model variables that are linearly dependent with the stratification variables.")
      x.mat<-x.mat[,-ind.na,drop=FALSE]
      x.centered<-x.centered[,-ind.na,drop=FALSE]
      beta1.Ohat<-beta1.Ohat[-ind.na]
      beta0.Ohat<-beta0.Ohat[-ind.na]
      U_CL<-sum(data.rev$I1*(data.rev$O.hat-t((x.centered)%*% beta1.Ohat )) -
                  data.rev$I0*(data.rev$O.hat-t((x.centered)%*% beta0.Ohat)))/n
    }

    score_logrank_anhecova_var<-function(data.rev,p_trt,fit1.Ohat,fit0.Ohat,beta1.Ohat,beta0.Ohat){
      my.var.null<-(wlogrank(data.simu)$var.alt/n -
                      # added by Ting to handle ties, uses log-rank variance estimator - variance reduction term
                      p_trt*(1-p_trt)*(beta1.Ohat+beta0.Ohat) %*% var(x.mat) %*% (beta1.Ohat+beta0.Ohat) )/n
      return(list(my.var.null=my.var.null))
    }
    tmp<-score_logrank_anhecova_var(data.rev,p_trt,fit1.Ohat,fit0.Ohat,beta1.Ohat,beta0.Ohat)
    se.null<-sqrt(tmp$my.var.null)
    return(list(U_CL=U_CL, se=se.null, res=U_CL/se.null))
  }

  return(covariate_adjusted_logrank(data))

}

DATA <- ovarian %>%
  rename(tte = futime, obs = fustat) %>%
  arrange(tte)

test_that("Logrank ovarian", {

  # (1) Survdiff function
  lr <- survdiff(Surv(tte, obs) ~ rx, data = DATA)
  lr_p <- lr$pvalue
  lr_num <- lr$obs[1] - lr$exp[1]
  lr_var <- lr$var[1, 1] / nrow(DATA)

  RC1 <- robincar_logrank(
    adj_method = "CL",
    df = DATA,
    treat_col = "rx",
    p_trt = 0.5,
    ref_arm = 1,
    response_col = "tte",
    event_col = "obs"
  )
  RC1res <- (2 * pnorm(abs(RC1$result$statistic), lower.tail = F))

  data <- DATA %>% select(t = tte)
  data$delta <- 1*(DATA$obs==T)
  data$I1 <- DATA$rx - 1
  data$I0 <- 1 - data$I1

  tingres <- c(2 * pnorm(abs(ting(data)$res), lower.tail = F))

  expect_equal(RC1res, lr_p)
  expect_equal(RC1res, tingres)

})

test_that("Logrank Ovarian with ties", {

  # Artificially create ties -----------------------------------
  # ------------------------------------------------------------
  # ------------------------------------------------------------
  DATA2 <- DATA

  DATA2$tte[3:5] <- DATA2$tte[3]
  DATA2$tte[10:12] <- DATA2$tte[10]
  DATA2$tte[19:21] <- DATA2$tte[19]

  # (1) Survdiff function
  lr <- survdiff(Surv(tte, obs) ~ rx, data = DATA2)
  lr_p <- lr$pvalue
  lr_num <- lr$obs[1] - lr$exp[1]
  lr_var <- lr$var[1, 1] / nrow(DATA)

  RC1 <- robincar_logrank(
    adj_method = "CL",
    df = DATA2,
    treat_col = "rx",
    p_trt = 0.5,
    ref_arm = 1,
    response_col = "tte",
    event_col = "obs"
  )
  RC1res <- (2 * pnorm(abs(RC1$result$statistic), lower.tail = F))

  data <- DATA2 %>% select(t = tte)
  data$delta <- 1*(DATA2$obs==T)
  data$I1 <- DATA2$rx - 1
  data$I0 <- 1 - data$I1

  tingres <- c(2 * pnorm(abs(ting(data)$res), lower.tail = F))

  expect_equal(RC1res, lr_p)
  expect_equal(RC1res, tingres)

})
