library(survival)
library(dplyr)

ting <- function(data){

  data.sort<-function(data.simu){
    n<-dim(data.simu)[1]
    data.simu.order<- data.simu[order(data.simu$t),]
    data.rev<- data.frame(data.simu.order,
                          Y=n:1,
                          Y1=cumsum(data.simu.order$I1[n:1])[n:1],
                          Y0=cumsum(data.simu.order$I0[n:1])[n:1])
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
    # T.seq<-data.rev$t
    # T.rep<-c(0,T.seq)[1:n]
    # T.diff<-T.seq-T.rep
    # same.ind<-which(T.diff==0)
    # n_col<-dim(data.rev)[2]
    # for(ind in same.ind){
    #   data.rev[ind,(n_col-2):n_col]<-data.rev[ind-1,(n_col-2):n_col]
    # }
    S.alt<-sum(data.rev$delta*(data.rev$I1*data.rev$Y0-data.rev$I0*data.rev$Y1)/data.rev$Y)
    var.alt<-sum(data.rev$delta*data.rev$Y0*data.rev$Y1*(data.rev$Y-data.rev$n.events)/data.rev$Y^2/(data.rev$Y-1),na.rm=TRUE)
    res<-S.alt/sqrt(var.alt)
    return(list(S.alt=S.alt, res=res, var.alt=var.alt))
  }

  return(wlogrank(data))

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

  tingres <- (2 * pnorm(abs(ting(data)$res), lower.tail = F))

  expect_equal(RC1res, lr_p)
  expect_equal(RC1res, tingres)

})

test_that("Logrank Ovarian with ties", {

  # Artificially create ties -----------------------------------
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

  tingres <- (2 * pnorm(abs(ting(data)$res), lower.tail = F))

  expect_equal(RC1res, lr_p)
  expect_equal(RC1res, tingres)

})
