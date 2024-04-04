###### Author: Ting Ye ######
###### Last modified: Jan 16, 2022 ######
# Use null variance #
data.sort<-function(data.simu){
  n<-dim(data.simu)[1]
  data.simu.order<- data.simu[order(data.simu$t),]
  data.rev<- data.frame(data.simu.order, Y=n:1, Y1=cumsum(data.simu.order$I1[n:1])[n:1], Y0=cumsum(data.simu.order$I0[n:1])[n:1])
  return(data.rev)
}

wlogrank<-function(data.simu){
  n<-dim(data.simu)[1]
  data.rev<-data.sort(data.simu)
  T.seq<-data.rev$t
  T.rep<-c(0,T.seq)[1:n]
  T.diff<-T.seq-T.rep
  same.ind<-which(T.diff==0)
  n_col<-dim(data.rev)[2]
  for(ind in same.ind){
    data.rev[ind,(n_col-2):n_col]<-data.rev[ind-1,(n_col-2):n_col]
  }
  S.alt<-sum(data.rev$delta*(data.rev$I1*data.rev$Y0-data.rev$I0*data.rev$Y1)/data.rev$Y)
  var.alt<-sum(data.rev$delta*data.rev$Y0*data.rev$Y1/data.rev$Y^2)
  res<-S.alt/sqrt(var.alt)
  return(list(S.alt=S.alt, res=res, var.alt=var.alt))
}

ind_to_factor<-function(data.simu){
  z<-data.simu[,grepl("car_strata",names(data.simu))]
  n_strata<-dim(z)[2]
  n<-dim(z)[1]
  fz<-numeric(n)
  for(i in 1:n_strata){
    fz<-fz+z[,i]*i
  }
  return(as.factor(fz))
}

wlogrank_var_calibrated<-function(data.simu,randomization,p_trt){ # try new version 19 May
  n<-dim(data.simu)[1]
  S.alt<-wlogrank(data.simu)$S.alt
  data.rev<-data.sort(data.simu)
  T.seq<-data.rev$t
  T.rep<-c(0,T.seq)[1:n]
  T.diff<-T.seq-T.rep
  same.ind<-which(T.diff==0)
  n_col<-dim(data.rev)[2]
  for(ind in same.ind){
    data.rev[ind,(n_col-2):n_col]<-data.rev[ind-1,(n_col-2):n_col]
  }

  strata_z<-data.rev[,grepl("car_strata",names(data.rev))]
  mean_at_risk<-data.rev$delta/data.rev$Y
  cumsum_at_risk<-cumsum(mean_at_risk)
  mean_at_risk_2<-data.rev$delta*data.rev$Y1/(data.rev$Y)^2
  cumsum_at_risk_2<-cumsum(mean_at_risk_2)
  mu_t<-data.rev$Y1/data.rev$Y
  O_i<-data.rev$delta*(data.rev$I1-data.rev$Y1/data.rev$Y)-data.rev$I1*cumsum_at_risk+cumsum_at_risk_2
  O_i1<-data.rev$I1*O_i
  O_i0<- -data.rev$I0*O_i
  strata_z<-data.rev[,grepl("car_strata",names(data.rev))]
  n_strata<-dim(strata_z)[2]
  z<-strata_z
  cond_var_1<-numeric(n_strata)
  cond_var_0<-numeric(n_strata)
  prob_z<-numeric(n_strata)
  cond_exp_1<-numeric(n_strata)
  cond_exp_0<-numeric(n_strata)
  for(i in 1:n_strata){
    cond_var_1[i]<-var(O_i1[z[,i]==1 & data.rev$I1==1])
    cond_var_0[i]<-var(O_i0[z[,i]==1 & data.rev$I0==1])
    prob_z[i]<-mean(z[,i])
    cond_exp_1[i]<-mean(O_i1[z[,i]==1 & data.rev$I1==1])
    cond_exp_0[i]<-mean(O_i0[z[,i]==1 & data.rev$I0==1])
  }
  na.ind<-is.na(cond_var_1+cond_var_0)
  if(length(which(na.ind))!=0) warning("conditional variance in log-rank equals 0")
  cond_var_1<-cond_var_1[!na.ind]
  cond_var_0<-cond_var_0[!na.ind]
  prob_z<-prob_z[!na.ind]
  cond_exp_1<-cond_exp_1[!na.ind]
  cond_exp_0<-cond_exp_0[!na.ind]
  sai_1<-(p_trt*(cond_var_1%*%prob_z)+(1-p_trt)*(cond_var_0%*%prob_z))
  sai_2<-0 # for cabc and permuted_block
  if(randomization=="SR"){
    sai_2<-p_trt*(1-p_trt)*sum((cond_exp_0+cond_exp_1)^2 * prob_z)
  }
  if(randomization =="urn"){
    if(p_trt!=1/2)  stop("For now, urn design is only designed for pi=1/2")
    sai_2<-1/12*sum((cond_exp_0+cond_exp_1)^2 * prob_z)
  }
  var_cal<-n*(sai_1+sai_2)
  T_cal<-S.alt/sqrt(var_cal)
  return(list(S.alt=S.alt, res=T_cal,var_cal=var_cal))
}

# # proposed T_RS in Ye and Shao, 2020, jrssb.
# score_robust<-function(data.simu,randomization,p_trt){
#   data.rev<-data.sort(data.simu)
#   x<-data.rev[,grepl("model",names(data.rev))] #only those has model
#   x<-cbind(x,data.rev[,grepl("car_strata",names(data.rev))][,-1])
#   n<-dim(data.rev)[1]
#   data_sub<-data.frame(t=data.rev$t,delta=data.rev$delta,x)
#   fit<-survival::coxph(survival::Surv(t,delta)~.,data=data_sub)
#   S_0_seq<-exp(fit$coefficients %*% t(x))
#   S_1_seq<-exp(fit$coefficients %*% t(x))*data.rev$I1
#   S_0<-cumsum(S_0_seq[n:1])[n:1]/n
#   S_1<-cumsum(S_1_seq[n:1])[n:1]/n
#   U<-sum(data.rev$delta*(data.rev$I1-S_1/S_0))
#   #mu_t<-data.rev$Y1/data.rev$Y # two choices, both valid
#   mu_t<-S_1/S_0
#   mean_at_risk<-data.rev$delta/(n*S_0)
#   cumsum_at_risk<-cumsum(mean_at_risk)
#   mean_at_risk_2<-data.rev$delta/(n*S_0)*mu_t
#   cumsum_at_risk_2<-cumsum(mean_at_risk_2)
#   O_i<-data.rev$delta*(data.rev$I1-mu_t)-S_1_seq*cumsum_at_risk+S_0_seq*cumsum_at_risk_2
#
#   O_i1<-data.rev$I1*O_i
#   O_i0<- -data.rev$I0*O_i
#
#   strata_z<-data.rev[,grepl("car_strata",names(data.rev))]
#   n_strata<-dim(strata_z)[2]
#   z<-strata_z
#   cond_var_1<-numeric(n_strata)
#   cond_var_0<-numeric(n_strata)
#   prob_z<-numeric(n_strata)
#   cond_exp_1<-numeric(n_strata)
#   cond_exp_0<-numeric(n_strata)
#   for(i in 1:n_strata){
#     cond_var_1[i]<-var(O_i1[z[,i]==1 & data.rev$I1==1])
#     cond_var_0[i]<-var(O_i0[z[,i]==1 & data.rev$I0==1])
#     prob_z[i]<-mean(z[,i])
#     cond_exp_1[i]<-mean(O_i1[z[,i]==1 & data.rev$I1==1])
#     cond_exp_0[i]<-mean(O_i0[z[,i]==1 & data.rev$I0==1])
#   }
#   na.ind<-is.na(cond_var_1+cond_var_0)
#   if(length(which(na.ind))!=0) warning("conditional variance in log-rank equals 0")
#   cond_var_1<-cond_var_1[!na.ind]
#   cond_var_0<-cond_var_0[!na.ind]
#   prob_z<-prob_z[!na.ind]
#   cond_exp_1<-cond_exp_1[!na.ind]
#   cond_exp_0<-cond_exp_0[!na.ind]
#   sai_1<-(p_trt*(cond_var_1%*%prob_z)+(1-p_trt)*(cond_var_0%*%prob_z))
#   sai_2<-0 # for cabc and permuted_block
#   if(randomization=="SR"){
#     sai_2<-p_trt*(1-p_trt)*sum((cond_exp_0+cond_exp_1)^2 * prob_z)
#   }
#   if(randomization =="urn"){
#     if(p_trt!=1/2)  stop("For now, urn design is only designed for pi=1/2")
#     sai_2<-1/12*sum((cond_exp_0+cond_exp_1)^2 * prob_z)
#   }
#   var_cal<-n*(sai_1+sai_2)
#   T_SR<-U/sqrt(var_cal)
#   return(T_SR)
# }


# proposed T_RS in Ye and Shao, 2020, jrssb.
score_robust<-function(data.simu,randomization,p_trt){
  data.rev<-data.sort(data.simu)
  x<-data.rev[,grepl("model",names(data.rev))] #only those has model
  n<-dim(data.rev)[1]
  data_sub<-data.frame(t=data.rev$t,delta=data.rev$delta,x)
  fit<-survival::coxph(survival::Surv(t,delta)~.,data=data_sub)
  S_0_seq<-exp(fit$coefficients %*% t(x))
  S_1_seq<-exp(fit$coefficients %*% t(x))*data.rev$I1
  S_0<-cumsum(S_0_seq[n:1])[n:1]/n
  S_1<-cumsum(S_1_seq[n:1])[n:1]/n
  U<-sum(data.rev$delta*(data.rev$I1-S_1/S_0))
  #mu_t<-data.rev$Y1/data.rev$Y # two choices, both valid
  mu_t<-S_1/S_0
  mean_at_risk<-data.rev$delta/(n*S_0)
  cumsum_at_risk<-cumsum(mean_at_risk)
  mean_at_risk_2<-data.rev$delta/(n*S_0)*mu_t
  cumsum_at_risk_2<-cumsum(mean_at_risk_2)
  O_i<-data.rev$delta*(data.rev$I1-mu_t)-S_1_seq*cumsum_at_risk+S_0_seq*cumsum_at_risk_2

  O_i1<-data.rev$I1*O_i
  O_i0<- -data.rev$I0*O_i

  strata_z<-data.rev[,grepl("car_strata",names(data.rev))]
  n_strata<-dim(strata_z)[2]
  z<-strata_z
  cond_var_1<-numeric(n_strata)
  cond_var_0<-numeric(n_strata)
  prob_z<-numeric(n_strata)
  cond_exp_1<-numeric(n_strata)
  cond_exp_0<-numeric(n_strata)
  for(i in 1:n_strata){
    cond_var_1[i]<-var(O_i1[z[,i]==1 & data.rev$I1==1])
    cond_var_0[i]<-var(O_i0[z[,i]==1 & data.rev$I0==1])
    prob_z[i]<-mean(z[,i])
    cond_exp_1[i]<-mean(O_i1[z[,i]==1 & data.rev$I1==1])
    cond_exp_0[i]<-mean(O_i0[z[,i]==1 & data.rev$I0==1])
  }
  na.ind<-is.na(cond_var_1+cond_var_0)
  if(length(which(na.ind))!=0) warning("conditional variance in log-rank equals 0")
  cond_var_1<-cond_var_1[!na.ind]
  cond_var_0<-cond_var_0[!na.ind]
  prob_z<-prob_z[!na.ind]
  cond_exp_1<-cond_exp_1[!na.ind]
  cond_exp_0<-cond_exp_0[!na.ind]
  sai_1<-(p_trt*(cond_var_1%*%prob_z)+(1-p_trt)*(cond_var_0%*%prob_z))
  sai_2<-0 # for cabc and permuted_block
  if(randomization=="SR"){
    # if(p_trt!=1/2)  stop("For now, SR is only designed for pi=1/2")
    # sai_2<-1/4*sum((cond_exp_0+cond_exp_1)^2 * prob_z)
    sai_2<-p_trt*(1-p_trt)*sum((cond_exp_0+cond_exp_1)^2 * prob_z)
  }
  # if(randomization =="urn"){
  #   if(p_trt!=1/2)  stop("For now, urn design is only designed for pi=1/2")
  #   sai_2<-1/12*sum((cond_exp_0+cond_exp_1)^2 * prob_z)
  # }
  var_cal<-n*(sai_1+sai_2)
  T_SR<-U/sqrt(var_cal)
  return(list(T_SR=T_SR,
              se=sqrt(var_cal),
              U=U))
}


# score test T_S in Ye and Shao, 2020, jrssb.
score<-function(data.simu){
  data.rev<-data.sort(data.simu)
  x<-data.rev[,grepl("model",names(data.rev))] #only those has model
  x<-cbind(x,data.rev[,grepl("car_strata",names(data.rev))][,-1])
  n<-dim(data.rev)[1]
  data_sub<-data.frame(t=data.rev$t,delta=data.rev$delta,x)
  fit<-survival::coxph(survival::Surv(t,delta)~.,data=data_sub)
  S_0_seq<-exp(fit$coefficients %*% t(x))
  S_1_seq<-exp(fit$coefficients %*% t(x))*data.rev$I1
  S_0<-cumsum(S_0_seq[n:1])[n:1]/n
  S_1<-cumsum(S_1_seq[n:1])[n:1]/n
  U<-sum(data.rev$delta*(data.rev$I1-S_1/S_0))
  #mu_t<-data.rev$Y1/data.rev$Y # two choices, both valid
  mu_t<-S_1/S_0
  mean_at_risk<-data.rev$delta/(n*S_0)
  cumsum_at_risk<-cumsum(mean_at_risk)
  mean_at_risk_2<-data.rev$delta/(n*S_0)*mu_t
  cumsum_at_risk_2<-cumsum(mean_at_risk_2)
  O_i<-data.rev$delta*(data.rev$I1-mu_t)-S_1_seq*cumsum_at_risk+S_0_seq*cumsum_at_risk_2
  var_cal<-mean(O_i^2)*n
  T_WR<-U/sqrt(var_cal)
  return(T_WR)
}

# numerator of T_S
score_numerator<-function(data.simu){
  data.rev<-data.sort(data.simu)
  x<-data.rev[,grepl("model",names(data.rev))] #only those has model
  x<-cbind(x,data.rev[,grepl("car_strata",names(data.rev))][,-1])
  n<-dim(data.rev)[1]
  data_sub<-data.frame(t=data.rev$t,delta=data.rev$delta,x)
  fit<-survival::coxph(survival::Surv(t,delta)~.,data=data_sub)
  S_0_seq<-exp(fit$coefficients %*% t(x))
  S_1_seq<-exp(fit$coefficients %*% t(x))*data.rev$I1
  S_0<-cumsum(S_0_seq[n:1])[n:1]/n
  S_1<-cumsum(S_1_seq[n:1])[n:1]/n
  U<-sum(data.rev$delta*(data.rev$I1-S_1/S_0))
  return(U)
}

#########################
#   Added Sep 5, 2021   #
#########################

# Using log-rank score as the derived outcome (Tangen-Koch, 1999)
logrank_TK_score<-function(data.simu,p_trt,randomization){
  n<-dim(data.simu)[1]
  data.simu$fz<-ind_to_factor(data.simu)
  data.rev<-data.sort(data.simu)
  T.seq<-data.rev$t
  T.rep<-c(0,T.seq)[1:n]
  T.diff<-T.seq-T.rep
  same.ind<-which(T.diff==0)
  n_col<-dim(data.rev)[2]
  for(ind in same.ind){
    data.rev[ind,(n_col-2):n_col]<-data.rev[ind-1,(n_col-2):n_col]
  }
  n<-dim(data.rev)[1]
  data.rev$R<-with(data.rev,delta-cumsum(data.rev$delta/data.rev$Y))

  x.model<-as.matrix(data.rev[,grepl("model",names(data.rev)),drop=FALSE]) #only those has model
  x.mat<-model.matrix(~factor(fz)+x.model,data=data.rev)[,-1]
  x.centered<-sweep(x.mat, 2, colMeans(x.mat))
  fit1<-lm(R~factor(fz)+x.model[data.rev$I1==1,],data=data.rev[data.rev$I1==1,])
  fit0<-lm(R~factor(fz)+x.model[data.rev$I1==0,],data=data.rev[data.rev$I1==0,])
  b1<-fit1$coefficients[-1]
  b0<-fit0$coefficients[-1]
  # calculating beta0 and beta1 for variance estimation
  mu_t<-data.rev$Y1/data.rev$Y
  data.rev$O.hat<-data.rev$delta*(data.rev$I1-data.rev$Y1/data.rev$Y)-
    data.rev$I1*cumsum(data.rev$delta/data.rev$Y)+cumsum(data.rev$delta*data.rev$Y1/(data.rev$Y)^2)
  data.rev$O.hat[data.rev$I0==1]<- (-1)*data.rev$O.hat[data.rev$I0==1]
  fit1.Ohat<-lm(O.hat~factor(fz)+x.model[data.rev$I1==1,],data=data.rev[data.rev$I1==1,])
  fit0.Ohat<-lm(O.hat~factor(fz)+x.model[data.rev$I1==0,],data=data.rev[data.rev$I1==0,])
  beta1.Ohat<-fit1.Ohat$coefficients[-1]
  beta0.Ohat<-fit0.Ohat$coefficients[-1]

  ind.na<-which(is.na(b1))
  if(length(ind.na)==0){
    U_TK<-sum(data.rev$I1*(data.rev$R-t((x.centered)%*% b1 )) -
                data.rev$I0*(data.rev$R-t((x.centered)%*% b0)))/n/2
  }else{
    warning("Removing model variables that are linearly dependent with the stratification variables.")
    b1<-b1[-ind.na]
    b0<-b0[-ind.na]
    x.mat<-x.mat[,-ind.na]
    x.centered<-x.centered[,-ind.na]
    beta1.Ohat<-beta1.Ohat[-ind.na]
    beta0.Ohat<-beta0.Ohat[-ind.na]
    U_TK<-sum(data.rev$I1*(data.rev$R-t((x.centered)%*% b1 )) -
                data.rev$I0*(data.rev$R-t((x.centered)%*% b0)))/n/2
  }
  var_CL.null<-(sum(data.rev$delta*data.rev$Y0*data.rev$Y1/data.rev$Y^2)/n-
                  p_trt*(1-p_trt)*(beta1.Ohat+beta0.Ohat) %*% var(x.mat) %*% (beta1.Ohat+beta0.Ohat) )/n
  var_TK.null<-var_CL.null+p_trt*(1-p_trt)*(beta1.Ohat-b1/2+beta0.Ohat-b0/2)%*%var(x.mat)%*%(beta1.Ohat-b1/2+beta0.Ohat-b0/2)/n
  if(randomization%in% c("permuted_block","SR","CABC")){
    se<-sqrt(var_TK.null)
  }else{
    se<-NA
  }
  return(list(U_TK=U_TK,se=se))
}


# covariate adjusted logrank proposed in Ye, Yi, Shao (2022)
#'@title Ting Ye's original code
#'@export
covariate_adjusted_logrank<-function(data.simu,p_trt){
  n<-dim(data.simu)[1]
  data.simu$fz<-ind_to_factor(data.simu)
  data.rev<-data.sort(data.simu)
  T.seq<-data.rev$t
  T.rep<-c(0,T.seq)[1:n]
  T.diff<-T.seq-T.rep
  same.ind<-which(T.diff==0)
  n_col<-dim(data.rev)[2]
  for(ind in same.ind){
    data.rev[ind,(n_col-2):n_col]<-data.rev[ind-1,(n_col-2):n_col]
  }
  n<-dim(data.rev)[1]

  x.model<-as.matrix(data.rev[,grepl("model",names(data.rev)),drop=FALSE]) #only those has model
  x.mat<-model.matrix(~factor(fz)+x.model,data=data.rev)[,-1]
  x.centered<-sweep(x.mat, 2, colMeans(x.mat))
  # calculating beta0 and beta1 for variance estimation
  mu_t<-data.rev$Y1/data.rev$Y
  data.rev$O.hat<-data.rev$delta*(data.rev$I1-data.rev$Y1/data.rev$Y)-
    data.rev$I1*cumsum(data.rev$delta/data.rev$Y)+cumsum(data.rev$delta*data.rev$Y1/(data.rev$Y)^2)
  data.rev$O.hat[data.rev$I0==1]<- (-1)*data.rev$O.hat[data.rev$I0==1]
  fit1.Ohat<-lm(O.hat~factor(fz)+x.model[data.rev$I1==1,],data=data.rev[data.rev$I1==1,])
  fit0.Ohat<-lm(O.hat~factor(fz)+x.model[data.rev$I1==0,],data=data.rev[data.rev$I1==0,])
  beta1.Ohat<-fit1.Ohat$coefficients[-1]
  beta0.Ohat<-fit0.Ohat$coefficients[-1]

  ind.na<-which(is.na(beta1.Ohat))
  if(length(ind.na)==0){
    U_CL<-sum(data.rev$I1*(data.rev$O.hat-t((x.centered)%*% beta1.Ohat )) -
                data.rev$I0*(data.rev$O.hat-t((x.centered)%*% beta0.Ohat)))/n
  }else{
    warning("Removing model variables that are linearly dependent with the stratification variables.")
    x.mat<-x.mat[,-ind.na]
    x.centered<-x.centered[,-ind.na]
    beta1.Ohat<-beta1.Ohat[-ind.na]
    beta0.Ohat<-beta0.Ohat[-ind.na]
    U_CL<-sum(data.rev$I1*(data.rev$O.hat-t((x.centered)%*% beta1.Ohat )) -
                data.rev$I0*(data.rev$O.hat-t((x.centered)%*% beta0.Ohat)))/n
  }

  score_logrank_anhecova_var<-function(data.rev,p_trt,fit1.Ohat,fit0.Ohat,beta1.Ohat,beta0.Ohat){
    my.var<-(p_trt*var(fit1.Ohat$residuals)+(1-p_trt)*var(fit0.Ohat$residuals)+
               (p_trt*beta1.Ohat-(1-p_trt)*beta0.Ohat) %*% var(x.mat) %*% (p_trt*beta1.Ohat-(1-p_trt)*beta0.Ohat))/n
    my.var.null<-(sum(data.rev$delta*data.rev$Y0*data.rev$Y1/data.rev$Y^2)/n-
                    p_trt*(1-p_trt)*(beta1.Ohat+beta0.Ohat) %*% var(x.mat) %*% (beta1.Ohat+beta0.Ohat) )/n
    return(list(my.var=my.var,my.var.null=my.var.null))
  }
  tmp<-score_logrank_anhecova_var(data.rev,p_trt,fit1.Ohat,fit0.Ohat,beta1.Ohat,beta0.Ohat)
  se<-sqrt(tmp$my.var)
  se.null<-sqrt(tmp$my.var.null)
  return(list(U_CL=U_CL,se=se.null,se.orig=se))
}


# covariate adjusted stratified logrank proposed in Ye, Yi, Shao (2022)
covariate_adjusted_stratified_logrank<-function(data.simu,p_trt){
  n<-dim(data.simu)[1]
  data.simu$fz<-ind_to_factor(data.simu)
  fz.n<-unique(data.simu$fz)
  data.rev<-data.sort(data.simu[data.simu$fz==fz.n[1],])
  T.seq<-data.rev$t
  T.rep<-c(0,T.seq)[1:length(T.seq)]
  T.diff<-T.seq-T.rep
  same.ind<-which(T.diff==0)
  n_col<-dim(data.rev)[2]
  for(ind in same.ind){
    data.rev[ind,(n_col-2):n_col]<-data.rev[ind-1,(n_col-2):n_col]
  }
  mu_t<-data.rev$Y1/data.rev$Y
  data.rev$O.hat<-data.rev$delta*(data.rev$I1-data.rev$Y1/data.rev$Y)-
    data.rev$I1*cumsum(data.rev$delta/data.rev$Y)+cumsum(data.rev$delta*data.rev$Y1/(data.rev$Y)^2)
  data.rev$O.hat[data.rev$I0==1]<- (-1)*data.rev$O.hat[data.rev$I0==1]

  # calculate stratum-specific O_zij
  for(z in fz.n[-1]){
    data.tmp<-data.sort(data.simu[data.simu$fz==z,])
    T.seq<-data.tmp$t
    T.rep<-c(0,T.seq)[1:length(T.seq)]
    T.diff<-T.seq-T.rep
    same.ind<-which(T.diff==0)
    n_col<-dim(data.tmp)[2]
    for(ind in same.ind){
      data.tmp[ind,(n_col-2):n_col]<-data.tmp[ind-1,(n_col-2):n_col]
    }
    mu_t<-data.tmp$Y1/data.tmp$Y
    data.tmp$O.hat<-data.tmp$delta*(data.tmp$I1-data.tmp$Y1/data.tmp$Y)-
      data.tmp$I1*cumsum(data.tmp$delta/data.tmp$Y)+cumsum(data.tmp$delta*data.tmp$Y1/(data.tmp$Y)^2)
    data.tmp$O.hat[data.tmp$I0==1]<- (-1)*data.tmp$O.hat[data.tmp$I0==1]
    data.rev<-rbind(data.rev,data.tmp)
  }
  n<-dim(data.rev)[1]

  x.model<-as.matrix(data.rev[,grepl("model",names(data.rev)),drop=FALSE]) #only those has model

  fit.anhc1<-lm(O.hat~factor(fz)+x.model[data.rev$I1==1,],data=data.rev[data.rev$I1==1,])
  fit.anhc0<-lm(O.hat~factor(fz)+x.model[data.rev$I1==0,],data=data.rev[data.rev$I1==0,])
  which.model<-which(grepl("x.model",names(fit.anhc1$coefficients)))
  which.model.notNA<-which(!is.na(fit.anhc1$coefficients[which.model]))
  if(length(which.model.notNA)==0){
    stop("The model variables are linearly dependent on Z. The covariate-adjusted stratified log-rank is the same as the stratified log-rank.")
  }else{
    fit.anhc1.coefx<-fit.anhc1$coefficients[which.model[which.model.notNA]]
    fit.anhc0.coefx<-fit.anhc0$coefficients[which.model[which.model.notNA]]

    x.coefx<-model.matrix(~factor(fz)+as.matrix(x.model),data=data.rev)[,which.model[which.model.notNA],drop=FALSE]
    for(z in fz.n){
      x.coefx[data.rev$fz==z,]<-sweep(x.coefx[data.rev$fz==z,,drop=FALSE], 2, colMeans(x.coefx[data.rev$fz==z,,drop=FALSE]))
    }
    x.coefx.centered.byZ<-x.coefx
    U_CSL<- sum(data.rev$I1*data.rev$O.hat- data.rev$I0*data.rev$O.hat)-
      sum(data.rev$I1*as.vector(x.coefx.centered.byZ %*% fit.anhc1.coefx) - data.rev$I0*as.vector(x.coefx.centered.byZ %*% fit.anhc0.coefx))
    U_CSL<-U_CSL/n
  }
  ### calculate variance ###
  var_CSL<-0
  for(ind in 1:length(fz.n)){
    z<-fz.n[ind]
    Bz1<-var(data.rev$O.hat[data.rev$fz==z & data.rev$I1==1]-x.coefx[data.rev$fz==z & data.rev$I1==1,,drop=FALSE]%*%fit.anhc1.coefx)
    Bz0<-var(data.rev$O.hat[data.rev$fz==z & data.rev$I1==0]-x.coefx[data.rev$fz==z & data.rev$I1==0,,drop=FALSE]%*%fit.anhc0.coefx)
    if(is.na(Bz1) | is.na(Bz0)){
      warning("conditional variance equals 0")
      Bz1<-ifelse(is.na(Bz1),0,Bz1)
      Bz0<-ifelse(is.na(Bz0),0,Bz0)
    }
    SigmaXz<-var(x.coefx[data.rev$fz==z,,drop=FALSE])
    pz<-length(which(data.rev$fz==z))/n
    var_CSL<-var_CSL+pz*(p_trt*Bz1+(1-p_trt)*Bz0)+
      pz*(p_trt*fit.anhc1.coefx-(1-p_trt)*fit.anhc0.coefx) %*% SigmaXz %*% (p_trt*fit.anhc1.coefx-(1-p_trt)*fit.anhc0.coefx)
  }


  tmp_SL_var<-0
  tmp_CSL_var<-0
  for(z in fz.n){
    tmp_SL_var<-tmp_SL_var+wlogrank(data.rev[data.rev$fz==z,])$var.alt
    pz<-length(which(data.rev$fz==z))/n
    SigmaXz<-var(x.coefx[data.rev$fz==z,,drop=FALSE])
    tmp_CSL_var<-tmp_CSL_var+wlogrank(data.rev[data.rev$fz==z,])$var.alt/n -
      p_trt*(1-p_trt)*(fit.anhc1.coefx+fit.anhc0.coefx) %*% SigmaXz %*%(fit.anhc1.coefx+fit.anhc0.coefx)*pz
  }
  se<-sqrt(var_CSL/n)
  se.null<-sqrt(tmp_CSL_var/n)
  return(list(U_CSL=U_CSL,se=se.null,se.orig=se))
}

# Covariate adjusted logrank proposed in Ye, Yi, Shao (2022)
#' @title Ting Ye's original code
#' @export
covariate_adjusted_logrank <- function(data.simu, p_trt, Z=TRUE) {
  n <- dim(data.simu)[1]
  if(Z) data.simu$fz <- ind_to_factor(data.simu)
  data.rev <- data.sort(data.simu)
  T.seq <- data.rev$t
  T.rep <- c(0, T.seq)[1:n]
  T.diff <- T.seq - T.rep
  same.ind <- which(T.diff == 0)
  n_col <- dim(data.rev)[2]
  for (ind in same.ind) {
    data.rev[ind, (n_col - 2):n_col] <- data.rev[ind - 1, (n_col - 2):n_col]
  }
  n <- dim(data.rev)[1]

  x.model <-
    as.matrix(data.rev[, grepl("model", names(data.rev)), drop = FALSE]) #only those has model
  if(Z){
    x.mat <- model.matrix( ~ factor(fz) + x.model, data = data.rev)[, -1]
  } else {
    x.mat <- model.matrix( ~ x.model, data = data.rev)[, -1]
  }

  # Need this if there's only one covariate and no car_strata
  x.mat <- as.matrix(x.mat)
  x.centered <- sweep(x.mat, 2, colMeans(x.mat))
  # calculating beta0 and beta1 for variance estimation
  mu_t <- data.rev$Y1 / data.rev$Y
  data.rev$O.hat <-
    data.rev$delta * (data.rev$I1 - data.rev$Y1 / data.rev$Y) -
    data.rev$I1 * cumsum(data.rev$delta / data.rev$Y) + cumsum(data.rev$delta *
                                                                 data.rev$Y1 / (data.rev$Y) ^ 2)
  data.rev$O.hat[data.rev$I0 == 1] <-
    (-1) * data.rev$O.hat[data.rev$I0 == 1]
  if(Z){
    fit1.Ohat <-
      lm(O.hat ~ factor(fz) + x.model[data.rev$I1 == 1, ], data = data.rev[data.rev$I1 ==
                                                                             1, ])
    fit0.Ohat <-
      lm(O.hat ~ factor(fz) + x.model[data.rev$I1 == 0, ], data = data.rev[data.rev$I1 ==
                                                                             0, ])
  } else {
    fit1.Ohat <-
      lm(O.hat ~ x.model[data.rev$I1 == 1, ], data = data.rev[data.rev$I1 ==
                                                                1, ])
    fit0.Ohat <-
      lm(O.hat ~ x.model[data.rev$I1 == 0, ], data = data.rev[data.rev$I1 ==
                                                                0, ])
  }

  beta1.Ohat <- fit1.Ohat$coefficients[-1]
  beta0.Ohat <- fit0.Ohat$coefficients[-1]

  ind.na <- which(is.na(beta1.Ohat))
  if (length(ind.na) == 0) {
    U_CL <-
      sum(data.rev$I1 * (data.rev$O.hat - t((x.centered) %*% beta1.Ohat)) -
            data.rev$I0 * (data.rev$O.hat - t((x.centered) %*% beta0.Ohat))) /
      n
  } else{
    warning(
      "Removing model variables that are linearly dependent with the stratification variables."
    )
    x.mat <- x.mat[, -ind.na]
    x.centered <- x.centered[, -ind.na]
    beta1.Ohat <- beta1.Ohat[-ind.na]
    beta0.Ohat <- beta0.Ohat[-ind.na]
    U_CL <-
      sum(data.rev$I1 * (data.rev$O.hat - t((x.centered) %*% beta1.Ohat)) -
            data.rev$I0 * (data.rev$O.hat - t((x.centered) %*% beta0.Ohat))) /
      n
  }

  score_logrank_anhecova_var <-
    function(data.rev,
             p_trt,
             fit1.Ohat,
             fit0.Ohat,
             beta1.Ohat,
             beta0.Ohat) {
      my.var <-
        (
          p_trt * var(fit1.Ohat$residuals) + (1 - p_trt) * var(fit0.Ohat$residuals) +
            (p_trt * beta1.Ohat - (1 - p_trt) * beta0.Ohat) %*% var(x.mat) %*% (p_trt *
                                                                                  beta1.Ohat - (1 - p_trt) * beta0.Ohat)
        ) / n
      my.var.null <-
        (
          sum(data.rev$delta * data.rev$Y0 * data.rev$Y1 / data.rev$Y ^ 2) / n -
            p_trt * (1 - p_trt) * (beta1.Ohat + beta0.Ohat) %*% var(x.mat) %*% (beta1.Ohat +
                                                                                  beta0.Ohat)
        ) / n
      return(list(my.var = my.var, my.var.null = my.var.null))
    }
  tmp <-
    score_logrank_anhecova_var(data.rev, p_trt, fit1.Ohat, fit0.Ohat, beta1.Ohat, beta0.Ohat)
  se <- sqrt(tmp$my.var)
  se.null <- sqrt(tmp$my.var.null)
  return(list(
    U_CL = U_CL,
    se = se.null,
    se.orig = se
  ))
}

