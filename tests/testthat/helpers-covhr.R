###### Author: Ting Ye ######
###### Last modified: Jan 1, 2023 ######
# Use null variance #
# Add estimation #

#########################
#   Added Jan 1, 2023   #
#########################

# covariate adjusted hazard ratio estimation proposed in Ye, Yi, Shao (2022)
covariate_adjusted_logrank_estimation<-function(data.simu,p_trt){
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
  score<-function(theta){
    S.alt<-sum(data.rev$delta*(data.rev$I1-exp(theta)*data.rev$Y1/(exp(theta)*data.rev$Y1+data.rev$Y0)))/n
    return(S.alt)
  }
  theta_L<-uniroot(score,interval=c(-10,10))$root # get the unadjusted estimator

  x.model<-as.matrix(data.rev[,grepl("model",names(data.rev)),drop=FALSE]) #only those has model
  if(length(levels(data.rev$fz))==1){
    x.mat<-model.matrix(~x.model,data=data.rev)[,-1,drop=FALSE]
  }else{
    x.mat<-model.matrix(~factor(fz)+x.model,data=data.rev)[,-1,drop=FALSE]
  }
  x.centered<-sweep(x.mat, 2, colMeans(x.mat))
  # calculating beta0 and beta1
  data.rev$O.hat1<-data.rev$delta*data.rev$Y0/(exp(theta_L)*data.rev$Y1+data.rev$Y0) -
    exp(theta_L)*cumsum(data.rev$delta*data.rev$Y0/(exp(theta_L)*data.rev$Y1+data.rev$Y0)^2)
  data.rev$O.hat0<-exp(theta_L)*data.rev$delta*data.rev$Y1/(exp(theta_L)*data.rev$Y1+data.rev$Y0) -
    exp(theta_L)*cumsum(data.rev$delta*data.rev$Y1/(exp(theta_L)*data.rev$Y1+data.rev$Y0)^2)

  fit1.Ohat<-lm(data.rev$O.hat1[data.rev$I1==1]~x.mat[data.rev$I1==1,])
  fit0.Ohat<-lm(data.rev$O.hat0[data.rev$I1==0]~x.mat[data.rev$I1==0,])
  beta1.Ohat<-fit1.Ohat$coefficients[-1]
  beta0.Ohat<-fit0.Ohat$coefficients[-1]

  ind.na<-which(is.na(beta1.Ohat))
  if(length(ind.na)!=0){
    warning("Removing model variables that are linearly dependent with the stratification variables.")
    x.mat<-x.mat[,-ind.na,drop=FALSE]
    x.centered<-x.centered[,-ind.na,drop=FALSE]
    beta1.Ohat<-beta1.Ohat[-ind.na]
    beta0.Ohat<-beta0.Ohat[-ind.na]
  }
  score_CL<-function(theta){
    res<-score(theta)-sum(data.rev$I1*(t((x.centered)%*% beta1.Ohat )) -
                            data.rev$I0*(t((x.centered)%*% beta0.Ohat)))/n
    return(res)
  }
  theta_CL<-uniroot(score_CL,interval=c(-10,10))$root # get the adjusted estimator

  # Modification for making consistent with current code for testing (theta_L in place of theta_CL)
  sigma2_L<-exp(theta_L)*sum(data.rev$delta*data.rev$Y0*data.rev$Y1/(exp(theta_L)*data.rev$Y1+data.rev$Y0)^2)/n
  sigma2_CL<-sigma2_L-p_trt*(1-p_trt)*(beta1.Ohat+beta0.Ohat) %*% var(x.mat) %*% (beta1.Ohat+beta0.Ohat)
  var_est<-sigma2_CL/sigma2_L^2/n
  return(list(theta_L=theta_L,se.theta_L=sqrt(1/sigma2_L/n),
              theta_CL=theta_CL,se.theta_CL=sqrt(var_est)))
}

covariate_adjusted_stratified_logrank_estimation<-function(data.simu,p_trt){
  n<-dim(data.simu)[1]
  data.simu$fz<-ind_to_factor(data.simu)
  fz.n<-unique(data.simu$fz)
  # obtain U_SL
  score_SL<-function(theta){
    S.alt<-0
    for(z in fz.n){
      data.tmp<-data.sort(data.simu[data.simu$fz==z,])
      T.seq<-data.tmp$t
      T.rep<-c(0,T.seq)[1:length(T.seq)]
      T.diff<-T.seq-T.rep
      same.ind<-which(T.diff==0)
      n_col<-dim(data.tmp)[2]
      for(ind in same.ind){
        data.tmp[ind,(n_col-2):n_col]<-data.tmp[ind-1,(n_col-2):n_col]
      }
      S.alt<-S.alt+sum(data.tmp$delta*(data.tmp$I1-exp(theta)*data.tmp$Y1/(exp(theta)*data.tmp$Y1+data.tmp$Y0)))/n
    }
    return(S.alt)
  }
  theta_SL<-uniroot(score_SL,interval=c(-10,10))$root # get the unadjusted stratified estimator

  data.rev<-data.sort(data.simu[data.simu$fz==fz.n[1],])
  T.seq<-data.rev$t
  T.rep<-c(0,T.seq)[1:length(T.seq)]
  T.diff<-T.seq-T.rep
  same.ind<-which(T.diff==0)
  n_col<-dim(data.rev)[2]
  for(ind in same.ind){
    data.rev[ind,(n_col-2):n_col]<-data.rev[ind-1,(n_col-2):n_col]
  }
  data.rev$O.hat1<-data.rev$delta*data.rev$Y0/(exp(theta_SL)*data.rev$Y1+data.rev$Y0) -
    exp(theta_SL)*cumsum(data.rev$delta*data.rev$Y0/(exp(theta_SL)*data.rev$Y1+data.rev$Y0)^2)
  data.rev$O.hat0<-exp(theta_SL)*data.rev$delta*data.rev$Y1/(exp(theta_SL)*data.rev$Y1+data.rev$Y0) -
    exp(theta_SL)*cumsum(data.rev$delta*data.rev$Y1/(exp(theta_SL)*data.rev$Y1+data.rev$Y0)^2)

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
    data.tmp$O.hat1<- data.tmp$O.hat1<-data.tmp$delta*data.tmp$Y0/(exp(theta_SL)*data.tmp$Y1+data.tmp$Y0) -
      exp(theta_SL)*cumsum(data.tmp$delta*data.tmp$Y0/(exp(theta_SL)*data.tmp$Y1+data.tmp$Y0)^2)
    data.tmp$O.hat0<-exp(theta_SL)*data.tmp$delta*data.tmp$Y1/(exp(theta_SL)*data.tmp$Y1+data.tmp$Y0) -
      exp(theta_SL)*cumsum(data.tmp$delta*data.tmp$Y1/(exp(theta_SL)*data.tmp$Y1+data.tmp$Y0)^2)
    data.rev<-rbind(data.rev,data.tmp)
  }
  n<-dim(data.rev)[1]

  x.model<-as.matrix(data.rev[,grepl("model",names(data.rev)),drop=FALSE]) #only those has model
  fit.anhc1<-lm(O.hat1~factor(fz)+x.model[data.rev$I1==1,],data=data.rev[data.rev$I1==1,])
  fit.anhc0<-lm(O.hat0~factor(fz)+x.model[data.rev$I1==0,],data=data.rev[data.rev$I1==0,])
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
    score_CSL<-function(theta){
      res<-score_SL(theta)-
        sum(data.rev$I1*as.vector(x.coefx.centered.byZ %*% fit.anhc1.coefx) - data.rev$I0*as.vector(x.coefx.centered.byZ %*% fit.anhc0.coefx))/n
      return(res)
    }
    theta_CSL<-uniroot(score_CSL,interval=c(-10,10))$root # get the adjusted estimator
  }
  ### calculate variance ###
  # Modification for making consistent with current code for testing (theta_SL in place of theta_CSL)
  sigma2_SL<-exp(theta_SL)*sum(data.rev$delta*data.rev$Y0*data.rev$Y1/(exp(theta_SL)*data.rev$Y1+data.rev$Y0)^2)/n
  sigma2_CSL<-sigma2_SL
  for(z in fz.n){
    pz<-length(which(data.rev$fz==z))/n
    SigmaXz<-var(x.coefx[data.rev$fz==z,,drop=FALSE])
    sigma2_CSL<-sigma2_CSL-  p_trt*(1-p_trt)*(fit.anhc1.coefx+fit.anhc0.coefx) %*% SigmaXz %*%(fit.anhc1.coefx+fit.anhc0.coefx)*pz
  }
  var_est<-sigma2_CSL/sigma2_SL^2/n
  return(list(theta_SL=theta_SL,se.theta_SL=sqrt(1/sigma2_SL/n),
              theta_CSL=theta_CSL,se.theta_CSL=sqrt(var_est)))
}



