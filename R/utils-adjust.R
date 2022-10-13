# Gets the vcovHC weights for the sample size and number of parameters
get.vcovHC <- function(vcovHC, n, p){
  if(vcovHC == "HC0"){
    wgt <- 1
  } else if(vcovHC == "HC1") {
    wgt <- n / (n - p)
  } else if(vcovHC == "HC3") {
    wgt <- (n / (n - p)) ** 2
  } else {
    stop("Unrecognized vcovHC method.")
  }
  return(wgt)
}

# Get design matrix for the specified adjustment variables
# using the data stored.
get.dmat <- function(data, adj_vars){
  if(is.null(adj_vars)){
    dmat <- NULL
  } else if(adj_vars == "x"){
    dmat <- data$covariate
  } else if(adj_vars == "z"){
    dmat <- data$strata
  } else if(adj_vars == "joint_z"){
    dmat <- as.matrix(data$joint_strata, ncol=1)
    colnames(dmat) <- "joint_strata"
  } else if(adj_vars == "joint_z_x"){
    joint_strata <- data$joint_strata
    dmat <- cbind(data$covariate, joint_strata)
  } else if(adj_vars == "formula"){
    dmat <- data$formula_vars
  } else {
    stop(paste("Unrecognized adjustment variable type ", adj_vars))
  }
  return(dmat)
}

#' @import stats
.center.dmat <- function(dmat){
  modmat <- stats::model.matrix(~ 0 + ., data=data.frame(dmat))
  modmat <- t(t(modmat) - colMeans(modmat))
  return(modmat)
}

