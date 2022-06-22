
diff <- function(est){
  base <- est[1]
  cont <- est[-1] - base
  return(cont)
}

ratio <- function(est){
  base <- est[1]
  cont <- est[-1] / base
  return(cont)
}

jacobian <- function (settings, ...) {
  UseMethod("jacobian", settings)
}

jacobian.DIFF <- function(settings, est){
  k <- length(est)
  jacob <- cbind(rep(-1, k-1), diag(1, k-1))
  return(jacob)
}

jacobian.RATIO <- function(settings, est){
  col1 <- -1 * est[-1] / est[1]^2
  mat1 <- diag(1/est[1], length(est)-1)
  jacob <- cbind(col1, mat1)
  return(jacob)
}
    
jacobian.CUSTOM <- function(settings, est){
  if(!is.null(settings$dh)){
    jacob <- settings$dh(est)
  } else {
    jacob <- numDeriv::jacobian(
      func=settings$h, 
      x=est
    )
  }
  return(jacob)
}

# Create contrast settings based on a contrast function
contrast.settings <- function(k, contrast_h, contrast_dh=NULL){
  
  if(is.function(contrast_h)){
    h <- contrast_h
    type <- "CUSTOM"
    d <- length(h(rep(1, k))) # size of output
    name <- function(tl) paste0("contrast ", 1:d)
  } else if(contrast_h == "diff"){
    h <- diff
    type <- "DIFF"
    name <- function(tl) paste0("treat ", tl[-1], " - " , tl[1])
  } else if(contrast_h == "ratio"){
    h <- ratio
    type <- "RATIO"
    name <- function(tl) paste0("treat ", tl[-1], " / ", tl[1])
  } else {
    stop("Unrecognized contrast function name.")
  }
  
  obj <- structure(
    list(
      h=h,
      dh=contrast_dh,
      name=name
    ),
    class=c("CONTRAST", type)
  )
  return(obj)
}

# Run a contrast using settings specified and estimates from a model
contrast <- function(settings, treat, est, varcov){
  
  # Get labels
  lab <- settings$name(treat)
  
  # Get transformed estimate
  c_est <- settings$h(est)
  
  # Get Jacobian matrix
  jac <- jacobian(settings, est)
  
  # Get new variance covariance matrix based on old one & Jacobian
  vcv <- quad.tform(varcov, jac)
  
  # Format results
  result <- format.results(lab, c_est, vcv, label_name="contrast")
  
  return(list(result=result, varcov=vcv, settings=settings))
}

# Create a contrast using a result object (that has both result and varcov)
# 
#' @import dplyr
#' @import numDeriv
#' @import emulator
#' @import tidyverse
#' @importFrom rlang .data
#' @export
robincar_contrast <- function(result, contrast_h, contrast_dh=NULL){

  # Get input dimension
  k <- nrow(result$result)
  
  # Create contrast settings
  settings <- contrast.settings(k, contrast_h, contrast_dh)
  
  # Run the contrast function
  c_result <- contrast(
    settings, 
    treat=result$result$treat, 
    est=result$result$estimate, 
    varcov=result$varcov
  )
  return(
    structure(
      class="ContrastResult",
      list(result=c_result$result, 
           varcov=c_result$varcov,
           settings=settings)
    )
  )
}
