
diff <- function(est){
  base <- est[1]
  cont <- est[-1] - base
  return(cont)
}

odds_ratio <- function(est){
  base <- est[1] / (1-est[1])
  cont <- est[-1] / (1-est[-1]) / base
  return(cont)
}

log_odds_ratio <- function(est){
  base <- est[1] / (1-est[1])
  cont <- log(est[-1] / (1-est[-1]) / base)
  return(cont)
}

ratio <- function(est){
  base <- est[1]
  cont <- est[-1] / base
  return(cont)
}

log_ratio <- function(est){
  base <- est[1]
  cont <- log(est[-1] / base)
  return(cont)
}

jacobian <- function (settings, est) {
  UseMethod("jacobian", settings)
}

#' @exportS3Method
jacobian.DIFF <- function(settings, est){
  k <- length(est)
  jacob <- cbind(rep(-1, k-1), diag(1, k-1))
  return(jacob)
}

#' @exportS3Method
jacobian.RATIO <- function(settings, est){
  col1 <- -1 * est[-1] / est[1]^2
  mat1 <- diag(1/est[1], length(est)-1)
  jacob <- cbind(col1, mat1)
  return(jacob)
}

#' @exportS3Method
jacobian.LOGRATIO <- function(settings, est){
  col1 <- -1 / est[1]
  mat1 <- diag(1/est[-1], length(est)-1)
  jacob <- cbind(col1, mat1)
  return(jacob)
}

#' @exportS3Method
jacobian.ODDSRATIO <- function(settings, est){
  col1 <- -1 / est[1] ** 2 * est[-1] / (1-est[-1])
  mat1 <- diag((1-est[-1]) / est[-1] / (1-est[1]) ** 2, length(est)-1)
  jacob <- cbind(col1, mat1)
  return(jacob)
}

#' @exportS3Method
jacobian.LOGODDSRATIO <- function(settings, est){
  col1 <- -1 / (est[1] * (1 - est[1]))
  mat1 <- diag(1 / (est[-1] * (1 - est[-1])), length(est)-1)
  jacob <- cbind(col1, mat1)
  return(jacob)
}

#' @exportS3Method
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
  } else if(contrast_h == "log_ratio"){
    h <- log_ratio
    type <- "LOGRATIO"
    name <- function(tl) paste0("log(treat ", tl[-1], " / ", tl[1], ")")
  } else if(contrast_h == "odds_ratio"){
    h <- odds_ratio
    type <- "ODDSRATIO"
    name <- function(tl) paste0("odds(treat ", tl[-1], ") / odds(", tl[1], ")")
  } else if(contrast_h == "log_odds_ratio"){
    h <- log_odds_ratio
    type <- "LOGODDSRATIO"
    name <- function(tl) paste0("log(odds(treat ", tl[-1], ")) - log(odds(", tl[1], "))")
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

#' @importFrom emulator quad.tform
contrast <- function(settings, treat, est, varcov){

  # Get labels
  lab <- settings$name(treat)

  # If using an odds ratio contrast, make sure
  # give a warning if the estimates are not between 1 and 0.
  if(class(settings)[2] %in% c("ODDSRATIO", "LOGODDSRATIO")){
    if(any(est < 0) || any(est > 1)){
      warning("Estimates are not between 0 and 1: are you sure you want an odds ratio?")
    }
  }

  # If using a ratio or odds ratio, give a warning
  # that the performance is better with logs.
  if(class(settings)[2] %in% c("ODDSRATIO", "RATIO")){
    warning("Using a ratio or odds ratio is not recommended.
            Consider using log_ratio or log_odds_ratio for better
            performance of the variance estimator.")
  }

  # Get transformed estimate
  c_est <- settings$h(est)

  # Get Jacobian matrix
  jac <- jacobian(settings, est)

  # Get new variance covariance matrix based on old one & Jacobian
  vcv <- quad.tform(varcov, jac)

  # Format results
  result <- format_results(lab, c_est, vcv, label_name="contrast")

  return(list(result=result, varcov=vcv, settings=settings))
}

#' Estimate a treatment contrast
#'
#' Estimate a treatment contrast using the result of \link[RobinCar:robincar_linear]{RobinCar::robincar_linear()}, \link[RobinCar:robincar_glm]{RobinCar::robincar_glm()}, or \link[RobinCar:robincar_SL]{RobinCar::robincar_SL()} using
#' the delta method.
#'
#' @param result A LinModelResult or GLMModelResult
#' @param contrast_h An optional function to specify a desired contrast
#' @param contrast_dh An optional jacobian function for the contrast
#' @importFrom rlang .data
#' @export
#'
#' @returns A contrast object which has the following attributes:
#'
#'  \item{result}{A \link[dplyr:tibble]{dplyr::tibble()} with the label of the treatment contrast (e.g., 1 vs. 0), the estimate of the treatment contrast, estimated SE, and p-value based on a z-test with estimate and SE.}
#'  \item{varcov}{The variance-covariance matrix for the treatment contrast estimates.}
#'  \item{settings}{List of model settings used for the contrast.}
#'
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
