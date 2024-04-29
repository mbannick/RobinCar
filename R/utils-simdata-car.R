#' Generate simple randomization treatment assignments
#'
#' @param n Number of observations
#' @param p_trt Proportion allotted to treatment
#' @export
#'
#' @importFrom stats rbinom
#'
#' @examples
#' car_sr(10, p_trt=0.4)
#'
#' @returns A vector of treatment assignments as 0's and 1's based on simple randomization.
car_sr <- function(n, p_trt){
  I <- stats::rbinom(n, 1, p_trt)
  return(I)
}

#' Generate permuted block treatment assignments
#'
#' @param z The car_strata design matrix, as a data frame with factor variables
#' @param trt_label Treatment label
#' @param trt_alc Treatment allocation vector
#' @param blocksize Permuted block blocksize
#'
#' @export
#' @importFrom dplyr summarise_all group_by_all cur_group_id mutate
#'
#' @examples
#' # Create car_strata variables
#' library(fastDummies)
#' library(dplyr)
#'
#' x <- runif(100)
#' z <- cut(x, breaks=c(0, 0.25, 0.5, 0.75, 1.0))
#' z <- dummy_cols(z) %>%
#'      mutate(across(where(is.numeric), as.factor))
#'
#' car_pb(z[, 2:5], c(0, 1, 2), trt_alc=c(1/4, 1/2, 1/4), blocksize=4L)
#'
#' @returns A vector of treatment assignments with labels from the `trt_label` argument, based on stratified permuted block randomization.
car_pb <- function(z, trt_label, trt_alc, blocksize=4L) {

  car_strata <- z

  if(!all(summarise_all(car_strata, is.factor) == TRUE)){
   stop("Please make sure every feature(column) of z is a factor!")
  }
  if(!is.integer(blocksize)){
    stop("blocksize should be an integer, e.g., blocksize = 4L")
  }

  ntrt_blk <- as.integer(blocksize*trt_alc/sum(trt_alc))
  A_blk <- rep(trt_label, ntrt_blk)

  if(length(A_blk) != blocksize){
    stop("car.blocksize*car.trt_alc/sum(car.trt_alc) should be integers!")
  }

  strata_cross <- car_strata %>%
    group_by_all %>%
    dplyr::mutate(strata_ind = dplyr::cur_group_id()) %>%
    ungroup()

  A <- rep(NA, nrow(car_strata))
  for (i in 1 : max(strata_cross$strata_ind)){
    car_strata.index.temp <- which(strata_cross$strata_ind == i)
    car_strata.size.temp <- length(car_strata.index.temp)
    A[car_strata.index.temp] <- as.vector(
      replicate(n = ceiling(car_strata.size.temp/blocksize),
                expr = sample(x = A_blk, size = blocksize))[1:car_strata.size.temp]
    )
  }

  return(A)
}

#' Generate Pocock-Simon minimization treatment assignments
#'
#' @param z The car_strata design matrix
#' @param treat A vector of length k (the number of treatment arms), which labels the treatment arms being compared.
#' @param ratio A vector of length k (the number of treatment arms), which indicates the allocation ratio, e.g., c(1,1,1) for equal allocation with three treatment arms.
#' @param imb_measure What measure of imbalance should be minimzed during randomization -- either "Range" or "SD"
#' @param p_bc The biased probability, i.e., the probability of assigning each patient to the arm that minimizes the imbalance. Default is 0.8
#' @return
#' \describe{
#' \item{res}{treatment assignment vector}
#' }
#'
#' @author Ting Ye Yanyao Yi
#' @import tidyverse
#' @import fastDummies
#' @importFrom stats sd
#' @export
#'
#' @examples
#'
#' # Create car_strata variables
#' library(fastDummies)
#' library(dplyr)
#'
#' x <- runif(100)
#' z <- cut(x, breaks=c(0, 0.25, 0.5, 0.75, 1.0))
#' z <- dummy_cols(z)
#' A <- car_ps(
#'   z=z[, 2:5],
#'   treat=c(0, 1, 2),
#'   ratio=c(1, 1, 1),
#'   imb_measure="Range"
#' )
#'
#' @returns A vector of treatment assignments with labels from the `treat` argument, based on Pocock-Simon's minimization.
car_ps <- function(z, treat, ratio, imb_measure, p_bc=0.8){

  if(!imb_measure %in% c("Range", "SD")) stop("Unrecognized imbalance measure.")
  z <- data.frame(z)

  miniRand <- function(j, result){

    ns <- ncol(z)
    ntrt <- length(treat)
    covwt <- rep(1 / ns, ns)

    if(j > 1){
      applydat <- z[1:(j - 1), , drop = FALSE]
      matchx <- apply(applydat, 1,
                      function(x, xrow) {
                        as.numeric(x == xrow)
                      }, z[j,])

      if(is.vector(matchx)){
        matchx <- matrix(matchx, nrow=1)
        n_matchx <- matrix(0, 1, ntrt)
        for(k in 1:ntrt){
          n_matchx[, k] <- sum(matchx[1, result[1:(j - 1)] == treat[k]])
        }
      } else{
        n_matchx <- matrix(0, ncol(z), ntrt)
        for (k in 1:ntrt) {
          kdat <- as.matrix(matchx[, result[1:(j - 1)] == treat[k]])
          n_matchx[, k] <- apply(kdat, 1, sum)
        }
      }

      imbalance <- rep(0, ntrt)
      for (i in 1:ntrt) {
        temp <- n_matchx
        temp[, i] <- temp[, i] + 1
        num_level <- temp %*% diag(1 / ratio)
        if(imb_measure == "Range"){
          range_level <- apply(num_level, 1, range)
          imb_margin <- range_level[2,] - range_level[1,]
        } else {
          sd_level <- apply(num_level, 1, stats::sd)
          imb_margin <- sd_level
        }
        imbalance[i] <- sum(covwt %*% imb_margin)
      }

      trt.highprob <- treat[imbalance == min(imbalance)]
      trt.lowprob <- treat[imbalance != min(imbalance)]

      res <- ifelse(
        length(trt.highprob) < ntrt,
        sample(
          c(trt.highprob, trt.lowprob),
          1,
          replace=TRUE,
          prob=c(rep(
            p_bc / length(trt.highprob), length(trt.highprob)
          ), rep((1 - p_bc) / length(trt.lowprob), length(trt.lowprob)
          ))
        ),
        sample(
          treat,
          1,
          replace=TRUE,
          prob=rep(1 / ntrt, ntrt)
        )
      )
    }
    return(res)
  }

  # Create result vector
  A <- rep(NA, nrow(z))
  A[1] <- sample(treat, 1, replace=TRUE,
                 prob=ratio / sum(ratio))

  # Loop through rows to create treatment assignments
  for(row in 2:nrow(z)){
    A[row] <- miniRand(j=row, result=A)
  }

  return(A)
}



