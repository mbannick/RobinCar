#' Generate simple randomization treatment assignments
#' Currently only works for two treatment groups
#'
#' @param n Number of observations
#' @param p_trt Proportion allotted to treatment
#' @export
#'
#' @example
#' SR(10, p_trt=0.4)
SR <- function(n, p_trt){
  I <- rbinom(n, 1, p_trt)
  return(I)
}

#' Generate biased coin treatment groups
#' Efron's biased coin randomization (Efron, 1971)
#'
#' @param z The strata design matrix
#' @param BCD_p The proportion for the biased coin (only two treatment groups)
#'
#' @export
#' @examples
#'
#' # Create strata variables
#' x <- runif(100)
#' z <- cut(x, breaks=c(0, 0.25, 0.5, 0.75, 1.0))
#' z <- model.matrix(~0 + factor(z))
#'
#' CABC(z, BCD_p=2/3)
CABC <- function(z, BCD_p=2/3){

  BC <- function(n_within_strata){

    stopifnot(BCD_p > 1 / 2)
    I <- numeric(n_within_strata)
    I[1] <- rbinom(1, 1, 1 / 2)
    D <- 2 * I[1] - 1

    for (i in 2:n_within_strata) {
      if (D > 0)
        I[i] <- rbinom(1, 1, (1 - BCD_p))
      else if (D < 0)
        I[i] <- rbinom(1, 1, BCD_p)
      else if (D == 0)
        I[i] <- rbinom(1, 1, 1 / 2)
      D <- D + 2 * I[i] - 1
    }
    return(I)
  }

  n_strata <- dim(z)[2]
  n <- dim(z)[1]
  n_each_strata <- apply(z, 2, sum)
  I <- numeric(n)
  for (s in 1:n_strata) {
    if (n_each_strata[s] == 0) {
      next
    }
    else if (n_each_strata[s] == 1) {
      I[z[, s] == 1] <- rbinom(1, 1, 0.5)
    }
    else{
      I[z[, s] == 1] <- BC(n_each_strata[s])
    }
  }
  return(I)
}

#' Generate permuted block treatment assignments
#'
#' @param z The strata design matrix, as a data frame with factor variables
#' @param trt_label Treatment label
#' @param trt_alc Treatment allocation vector
#' @param blocksize Permuted block blocksize
#'
#' @export
#' @importFrom dplyr summarise_all
#'
#' @examples
#' # Create strata variables
#' library(fastDummies)
#' library(dplyr)
#'
#' x <- runif(100)
#' z <- cut(x, breaks=c(0, 0.25, 0.5, 0.75, 1.0))
#' z <- dummy_cols(z) %>%
#'      mutate(across(where(is.numeric), as.factor))
#'
#' PB(z[, 2:5], c(0, 1, 2), trt_alc=c(1/4, 1/2, 1/4), blocksize=4L)
PB <- function(z, trt_label, trt_alc, blocksize=4L) {

  strata <- z

  if(!all(summarise_all(strata, is.factor) == TRUE)){
   stop("Please make sure every feature(column) of z is a factor!")
  }
  if(!is.integer(blocksize)){
    stop("blocksize should be an integer, e.g., blocksize = 4L")
  }

  ntrt_blk <- as.integer(blocksize*trt_alc/sum(trt_alc))
  A_blk <- rep(trt_label, ntrt_blk)

  if(length(A_blk) != blocksize){
    stop(print("car.blocksize*car.trt_alc/sum(car.trt_alc) should be integers!"))
  }

  strata_cross <- strata %>%
    group_by_all %>%
    mutate(strata_ind = cur_group_id()) %>%
    ungroup()

  A <- rep(NA, nrow(strata))
  for (i in 1 : max(strata_cross$strata_ind)){
    strata.index.temp <- which(strata_cross$strata_ind == i)
    strata.size.temp <- length(strata.index.temp)
    A[strata.index.temp] <- as.vector(
      replicate(n = ceiling(strata.size.temp/blocksize),
                expr = sample(x = A_blk, size = blocksize))[1:strata.size.temp]
    )
  }

  return(A)
}

#' Generate Pocock-Simon minimization treatment assignments
#'
#' @param z The strata design matrix
#' @param treat A vector of length k (the number of treatment arms), which labels the treatment arms being compared.
#' @param pi A vector of length k (the number of treatment arms), which indicates the allocation ratio, e.g., c(1,1,1) for equal allocation with three treatment arms.
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
#' @export
#'
#' @examples
#'
#' # Create strata variables
#' library(fastDummies)
#' library(dplyr)
#'
#' x <- runif(100)
#' z <- cut(x, breaks=c(0, 0.25, 0.5, 0.75, 1.0))
#' z <- dummy_cols(z)
#' A <- PS(
#'   z=z[, 2:5],
#'   treat=c(0, 1, 2),
#'   ratio=c(1, 1, 1),
#'   imb_measure="Range"
#' )
PS <- function(z, treat, ratio, imb_measure, p_bc=0.8){

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
          sd_level <- apply(num_level, 1, sd)
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



