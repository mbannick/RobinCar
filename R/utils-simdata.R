#' @title Data generation function from JRSS-B paper
#' @param n total number of subjects to be generated
#' @param theta true treatment effect
#' @param randomization randomization method in c("SR","CABC","permuted_block","minimization","urn")
#' @param p_trt proportion of treatment arm
#' @param case simulation case in the paper
#' @export
#' @returns A data frame with the following columns:
#'
#' \item{t}{event time}
#' \item{delta}{event indicator}
#' \item{I1}{assignment to treatment group 1}
#' \item{I0}{assignment to treatment group 0}
#' \item{model_z1, model_z2}{covariates}
#' \item{car_strata1, ...}{strata variables}
data_gen <- function(n, theta, randomization, p_trt,
                     case=c("case1", "case2", "case3", "case4", "case5")){
  if(case=="case1") {
    # multiple z
    z1 <- stats::rmultinom(n, 1, c(0.5, 0.5))
    z2 <- stats::rmultinom(n, 1, c(0.4, 0.3, 0.3))
    minimization_z <- t(rbind(z1, z2))
    n_strata <- dim(z1)[1] * dim(z2)[1]
    strata_z <- matrix(nrow = n, ncol = n_strata)
    ind <- 1
    for (i in 1:dim(z1)[1]) {
      for (j in 1:dim(z2)[1]) {
        strata_z[, ind] <- z1[i, ] * z2[j, ]
        ind <- ind + 1
      }
    }
    I <-
      treatment_assignment(n, strata_z, minimization_z, randomization, p_trt)
    lambda0 <- log(2) / 12
    cumuhazard <- stats::rexp(n)
    HR <- exp(I * theta + 1.5 * z1[1, ] - z2[1, ] - 0.5 * z2[2, ])
    t.star <- cumuhazard / lambda0 / HR
    C <- stats::runif(n, 20, 40)
    t <- pmin(t.star, C)
    delta <- 1 * (t.star <= C)
    data.simu <-
      data.frame(t, delta, I, 1 - I, z1[1, ], z2[1, ], z2[2, ], strata_z, minimization_z)[1:n, ]
    names(data.simu) <-
      c(
        "t",
        "delta",
        "I1",
        "I0",
        "model_Z1",
        "model_Z21",
        "model_Z22",
        paste0("car_strata", 1:n_strata),
        paste0("minimization", 1:5)
      )
    return(data.simu)
  }
  else if (case == "case2") {
    # continuous z
    z1 <- stats::rmultinom(n, 1, c(0.5, 0.5))
    row.names(z1) <- paste0("z1_", 1:2)
    n_category <- 16
    tmp <- discretize_z(n, 0, 1, n_category)
    z2_continous <- tmp$Z
    z2 <- t(tmp$Z_dis)

    minimization_z <- t(rbind(z1, z2))
    minimization_n_col <- ncol(minimization_z)

    n_strata <- dim(z1)[1] * dim(z2)[1]
    strata_z <- matrix(nrow = n, ncol = n_strata)
    ind <- 1
    for (i in 1:dim(z1)[1]) {
      for (j in 1:dim(z2)[1]) {
        strata_z[, ind] <- z1[i, ] * z2[j, ]
        ind <- ind + 1
      }
    }
    I <-
      treatment_assignment(n, strata_z, minimization_z, randomization, p_trt)
    lambda0 <- log(2) / 12
    cumuhazard <- stats::rexp(n)
    HR <- exp(I * theta - 1.5 * z1[1, ] + 0.5 * z2_continous ^ 2)
    t.star <- cumuhazard / lambda0 / HR
    C <- stats::runif(n, 10, 40)
    t <- pmin(t.star, C)
    delta <- 1 * (t.star <= C)
    data.simu <-
      data.frame(t, delta, I, 1 - I, z1[1, ], z2_continous ^ 2, strata_z, minimization_z)[1:n, ]
    names(data.simu) <-
      c(
        "t",
        "delta",
        "I1",
        "I0",
        paste0("model_z", 1:2),
        paste0("car_strata", 1:n_strata),
        paste("minimization", 1:minimization_n_col)
      )
    return(data.simu)
  }
  else if (case == "case3") {
    # square (mis)
    z1 <- stats::rmultinom(n, 1, c(0.5, 0.5))
    row.names(z1) <- paste0("z1_", 1:2)
    n_category <- 4
    tmp <- discretize_z(n, 0, 1, n_category)
    z2_continous <- tmp$Z
    z2 <- t(tmp$Z_dis)

    minimization_z <- t(rbind(z1, z2))
    minimization_n_col <- ncol(minimization_z)

    n_strata <- dim(z1)[1] * dim(z2)[1]
    strata_z <- matrix(nrow = n, ncol = n_strata)
    ind <- 1
    for (i in 1:dim(z1)[1]) {
      for (j in 1:dim(z2)[1]) {
        strata_z[, ind] <- z1[i, ] * z2[j, ]
        ind <- ind + 1
      }
    }
    I <-
      treatment_assignment(n, strata_z, minimization_z, randomization, p_trt)
    lambda0 <- log(2) / 12
    cumuhazard <- stats::rexp(n)
    HR <- exp(I * theta - 0.5 * z1[1, ] + 1.5 * z2_continous ^ 2)
    t.star <- cumuhazard / lambda0 / HR
    C <- 10 + stats::rexp(2 * z1[1, ])
    t <- pmin(t.star, C)
    delta <- 1 * (t.star <= C)
    data.simu <-
      data.frame(t, delta, I, 1 - I, z1[1, ], z2_continous, strata_z, minimization_z)[1:n, ]
    names(data.simu) <-
      c(
        "t",
        "delta",
        "I1",
        "I0",
        paste0("model_z", 1:2),
        paste0("car_strata", 1:n_strata),
        paste("minimization", 1:minimization_n_col)
      )
    return(data.simu)
  }
  else if (case == "case4") {
    # non-PH (mis), C unif
    n_category <- 4
    tmp <- discretize_z(n, 0, 1, n_category)
    z1_continuous <- tmp$Z
    z1 <- t(tmp$Z_dis)
    minimization_z <- t(z1)
    minimization_n_col <- ncol(minimization_z)
    n_strata <- dim(z1)[1]
    strata_z <- t(z1)
    I <-
      treatment_assignment(n, strata_z, minimization_z, randomization, p_trt)
    t.star <- exp(theta * I + 1.5 * z1_continuous) + stats::rexp(n, 1)
    C <- stats::runif(n, 10, 20)
    t <- pmin(t.star, C)
    delta <- 1 * (t.star <= C)
    data.simu <-
      data.frame(t, delta, I, 1 - I, z1_continuous, strata_z, minimization_z)[1:n, ]
    names(data.simu) <-
      c(
        "t",
        "delta",
        "I1",
        "I0",
        "model_Z1",
        paste0("car_strata", 1:n_strata),
        paste0("minimization", 1:minimization_n_col)
      )
    return(data.simu)
  }
  else if (case == "case5") {
    # non-PH (mis), C depends on I
    n_category <- 4
    tmp <- discretize_z(n, 0, 1, n_category)
    z1_continuous <- tmp$Z
    z1 <- t(tmp$Z_dis)
    minimization_z <- t(z1)
    minimization_n_col <- ncol(minimization_z)
    n_strata <- dim(z1)[1]
    strata_z <- t(z1)
    I <-
      treatment_assignment(n, strata_z, minimization_z, randomization, p_trt)
    t.star <- exp(theta * I + 1.5 * z1_continuous) + stats::rexp(n, 1)
    C <- 3 - 3 * I + stats::rexp(n, 1)
    t <- pmin(t.star, C)
    delta <- 1 * (t.star <= C)
    data.simu <-
      data.frame(t, delta, I, 1 - I, z1_continuous, strata_z, minimization_z)[1:n, ]
    names(data.simu) <-
      c(
        "t",
        "delta",
        "I1",
        "I0",
        "model_Z1",
        paste0("car_strata", 1:n_strata),
        paste0("minimization", 1:minimization_n_col)
      )
    return(data.simu)
  }
}

#' @title Data generation function from covariate adjusted log-rank paper
#' @inheritParams data_gen
#' @param blocksize block size for permuted block design
#' @export
#'
#' @returns A data frame with the following columns:
#'
#' \item{t}{event time}
#' \item{delta}{event indicator}
#' \item{I1}{assignment to treatment group 1}
#' \item{I0}{assignment to treatment group 0}
#' \item{model_w3}{covariates}
#' \item{car_strata1, ...}{strata variables}
data_gen2 <-
  function(n,
           theta,
           randomization,
           p_trt,
           case = c("case1", "case2", "case3", "case4"),
           blocksize = 4) {
    if (case == "case1") {
      # Cox failure time and uniform censoring
      tmp <- discretize_z(n, 0, 1, n_category = 2)
      w1 <- tmp$Z
      z1 <- t(tmp$Z_dis)
      tmp <- discretize_z(n, 0, 1, n_category = 3)
      w2 <- tmp$Z
      z2 <- t(tmp$Z_dis)
      w3 <- stats::rnorm(n)
      z1.f <- as.factor(z1[1, ] + 2 * z1[2, ])
      z2.f <- as.factor(z2[1, ] + 2 * z2[2, ] + 3 * z2[3, ])
      minimization_z <- data.frame(z1.f, z2.f)

      n_strata <- dim(z1)[1] * dim(z2)[1]
      strata_z <- matrix(nrow = n, ncol = n_strata)
      ind <- 1
      for (i in 1:dim(z1)[1]) {
        for (j in 1:dim(z2)[1]) {
          strata_z[, ind] <- z1[i, ] * z2[j, ]
          ind <- ind + 1
        }
      }
      I <-
        treatment_assignment(n,
                             strata_z,
                             minimization_z,
                             randomization,
                             p_trt,
                             blocksize)
      lambda0 <- log(2)
      cumuhazard <- stats::rexp(n)
      HR <- exp(I * theta + 0.5 * w1 + 0.5 * w2 + 0.5 * w3)
      t.star <- cumuhazard / lambda0 / HR
      C <- stats::runif(n, 10, 40)
      t <- pmin(t.star, C)
      delta <- 1 * (t.star <= C)
      data.simu <- data.frame(t, delta, I, 1 - I, w3, strata_z)[1:n, ]
      names(data.simu) <-
        c("t",
          "delta",
          "I1",
          "I0",
          "model_w3",
          paste0("car_strata", 1:n_strata))
      return(data.simu)
    } else if (case == "case2") {
      # Cox failure time and censoring depends on I
      tmp <- discretize_z(n, 0, 1, n_category = 2)
      w1 <- tmp$Z
      z1 <- t(tmp$Z_dis)
      tmp <- discretize_z(n, 0, 1, n_category = 3)
      w2 <- tmp$Z
      z2 <- t(tmp$Z_dis)
      w3 <- stats::rnorm(n)
      z1.f <- as.factor(z1[1, ] + 2 * z1[2, ])
      z2.f <- as.factor(z2[1, ] + 2 * z2[2, ] + 3 * z2[3, ])
      minimization_z <- data.frame(z1.f, z2.f)

      n_strata <- dim(z1)[1] * dim(z2)[1]
      strata_z <- matrix(nrow = n, ncol = n_strata)
      ind <- 1
      for (i in 1:dim(z1)[1]) {
        for (j in 1:dim(z2)[1]) {
          strata_z[, ind] <- z1[i, ] * z2[j, ]
          ind <- ind + 1
        }
      }
      I <-
        treatment_assignment(n,
                             strata_z,
                             minimization_z,
                             randomization,
                             p_trt,
                             blocksize)
      lambda0 <- log(2)
      cumuhazard <- stats::rexp(n)
      HR <- exp(I * theta + 0.5 * w1 + 0.5 * w2 + 0.5 * w3)
      t.star <- cumuhazard / lambda0 / HR
      C <- 3 - 3 * I + stats::rexp(n, 2)
      t <- pmin(t.star, C)
      delta <- 1 * (t.star <= C)
      data.simu <- data.frame(t, delta, I, 1 - I, w3, strata_z)[1:n, ]
      names(data.simu) <-
        c("t",
          "delta",
          "I1",
          "I0",
          "model_w3",
          paste0("car_strata", 1:n_strata))
      return(data.simu)
    } else if (case == "case3") {
      # Non-Cox failure time and uniform censoring
      tmp <- discretize_z(n, 0, 1, n_category = 2)
      w1 <- tmp$Z
      z1 <- t(tmp$Z_dis)
      tmp <- discretize_z(n, 0, 1, n_category = 3)
      w2 <- tmp$Z
      z2 <- t(tmp$Z_dis)
      w3 <- stats::rnorm(n)
      z1.f <- as.factor(z1[1, ] + 2 * z1[2, ])
      z2.f <- as.factor(z2[1, ] + 2 * z2[2, ] + 3 * z2[3, ])
      minimization_z <- data.frame(z1.f, z2.f)

      n_strata <- dim(z1)[1] * dim(z2)[1]
      strata_z <- matrix(nrow = n, ncol = n_strata)
      ind <- 1
      for (i in 1:dim(z1)[1]) {
        for (j in 1:dim(z2)[1]) {
          strata_z[, ind] <- z1[i, ] * z2[j, ]
          ind <- ind + 1
        }
      }
      I <-
        treatment_assignment(n,
                             strata_z,
                             minimization_z,
                             randomization,
                             p_trt,
                             blocksize)
      t.star <- exp(I * theta + 0.5 * w1 + 0.5 * w2 + 0.5 * w3) + stats::rexp(n, 1)
      C <- stats::runif(n, 10, 40)
      t <- pmin(t.star, C)
      delta <- 1 * (t.star <= C)
      data.simu <- data.frame(t, delta, I, 1 - I, w3, strata_z)[1:n, ]
      names(data.simu) <-
        c("t",
          "delta",
          "I1",
          "I0",
          "model_w3",
          paste0("car_strata", 1:n_strata))
      return(data.simu)
    } else if (case == "case4") {
      # Non-Cox failure time and censoring depends on I
      tmp <- discretize_z(n, 0, 1, n_category = 2)
      w1 <- tmp$Z
      z1 <- t(tmp$Z_dis)
      tmp <- discretize_z(n, 0, 1, n_category = 3)
      w2 <- tmp$Z
      z2 <- t(tmp$Z_dis)
      w3 <- stats::rnorm(n)
      z1.f <- as.factor(z1[1, ] + 2 * z1[2, ])
      z2.f <- as.factor(z2[1, ] + 2 * z2[2, ] + 3 * z2[3, ])
      minimization_z <- data.frame(z1.f, z2.f)

      n_strata <- dim(z1)[1] * dim(z2)[1]
      strata_z <- matrix(nrow = n, ncol = n_strata)
      ind <- 1
      for (i in 1:dim(z1)[1]) {
        for (j in 1:dim(z2)[1]) {
          strata_z[, ind] <- z1[i, ] * z2[j, ]
          ind <- ind + 1
        }
      }
      I <-
        treatment_assignment(n,
                             strata_z,
                             minimization_z,
                             randomization,
                             p_trt,
                             blocksize)
      t.star <- exp(I * theta + 0.5 * w1 + 0.5 * w2 + 0.5 * w3) + stats::rexp(n, 1)
      C <- 3 - 3 * I + stats::rexp(n, 2)
      t <- pmin(t.star, C)
      delta <- 1 * (t.star <= C)
      data.simu <- data.frame(t, delta, I, 1 - I, w3, strata_z)[1:n, ]
      names(data.simu) <-
        c("t",
          "delta",
          "I1",
          "I0",
          "model_w3",
          paste0("car_strata", 1:n_strata))
      return(data.simu)
    } else if (case == "case5") {
      # Non-Cox failure time and censoring depends on I
      tmp <- discretize_z(n, 0, 1, n_category = 2)
      w1 <- tmp$Z
      z1 <- t(tmp$Z_dis)
      z1.f <- as.factor(z1[1, ] + 2 * z1[2, ])
      minimization_z <- data.frame(z1.f)
      n_strata <- dim(z1)[1]
      strata_z <- t(z1)
      I <-
        treatment_assignment(n,
                             strata_z,
                             minimization_z,
                             randomization,
                             p_trt,
                             blocksize)
      t.star <- exp(theta * I + 1.5 * w1) + stats::rexp(n, 1)
      C <- stats::rgamma(n, shape = 5) * (1 - I) + stats::rgamma(n, shape = 1) * I
      t <- pmin(t.star, C)
      delta <- 1 * (t.star <= C)
      data.simu <- data.frame(t, delta, I, 1 - I, w1, strata_z)[1:n, ]
      names(data.simu) <-
        c("t",
          "delta",
          "I1",
          "I0",
          "model_w1",
          paste0("car_strata", 1:n_strata))
      return(data.simu)
    } else if (case == "case4_test") {
      # Non-Cox failure time and censoring depends on I
      tmp <- discretize_z(n, 0, 1, n_category = 2)
      w1 <- tmp$Z
      z1 <- t(tmp$Z_dis)
      tmp <- discretize_z(n, 0, 1, n_category = 3)
      w2 <- tmp$Z
      z2 <- t(tmp$Z_dis)
      w3 <- stats::rnorm(n)
      z1.f <- as.factor(z1[1, ] + 2 * z1[2, ])
      z2.f <- as.factor(z2[1, ] + 2 * z2[2, ] + 3 * z2[3, ])
      minimization_z <- data.frame(z1.f, z2.f)

      n_strata <- dim(z1)[1] * dim(z2)[1]
      strata_z <- matrix(nrow = n, ncol = n_strata)
      ind <- 1
      for (i in 1:dim(z1)[1]) {
        for (j in 1:dim(z2)[1]) {
          strata_z[, ind] <- z1[i, ] * z2[j, ]
          ind <- ind + 1
        }
      }

      fz <- numeric(n)
      for (i in 1:n_strata) {
        fz <- fz + strata_z[, i] * i
      }

      I <-
        treatment_assignment(n,
                             strata_z,
                             data.frame(as.factor(fz)),
                             randomization,
                             p_trt,
                             blocksize)
      t.star <- exp(I * theta + 0.5 * w1 + 0.5 * w2 + 0.5 * w3) + stats::rexp(n, 1)
      C <- 3 - 3 * I + stats::rexp(n, 2)
      t <- pmin(t.star, C)
      delta <- 1 * (t.star <= C)
      data.simu <- data.frame(t, delta, I, 1 - I, w3, strata_z)[1:n, ]
      names(data.simu) <-
        c("t",
          "delta",
          "I1",
          "I0",
          "model_w3",
          paste0("car_strata", 1:n_strata))
      return(data.simu)
    }
  }

discretize_z <- function(n, z_mean, z_sd, n_category) {
  # only normal
  Z <- stats::rnorm(n, z_mean, z_sd)
  p_cut <- seq(0, 1, length.out = (n_category + 1))
  q_cut <- sapply(p_cut, function(x)
    stats::qnorm(x, z_mean, z_sd))
  Z_dis_scale <- as.factor(sapply(Z, function(x)
    max(which(q_cut < x))))
  Z_model_matrix <- model.matrix( ~ 0 + Z_dis_scale)
  return(list(Z = Z, Z_dis = Z_model_matrix))
}

SR <- function(n, p_trt) {
  I <- rbinom(n, 1, p_trt)
  return(I)
}


BC <-
  function(n_within_strata, BCD_p = 2 / 3) {
    # Efron's biased coin randomization (Efron, 1971)
    stopifnot(BCD_p > 1 / 2)
    I <- numeric(n_within_strata)
    I[1] <- rbinom(1, 1, 1 / 2)
    D <- 2 * I[1] - 1
    for (i in 2:n_within_strata) {
      if (D > 0)
        I[i] <- stats::rbinom(1, 1, (1 - BCD_p))
      else if (D < 0)
        I[i] <- stats::rbinom(1, 1, BCD_p)
      else if (D == 0)
        I[i] <- stats::rbinom(1, 1, 1 / 2)
      D <- D + 2 * I[i] - 1
    }
    return(I)
  }

CABC <- function(z, BCD_p = 2 / 3) {
  n_strata <- dim(z)[2]
  n <- dim(z)[1]
  n_each_strata <- apply(z, 2, sum)
  I <- numeric(n)
  for (s in 1:n_strata) {
    if (n_each_strata[s] == 0) {
      next
    }
    else if (n_each_strata[s] == 1) {
      I[z[, s] == 1] <- stats::rbinom(1, 1, 0.5)
    }
    else{
      I[z[, s] == 1] <- BC(n_each_strata[s], BCD_p)
    }
  }
  return(I)
}

urn_design <-
  function(n_within_strata, omega, gamma, beta) {
    # Wei's urn design (Wei, 1978, JASA, Vol.73, No. 363)
    I <- numeric(n_within_strata)
    n_balls_0 <- omega
    n_balls_1 <- omega
    for (i in 1:n_within_strata) {
      p <- n_balls_1 / (n_balls_1 + n_balls_0)
      I[i] <- stats::rbinom(1, 1, p)
      n_balls_0 <- n_balls_0 + (1 - I[i]) * gamma + I[i] * beta
      n_balls_1 <- n_balls_1 + (1 - I[i]) * beta + I[i] * gamma
    }
    return(I)
  }

strata_urn_design <- function(z) {
  omega <- 1
  gamma <- 0
  beta <- 1
  n_strata <- dim(z)[2]
  n <- dim(z)[1]
  n_each_strata <- apply(z, 2, sum)
  I <- numeric(n)
  for (i in 1:n_strata) {
    if (n_each_strata[i] == 0) {
      next
    }
    else if (n_each_strata[i] == 1) {
      I[z[, i] == 1] <- stats::rbinom(1, 1, 0.5)
    }
    else{
      I[z[, i] == 1] <- urn_design(n_each_strata[i], omega, gamma, beta)
    }
  }
  return(I)
}


permuted_block <- function(z, blocksize, p_trt) {
  n_trt_blk <- blocksize * p_trt
  n_ctl_blk <- blocksize * (1 - p_trt)
  if (!(floor(n_ctl_blk) == n_ctl_blk &
        floor(n_trt_blk) == n_trt_blk)) {
    stop("number of treatment and control within block are not integers!")
  }
  if (blocksize)
    n_strata <- dim(z)[2]
  n <- dim(z)[1]
  n_each_strata <- apply(z, 2, sum)
  I <- numeric(n)
  for (i in 1:n_strata) {
    I_tmp <- 0
    n_block_tmp <- ceiling(n_each_strata[i] / blocksize)
    for (b in 1:n_block_tmp) {
      block_permuted_tmp <- sample(c(rep(0, n_ctl_blk), rep(1, n_trt_blk)))
      I_tmp <- c(I_tmp, block_permuted_tmp)
    }
    I_tmp <- I_tmp[-1]
    I[z[, i] == 1] <- I_tmp[1:n_each_strata[i]]
  }
  return(I)
}

treatment_assignment <-
  function(n,
           strata_z,
           minimization_z,
           randomization,
           p_trt,
           blocksize = 4) {
    if (randomization == "SR") {
      I <- SR(n, p_trt)
    } else if (randomization == "CABC") {
      if (p_trt != 1 / 2) {
        stop("For now, the proportion of treatment under CABC has to be 1/2.")
      }
      I <- CABC(strata_z)
    } else if (randomization == "permuted_block") {
      if (p_trt != 1 / 2) {
        stop("For now, the proportion of treatment under permuted_block has to be 1/2.")
      }
      I <- permuted_block(strata_z, blocksize, p_trt)
    } else if (randomization == "minimization") {
      if (p_trt != 1 / 2) {
        stop("For now, the proportion of treatment under CABC has to be 1/2.")
      }
      I <- minimization(minimization_z)
    } else if (randomization == "urn") {
      if (p_trt != 1 / 2) {
        stop("For now, the proportion of treatment under CABC has to be 1/2.")
      }
      I <- strata_urn_design(strata_z)
    }
    return(I)
  }

minimization <- function(z,
                         imbalance_measure = "square",
                         BCD_p = 2 / 3) {
  z <- as.matrix(z)
  n <- dim(z)[1]
  I <- numeric(n)
  I[1] <- rbinom(1, 1, 1 / 2)
  D <- matrix(0, nrow = dim(z)[1], ncol = dim(z)[2])
  D[1, ] <- (2 * I[1] - 1) * z[1, ]
  for (i in 2:n) {
    D_before_i <- D[i - 1, ]
    D_potential_1 <- D_before_i + z[i, ]
    D_potential_0 <- D_before_i - z[i, ]
    if (imbalance_measure == "square") {
      imb_potential_1 <- sum(D_potential_1 ^ 2)
      imb_potential_0 <- sum(D_potential_0 ^ 2)
    } else if (imbalance_measure == "absolute") {
      imb_potential_1 <- sum(abs(D_potential_1))
      imb_potential_0 <- sum(abs(D_potential_0))
    }

    if (imb_potential_1 > imb_potential_0) {
      # imbalance is larger if assinged to I=1
      I[i] <- rbinom(1, 1, (1 - BCD_p)) #BCD_p>1/2
    } else if (imb_potential_1 < imb_potential_0) {
      I[i] <- rbinom(1, 1, BCD_p)
    } else if (imb_potential_1 == imb_potential_0) {
      I[i] <- rbinom(1, 1, 0.5)
    }
    D[i, ] <- D_potential_1 * I[i] + D_potential_0 * (1 - I[i])
  }
  return(I)
}
