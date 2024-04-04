Robin_g <- function(robin.data = data.tmp,
                    car.scheme = "PB", # SR
                    car.z = c("z1", "z2"), # joint levels of 'z' will be used under CAR
                    robin.x = c("x1"),
                    robin.x_to_include_z = TRUE, # default
                    robin.adj_method = c("heterogeneous", "homogeneous"), # option name may be revised, homo ~ ANCOVA-type, hetero ~ ANHECOVA-type
                    robin.formula = as.formula("y ~ z1 + z2 + x1 + A + x1*A"), # if robin.formula is specified, robin.x, robin.x_to_include_z and robin.adj_method are ignored, robin.formula does not work for "minimization"
                    robin.g_family = binomial(link = "logit"), # usually only family is needed, i.e., default canonical link will be applied, similar to glm
                    robin.g_accuracy = 7,
                    robin.vcovHC = c("HC0", "HC1", "HC3")){

  library(tidyverse)

  robin.n <- nrow(robin.data)
  robin.trt_pie <- robin.data %>% group_by(A) %>% summarise(pie_A = n()/robin.n, .groups = "drop")

  # generate joint_z for CAR
  if(car.scheme == "SR"){
    car.z = NULL
    robin.x_to_include_z = FALSE
    robin.g_data <- robin.data
  }else{ # car.scheme == "PB"
    if(is.null(car.z)){ # under CAR, car.z cannot be NULL
      stop("under CAR, please specify 'car.z'.")
    }else{ # car.z is specified and car.scheme is not SR, generate joint levels of car.z
      robin.g_data <- robin.data %>% group_by(across(all_of(car.z))) %>%
        mutate(joint_z = factor(cur_group_id())) %>% ungroup()
    }
  }

  # customized formula or not?
  if(class(robin.formula) == "formula"){
    robin.x <- NULL
    robin.x_to_include_z <- F
    robin.adj_method <- NULL

  }else{

    if(car.scheme == "minimization"){
      # if robin.formula is NOT used, joint_z should at least be included for covariate adjustment under minimization
      if(!robin.x_to_include_z){
        robin.x_to_include_z <- TRUE
        warning("robin.x_to_include_z is set to be 'TRUE'! ")
      }
    }

    if(robin.x_to_include_z){
      robin.x <- c(robin.x, "joint_z")
      robin.g_data <- robin.g_data %>% select(all_of(c("y", "A", car.z, robin.x)))
    }else{
      if(car.scheme == "SR"){
        robin.g_data <- robin.g_data %>% select(all_of(c("y", "A", robin.x)))
      }else{ # car.scheme == "applicable CAR"
        robin.g_data <- robin.g_data %>% select(all_of(c("y", "A", car.z, robin.x, "joint_z")))
      }
    }

  }



  # make sure robin.x does not have elements starting with "g_pred_"
  if(sum(str_detect(colnames(robin.g_data), "g_pred"))>0){
    stop("Please make sure robin.data does not include column name starting with 'g_pred'.")
  }

  # fit glm and get g_pred_k(mu_k_hat)
  if(class(robin.formula) == "formula"){ # this part is not that mature

    robin.g_fit <- glm(as.formula(robin.formula), family = robin.g_family, data = robin.g_data)
    robin.x <- setdiff(as.character(attributes(robin.g_fit$terms)$predvars), c('list', 'y', 'A'))
    robin.g_pred_data <- do.call(cbind.data.frame,
                                 lapply(lapply(robin.trt_pie$A, function(x){cbind.data.frame(A=x, robin.g_data %>% select(all_of(robin.x)))}),
                                        function(x){predict(robin.g_fit, newdata = x, type = "response")}))
    colnames(robin.g_pred_data) <- paste0("g_pred_", robin.trt_pie$A)
    robin.g_data <- cbind.data.frame(robin.g_data, robin.g_pred_data)

  }else if(robin.adj_method == "homogeneous"){

    robin.g_fit <- glm(y ~ A+., family = robin.g_family, data = robin.g_data %>% select(all_of(c("y", "A", robin.x))))
    robin.g_pred_data <- do.call(cbind.data.frame,
                                 lapply(lapply(robin.trt_pie$A, function(x){cbind.data.frame(A=x, robin.g_data %>% select(all_of(robin.x)))}),
                                        function(x){predict(robin.g_fit, newdata = x, type = "response")}))
    colnames(robin.g_pred_data) <- paste0("g_pred_", robin.trt_pie$A)
    robin.g_data <- cbind.data.frame(robin.g_data, robin.g_pred_data)

  }else if(robin.adj_method == "heterogeneous"){

    robin.g_fit <- robin.g_data %>% group_by(A) %>%
      do(g_fit = glm(y ~ ., family = robin.g_family, data = .data %>% select(all_of(c("y", robin.x)))))
    robin.g_pred_data <- do.call(cbind.data.frame,
                                 lapply(robin.g_fit$g_fit,
                                        function(x){predict(x, newdata = robin.g_data %>% select(all_of(robin.x)), type = "response")}))
    names(robin.g_pred_data) <- paste0("g_pred_", robin.g_fit$A)
    robin.g_data <- cbind.data.frame(robin.g_data, robin.g_pred_data)

  }

  # - check predictive unbiasedness
  if(car.scheme == "minimization"){ # check predictive unbiasedness over joint_z
    if(
      all(
        diag(
          robin.g_data %>%
          select(all_of(c("y", "A", "joint_z", paste0("g_pred_", robin.trt_pie$A)))) %>%
          group_by(A, joint_z) %>%
          summarise(across(.cols = everything(), .fns = ~ round(mean(.x), digits = robin.g_accuracy)), .groups = "drop") %>%
          group_by(A) %>%
          summarise(across(.cols = paste0("g_pred_", robin.trt_pie$A), .fns = ~ sum(.x == y)), .groups = "drop") %>%
          select(all_of(paste0("g_pred_", robin.trt_pie$A))) %>%
          as.matrix()
        ) == length(levels(robin.g_data$joint_z))
      )){
      robin.predunbiased <- T
    }else{
      stop("minimization requires estimator to be predictive unbiased across every car_strata of joint z!")
    }

  }else{ # car.scheme != "minimization", i.e., car.scheme %in% SR or commonly used CAR
    # check predictive unbiasedness over all subjects
    if(
      all(
        diag(
          robin.g_data %>%
          select(all_of(c("y", "A", paste0("g_pred_", robin.trt_pie$A)))) %>%
          group_by(A) %>%
          summarise(across(.cols = everything(), .fns = ~ round(mean(.x), digits = robin.g_accuracy)), .groups = "drop") %>%
          group_by(A) %>%
          summarise(across(.cols = paste0("g_pred_", robin.trt_pie$A), .fns = ~ sum(.x == y)), .groups = "drop") %>%
          select(all_of(paste0("g_pred_", robin.trt_pie$A))) %>%
          as.matrix()
        ) == 1)){
      robin.predunbiased <- T
    }else{
      robin.predunbiased <- F
    }
  }

  # - AIPW estimator: g-estimator = AIPW if predictive unbiasedness holds
  robin.g_pred_center <- diag(
    robin.g_data %>%
      group_by(A) %>%
      summarise(across(.cols = paste0("g_pred_", robin.trt_pie$A), .fns = ~ mean(y - .x)), .groups = "drop") %>%
      mutate(A = paste0("g_pred_", A)) %>%
      column_to_rownames("A") %>%
      as.matrix()
  )
  robin.g_pred_center <- as.list(robin.g_pred_center)

  # centering "g_pred_k", if predictive unbiasedness holds, automatically centering by zero
  robin.g_data <- robin.g_data %>% mutate(across(all_of(names(robin.g_pred_center)), ~ .x + robin.g_pred_center[[cur_column()]]))

  robin.est <- apply(robin.g_data %>% select(all_of(paste0("g_pred_", robin.trt_pie$A))), 2, mean)

  if(! robin.predunbiased){ # can be placed in later place.
    warning("The g-estimator does not satisfy predictive unbiasedness, AIPW is reported.")
  }

  # - variance-covariance matrix for AIPW estimator
  robin.varcov.term1 <- robin.g_data %>%
    select(all_of(c("y", "A", paste0("g_pred_", robin.trt_pie$A)))) %>%
    group_by(A) %>%
    summarise(across(.cols = all_of(paste0("g_pred_", robin.trt_pie$A)), .fns = ~ var(y - .x)), .groups = "drop") %>%
    left_join(., robin.trt_pie, by = "A") %>%
    group_by(A) %>%
    summarise(across(.cols = all_of(paste0("g_pred_", robin.trt_pie$A)), .fns = ~ .x/pie_A), .groups = "drop") %>%
    mutate(A = paste0("g_pred_", A)) %>%
    column_to_rownames("A") %>%
    as.matrix()
  robin.varcov.term1 <- diag(diag(robin.varcov.term1))
  colnames(robin.varcov.term1) <- paste0("g_pred_", robin.trt_pie$A)
  rownames(robin.varcov.term1) <- paste0("g_pred_", robin.trt_pie$A)


  robin.varcov.term2 <- robin.g_data %>%
    select(all_of(c("y", "A", paste0("g_pred_", robin.trt_pie$A)))) %>%
    group_by(A) %>%
    summarise(across(.cols = all_of(paste0("g_pred_", robin.trt_pie$A)), .fns = ~ cov(y, .x)), .groups = "drop") %>%
    mutate(A = paste0("g_pred_", A)) %>%
    column_to_rownames("A") %>%
    as.matrix()

  robin.varcov.term3 <- t(robin.varcov.term2)

  robin.varcov.term4 <- robin.g_data %>%
    select(all_of(paste0("g_pred_", robin.trt_pie$A)))
  robin.varcov.term4 <- var(robin.varcov.term4)

  #- robin.varcov.term5 is determined by design
  if(car.scheme %in% c("SR", "minimization")){
    robin.varcov.term5 <- 0
    # by definition, SR > robin.varcov.term5 to be zero
    # by previous check, minimization requires predictive unbiasedness across all car_strata of joint_z

  }else if(car.scheme %in% c("PB", "stratified_biased_coin")){
    robin.varcov.term5.omiga.SR <- diag(robin.trt_pie$pie_A) - robin.trt_pie$pie_A %*% t(robin.trt_pie$pie_A)
    rownames(robin.varcov.term5.omiga.SR) <- paste0("g_pred_", robin.trt_pie$A)
    colnames(robin.varcov.term5.omiga.SR) <- paste0("g_pred_", robin.trt_pie$A)
    robin.varcov.term5.omiga.CAR_by_Z <- diag(0, nrow(robin.varcov.term5.omiga.SR))
    rownames(robin.varcov.term5.omiga.CAR_by_Z) <- paste0("g_pred_", robin.trt_pie$A)
    colnames(robin.varcov.term5.omiga.CAR_by_Z) <- paste0("g_pred_", robin.trt_pie$A)

    robin.varcov.term5.RZ.tmp <- robin.g_data %>%
      select(all_of(c("y", "A", "joint_z", paste0("g_pred_", robin.trt_pie$A)))) %>%
      left_join(., robin.trt_pie, by = "A") %>%
      mutate(across(.cols = paste0("g_pred_", robin.trt_pie$A), .fns = ~ (y - .x)/pie_A)) %>%
      group_by(A, joint_z) %>%
      summarise(p_joint_z = n()/robin.n,
                across(.cols = paste0("g_pred_", robin.trt_pie$A), .fns = ~ mean(.x)),
                .groups = "drop") %>%
      group_by(joint_z) %>%
      group_split()

    robin.varcov.term5.p_joint_z <-
      do.call(
        rbind.data.frame,
        lapply(robin.varcov.term5.RZ.tmp,
               function(x){x %>% group_by(joint_z) %>% summarise(p_joint_z = sum(p_joint_z))})
      )

    robin.varcov.term5 <-
      lapply(robin.varcov.term5.RZ.tmp,
             function(x){
               xI <- x %>% select(A, all_of(paste0("g_pred_", robin.trt_pie$A))) %>%
                 mutate(A = paste0("g_pred_", A)) %>%
                 column_to_rownames("A") %>%
                 as.matrix() # the diagonal terms of xI is robin.varcov.term5.RZ, which is diag(diag(xI))
               xII <- diag(diag(xI)) %*% (robin.varcov.term5.omiga.SR - robin.varcov.term5.omiga.CAR_by_Z) %*% diag(diag(xI))
               colnames(xII) <- paste0("g_pred_", robin.trt_pie$A)
               xII <- cbind.data.frame(joint_z = x$joint_z, as.data.frame(xII)) %>%
                 left_join(., robin.varcov.term5.p_joint_z, by = "joint_z") %>%
                 mutate(across(.cols = paste0("g_pred_", robin.trt_pie$A), .fns = ~ .x*p_joint_z)) %>%
                 select(all_of(paste0("g_pred_", robin.trt_pie$A)))
               return(xII)}
      )

    robin.varcov.term5 <- Reduce("+", robin.varcov.term5) %>% as.matrix() # cautious: do not use "reduce", capital "R" matters

  }else if(car.scheme %in% c("stratified urn design")){ # has not been implemented
    # robin.varcov.term5.omiga.SR <- diag(robin.trt_pie$pie_A) - robin.trt_pie$pie_A %*% t(robin.trt_pie$pie_A)
    # robin.varcov.term5.omiga.CAR_by_Z
  }else{
    stop("There is no such design") # this actually should be checked at the beginning.
  }

  #- put all robin.varcov.terms together
  robin.varcov <- robin.varcov.term1 + robin.varcov.term2 + robin.varcov.term3 - robin.varcov.term4 - robin.varcov.term5

  #--- adjust by HC1 or HC3
  if(length(robin.vcovHC)>1){
    robin.vcovHC.wgt <- 1
  }else if(robin.vcovHC=="HC0"){
    robin.vcovHC.wgt <- 1
  }else if(robin.vcovHC=="HC1"){
    robin.vcovHC.wgt <- 1 # not implemented yet
  }else if(robin.vcovHC=="HC3"){
    robin.vcovHC.wgt <- 1 # not implemented yet
  }else{
    stop("There is no such vcovHC method!")
  }

  robin.varcov <- robin.varcov*robin.vcovHC.wgt/robin.n

  # - summaries return
  robin.est <- data.frame(
    A = sub(".*g_pred_", "", names(robin.est)), # str_extract anything after "g_pred_"
    estimate = robin.est,
    se = sqrt(diag(robin.varcov))
  ) %>%
    mutate(`pval (2-sided)` = 2*pnorm(abs(estimate/se), lower.tail = F))


  robin.return <- list(estimation = robin.est,
                       varcov = robin.varcov,
                       n = robin.n,
                       robin.g_object = list(robin.g_data = robin.g_data,
                                             car.scheme = car.scheme,
                                             car.z = car.z))

  return(robin.return)
}
