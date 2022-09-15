
#' @import dplyr
#' @import stats
append.design.matrix <- function(df, covnames){

  formula <- as.formula(
    paste0("~ 0 + ", paste0(covnames, collapse="+"))
  )
  # Get design matrix based on formula
  mat <- stats::model.matrix(formula, df)

  # Center each column of the design matrix
  # This shouldn't do anything for the covariates, because
  # they should already be centered in the .make.data function,
  # but this will center the strata variables also, if they are included.
  mat <- data.frame(mat) %>%
    mutate(across(
      .cols=everything(),
      .fns=list(center=~scale(., center=TRUE, scale=FALSE)),
      .names="{col}_{fn}"
    ))
  browser()
  # Take only the centerd columns back into the data frame
  centered_names <- paste0(covnames, "_center")
  df <- cbind(df, mat %>% select(centered_names))

  return(df)
}

#' @import broom
regress.to.Ohat <- function(df){
  browser()
  res <- df %>%
    group_by(.data$trt1) %>%
    group_modify(~broom::tidy(
      lm(O.hat ~ 0 + ., data=.x %>%
           select(.data$O.hat, model_z1, model_z2))
    ))
}

adjust.LogRank <- function(model, data){

  # Creates data
  df <- create.tte.df(model, data)
  # Risk set and failure time times
  df <- process.tte.df(df, ref_arm=model$ref_arm)
  df <- df %>% mutate(lin_pred=0)

  df <- get.ordered.data(df, ref_arm=model$ref_arm)

  # Get names of covariates if desire adjustment
  if(model$adj_strata | model$adj_cov){
    # Add covariate names in order
    # with the strata covariates first
    # so that if there is collinearity, they
    # are chosen first over the x covariates.
    covnames <- c()
    if(model$adj_strata){
      covnames <- c(covnames, carcov_z)
    }
    if(model$adj_cov){
      covnames <- c(covnames, names(data$covariate))
    }
    # Get centered versions of variables
    df <- append.design.matrix(df, covnames)
    # Regress to Ohat
    df <- regress.to.Ohat(df)
  } else {
    df <- df %>%
      mutate(
        adjust1 = 0,
        adjust0 = 0,
        bsigb   = 0
      )
  }


  # If desire adjustment
  # (either with covariates or stratification variable),
  # regress to Ohat using linear regression.
  if(model$adj_strata | model$adj_cov){

  } else {

  }
  browser()

  result <- list(
    numerator=numerator,
    denominator=denominator,
    statistic=statistic
  )

  return(
    structure(
      class="CoxScoreResult",
      list(result=result, settings=model, data=data)
    )
  )

  result <- list()

  return(
    structure(
      class="LogRankResult",
      list(result=result, settings=model, data=data)
    )
  )
}
