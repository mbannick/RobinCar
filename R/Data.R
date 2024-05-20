# Functions for validating data and creating data classes
# to be used across all analysis methods.

.check.response <- function(data){
  if(!is.numeric(data$response)){
    return("Response column must be numeric.")
  }
}

.check.event <- function(data){
  if(!all(data$event %in% c(0, 1))){
    return("Event column must contain only 0 and 1.")
  }
}

.check.treat <- function(data){
  if(!is.factor(data$treat)){
    return("Treatment column must be a factor variable.")
  }
}

.check.attributes <- function(x, ...){
  required <- c(...)
  existing <- names(x)
  missing <- which(!required %in% existing)

  if(length(missing > 0)){
    return(paste0("Missing data attributes ", paste(missing)))
  }
}

validate <- function (data, ...) {
  UseMethod("validate", data)
}

#' @exportS3Method
validate.RoboDataLinear <- function(data, ...){

  errors <- character()
  errors <- c(errors, .check.attributes(data, "treat", "response"))
  errors <- c(errors, .check.response(data))
  errors <- c(errors, .check.treat(data))

  .return.error(errors)
}

#' @exportS3Method
validate.RoboDataGLM <- function(data, ...){

  errors <- character()
  .return.error(errors)

}

#' @exportS3Method
validate.RoboDataSL <- function(data, ...){

  errors <- character()
  .return.error(errors)

}

#' @exportS3Method
validate.RoboDataTTE <- function(data, ref_arm, ...){

  errors <- character()
  errors <- c(errors, .check.attributes(data, "treat", "response", "event"))
  errors <- c(errors, .check.response(data))
  errors <- c(errors, .check.event(data))
  errors <- c(errors, .check.treat(data))

  if(length(data$treat_levels) != 2){
    errors <- c(errors, "Need to have only two treatment levels.")
  }
  if(!is.null(ref_arm)){
    if(!ref_arm %in% data$treat_levels){
      errors <- c(
        errors,
        message(
          "The reference group %s is not in the treatment levels %s",
          ref_arm,
          data$treat_levels
        )
      )
    }
  }

  .return.error(errors)
}

.edit.formula <- function(form, response_col, treat_col){

  # Convert formula to character
  form <- deparse(stats::as.formula(form))

  # Replace formula with proper names for response and treatment
  # -- response should only be at beginning
  newform <- gsub(paste0("^", response_col), "response", form)
  # -- only replace for the entire word
  newform <- gsub(paste0("\\b", treat_col, "\\b"), "treat", newform)

  # Extract covariate names from formula
  vars <- strsplit(newform, c("(~|\\+|:|\\*)"))[[1]]
  vars <- sapply(vars, function(x) gsub(" ", "", x))
  vars <- sapply(vars, function(x) gsub("\\(", "", x))
  vars <- sapply(vars, function(x) gsub("\\)", "", x))

  # Check to make sure that response and treatment columns
  # are in the formula
  if(vars[1] != "response") stop("Must include response variable on LHS of formula.")
  if(!"treat" %in% vars) stop("Must include treatment variable in formula.")

  vars <- vars[!vars %in% c("response", "treat", "0", "1")]
  vars <- unique(vars)

  return(list(newform, vars))
}

# Creates formula for robincar_linear
# based on adj_method=c("ANOVA", "ANCOVA", "ANHECOVA").
.create.formula <- function(adj_method, response_col, treat_col, covariate_cols){

  if(adj_method == "ANOVA"){
    form <- paste0(response_col, " ~ ", treat_col)
  } else {

    if(is.null(covariate_cols)) stop("Must provide covariates if
                                     adjustment method is ANCOVA or ANHECOVA.")
    covariates <- paste(covariate_cols, collapse=" + ")

    if(adj_method == "ANCOVA"){
      form <- paste0(response_col, " ~ ", treat_col, " + ", covariates)
    } else if(adj_method == "ANHECOVA"){
      form <- paste0(response_col, " ~ ", treat_col, " * (", covariates, ")")
    } else {
      stop("Unrecognized adjustment method.")
    }

  }
  return(form)
}

.df.toclass <- function(df, classname, ...){

  data <- list()
  atts <- list(...)

  for(i in 1:length(atts)){

    att_name <- names(atts)[i]
    att <- atts[[i]]

    # save the original names of the attributes
    data[[att_name]] <- att

    if(att_name == "formula"){

      if(is.null(att)) next

      response_col <- atts["response_col"]
      treat_col <- atts["treat_col"]

      if(is.null(response_col)) stop("Must provide a response column name.")
      if(is.null(treat_col)) stop("Must provide a treatment column name.")

      # Edit formula and get names of formula variables
      # The function below returns a list with [[1]] formula and [[2]] formula variable names
      forms <- .edit.formula(att, response_col, treat_col)

      # Include the covariates that will be needed for the formula
      # in a formula vector
      data[["formula"]] <- forms[[1]]

    } else if(grepl("col", att_name)){

      att_name <- gsub("_cols", "", att_name)
      att_name <- gsub("_col", "", att_name)

      data[[att_name]] <- df[att]

    } else {

      stop(paste0("Unrecognized column arguments ",
                  att_name, ". All columns must have
                                            _col or _cols as a suffix."))

    }
  }
  argnames <- list(...)

  class(data) <- c("Data", classname)

  return(data)
}

.make.data <- function(df, classname, ...){

  # Convert data frame to object
  data <- .df.toclass(df=df, classname=classname, ...)

  if(!is.null(data$treat)){
    data$treat <- data$treat[[1]]
    if(!is.factor(data$treat)) data$treat <- as.factor(data$treat)
  }
  if(!is.null(data$response)){
    data$response <- data$response[[1]]
  }
  if(!is.null(data$event)){
    data$event <- data$event[[1]]
  }
  if(ncol(data$car_strata) == 0){
    data$car_strata <- NULL
  }
  if(!is.null(data$car_strata)){
    # Create joint car_strata levels
    data$joint_strata <- joint.car_strata(data$car_strata)
    data$joint_strata_levels <- levels(data$joint_strata)

    # Make sure all car_strata are factors
    for(col in colnames(data$car_strata)){
      data$car_strata[col] <- as.factor(data$car_strata[[col]])
    }
  }
  if(!is.null(data$covariate)){
    if(ncol(data$covariate) == 0){
      data$covariate <- NULL
    }
  }

  # Save original data frame
  data$df <- df

  # Add additional data attributes
  data$n <- nrow(df)
  data$k <- length(levels(data$treat))
  data$treat_levels <- levels(data$treat)

  # Estimate treatment allocation proportion
  data$pie <- table(data$treat) %>%
    proportions() %>%
    as.matrix()

  return(data)
}
