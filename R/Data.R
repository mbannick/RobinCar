# Functions for validating data and creating data classes
# to be used across all analysis methods.

.check.response <- function(data){
  if(class(data$response) != "numeric"){
    return("Response column must be numeric.")
  }
}

.check.event <- function(data){
  if(!all(data$event %in% c(0, 1))){
    return("Event column must contain only 0 and 1.")
  }
}

.check.treat <- function(data){
  if(class(data$treat) != "factor"){
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

validate <- function (data) {
  UseMethod("validate", data)
}

validate.RoboDataLinear <- function(data){

  errors <- character()
  errors <- c(errors, .check.attributes(data, "treat", "response"))
  errors <- c(errors, .check.response(data))
  errors <- c(errors, .check.treat(data))

  .return.error(errors)
}

validate.RoboDataGLM <- function(data){

  errors <- character()
  .return.error(errors)

}

validate.RoboDataTTE <- function(data){

  errors <- character()
  errors <- c(errors, .check.attributes(data, "treat", "response", "event"))
  errors <- c(errors, .check.response(data))
  errors <- c(errors, .check.event(data))
  errors <- c(errors, .check.treat(data))

  .return.error(errors)
}

.edit.formula <- function(form, response_col, treat_col){

  # Convert formula to character
  form <- deparse(as.formula(form))

  # Replace formula with proper names for response and treatment
  newform <- gsub(response_col, "response", form)
  newform <- gsub(treat_col, "treat", newform)

  # Extract covariate names from formula
  vars <- strsplit(newform, c("(~|\\+|:|\\*)"))[[1]]
  vars <- sapply(vars, function(x) gsub(" ", "", x))

  # Check to make sure that response and treatment columns
  # are in the formula
  if(vars[1] != "response") stop("Must include response variable on LHS of formula.")
  if(!"treat" %in% vars) stop("Must include treatment variable in formula.")

  vars <- vars[!vars %in% c("response", "treat", "0", "1")]
  vars <- unique(vars)

  return(list(newform, vars))
}

# Takes a data frame and converts it to a list with
# the attributes as specified by the names passed to ...
.df.toclass <- function(df, classname, ...){

  data <- list()
  atts <- list(...)
  for(i in 1:length(atts)){

    att_name <- names(atts)[i]
    att <- atts[[i]]

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
      data[["formula_vars"]] <- df[forms[[2]]]

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
  class(data) <- c("Data", classname)

  return(data)
}

.make.data <- function(df, classname, ...){

  # Convert data frame to object
  data <- .df.toclass(df=df, classname=classname, ...)

  if(!is.null(data$treat)){
    data$treat <- as.factor(as.vector(data$treat[[1]]))
  }
  if(!is.null(data$response)){
    data$response <- as.vector(data$response[[1]])
  }
  if(!is.null(data$strata)){
    # Create joint strata levels
    data$joint_strata <- joint.strata(data$strata)
    data$joint_strata_levels <- levels(data$joint_strata)

    # Make sure all strata are factors
    for(col in colnames(data$strata)){
      data$strata[col] <- as.factor(data$strata[[col]])
    }
  }
  if(!is.null(data$covariate)){
    # Center x variables
    for(col in colnames(data$covariate)){
      data$covariate[col] <- data$covariate[col] - mean(data$covariate[[col]])
    }
  }

  # Add additional data attributes
  data$n <- nrow(df)
  data$k <- length(levels(data$treat))
  data$treat_levels <- levels(data$treat) %>% as.numeric

  # Estimate treatment allocation proportion
  data$pie <- table(data$treat) %>%
    proportions() %>%
    as.matrix()

  return(data)
}
