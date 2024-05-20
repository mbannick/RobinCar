# Defines allowable arguments to the functions available to users

.check.options <- function(name, var, options){
  if(!var %in% options) stop(paste0(name, " must be one of: ",
                                    paste0(options, collapse = ", ")))
}

.check.car_scheme <- function(car_scheme, car_strata_cols){
  OPTIONS <- c("simple",
               "permuted-block",
               "pocock-simon",
               "biased-coin",
               "urn")
  .check.options("car_scheme", car_scheme, OPTIONS)

  # Check to make sure that existence of
  # strata matches with the car scheme
  if(!is.null(car_strata_cols)){
    if(car_scheme == "simple") .z.exist.warn()
  } else {
    if(car_scheme != "simple") .z.miss.err()
  }
}

.check.adj_method.logrank <- function(adj_method){
  OPTIONS <- c("CL",
               "CSL")
  .check.options("adj_method", adj_method, OPTIONS)
}

.check.adj_method.glm <- function(adj_method){
  OPTIONS <- c("heterogeneous",
               "homogeneous")
  .check.options("adj_method", adj_method, OPTIONS)
}

.check.adj_method.linear <- function(adj_method){
  OPTIONS <- c("ANOVA",
               "ANCOVA",
               "ANHECOVA")
  .check.options("adj_method", adj_method, OPTIONS)
}

.check.vcovHC <- function(vcovHC){
  OPTIONS <- c("HC0",
               "HC1",
               "HC3")
  .check.options("vcovHC", vcovHC, OPTIONS)
}

# .check.sl.libraries <- function(SL_libraries){
#   libs <- invisible(SuperLearner::listWrappers(what="SL"))
#   OPTIONS <- libs[grepl("^SL.", libs)]
#   # for(lib in SL_libraries){
#   #   .check.options("Super Learner libraries", lib, OPTIONS)
#   # }
# }
