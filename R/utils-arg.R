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

# Names of columns that are created or overwritten internally when
# building and processing the time-to-event dataset (see adjust-tte.R,
# adjust-logrank.R, adjust-coxscore.R, adjust-covhr.R). User-supplied
# covariates with these names would be silently overwritten or clash.
.TTE.RESERVED.NAMES <- c(
  # create.tte.df
  "treat", "response", "event", "nu_d", "id", "carcov_z", "car_strata",
  # process.tte.df
  "Y", "Y0", "Y1", "trt0", "trt1", "n.events",
  "cum.delta.Y0.over.Ysq", "cum.delta.Y1.over.Ysq",
  "cum.delta.Y0.over.Ysq.theta_L", "cum.delta.Y1.over.Ysq.theta_L",
  # get.ordered.data
  "mu_t", "O.hat", "O.hat0", "O.hat1", "lin_pred",
  "s0", "s1", "s0_seq", "s1_seq",
  "mean_at_risk", "cumsum_at_risk", "mean_at_risk2", "cumsum_at_risk2",
  "O_i", "u_i",
  # calculate.adjustment and test statistics
  "adjust1", "adjust0", "bsigb",
  "uu_cl", "ssig_l", "U_SL_z", "var_adj",
  "cond_var0", "cond_var1", "cond_mean0", "cond_mean1", "prob_z", "z_var"
)

# Prefixes of internally generated columns (design matrix and regression betas)
.TTE.RESERVED.PREFIXES <- c("xmat_", "trt1_")

.check.tte_reserved_names <- function(df, treat_col, response_col, event_col,
                                      car_strata_cols=NULL,
                                      covariate_cols=NULL, id_col=NULL){
  # Columns of df that are used as covariates are copied into the internal
  # dataset, so they must not collide with internally created columns.
  # Other user columns (treat, response, event, id, strata) are renamed
  # internally, but we still check every column of df that is used.
  used <- c(covariate_cols)
  bad <- used[used %in% .TTE.RESERVED.NAMES]
  prefix_pattern <- paste0("^(", paste0(.TTE.RESERVED.PREFIXES, collapse="|"), ")")
  bad <- unique(c(bad, used[grepl(prefix_pattern, used)]))
  if(length(bad) > 0){
    stop(paste0(
      "The following covariate column(s) in df use reserved names that are ",
      "overwritten internally and would produce incorrect results: ",
      paste0(bad, collapse=", "),
      ". Please rename these columns. Reserved names include: ",
      paste0(.TTE.RESERVED.NAMES, collapse=", "),
      "; and names beginning with: ",
      paste0(.TTE.RESERVED.PREFIXES, collapse=", "), "."
    ))
  }
  # Warn if the treatment/response/event/id/strata columns have reserved names
  # that coincide with a *different* role (e.g. response column named "event").
  roles <- list(treat=treat_col, response=response_col, event=event_col,
                id=id_col)
  for(role in names(roles)){
    col <- roles[[role]]
    if(!is.null(col) && col %in% .TTE.RESERVED.NAMES && col != role){
      stop(paste0("Column '", col, "' supplied as ", role, "_col has a reserved ",
                  "name used internally. Please rename it."))
    }
  }
  for(col in car_strata_cols){
    if(col %in% .TTE.RESERVED.NAMES){
      stop(paste0("car_strata column '", col, "' has a reserved name used ",
                  "internally. Please rename it."))
    }
  }
  invisible(TRUE)
}

.check.variance_type <- function(variance_type){
  OPTIONS <- c(1, 2, 3)
  .check.options("variance_type", variance_type, OPTIONS)
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

.check.estimand_ci <- function(estimand,ci_type){
  OPTIONS.estimand <- c("MH",
                        "ATE")
  .check.options("estimand", estimand, OPTIONS.estimand)
  
  OPTIONS.ci_type <- c("GR", 
                       "Sato", 
                       "mGR")
  .check.options("ci_type", ci_type, OPTIONS.ci_type)
  
  if(estimand=="ATE" & ci_type!="mGR"){
    stop(paste0('Using ', ci_type, ' variance estimator is invalid when estimating ATE.\n Please use modified GR estimator and set ci_type = "mGR".'))
  }
}
