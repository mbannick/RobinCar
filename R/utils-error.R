.return.error <- function(err){
  if(length(err) > 0) stop(paste0(err, sep="\n"))
}

.x.exist.warn <- function() warning("Covariates specified, but no adjustment desired. Ignoring them.")
.z.exist.warn <- function() warning("Strata specified, but no adjustment desired. Ignoring them.")
.z.exist.warn.simple <- function() warning("Strata specified, but simple randomization chosen. Ignoring Z in adjustment.")

.x.miss.warn <- function() warning("Covariates not specified, but adjustment desired. Changing method to ANOVA.")
.z.miss.err <- function() stop("No strata specified, but covariate-adaptive randomization desired. Please provide strata.")

.form.warn <- function() warning("Formula specified; ignoring arguments for covariate variables and overriding adjustment method.")

.car.min.err <- function() stop("Minimization is not compatible with desired adjustment method. Please use ANHECOVA instead.")

.z.include.warn <- function() warning("Must include strata variables Z as covariates, modifying the covariates to include...")

.pu.warn <- function() warning("Prediction unbiasedness does not hold.")
.pu.z.err <- function() stop("Prediction unbiasedness does not hold in all joint levels of Z.")
.pu.z.calibrate <- function() warning("We suggest you try out robin_calibrate.")

.homogeneous.min.error <- function() stop("Inapplicable settings. Try method = 'heterogeneous', include strata as covariates, and a g-family with canonical link.")