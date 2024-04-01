.return.error <- function(err){
  if(length(err) > 0) stop(paste0(err, sep="\n"))
}

.x.exist.warn <- function() warning("Covariates specified, but no adjustment desired. Ignoring them.")
.z.exist.warn <- function() warning("Strata specified, but no adjustment desired. Ignoring them.")
.z.exist.warn.simple <- function() warning("Strata specified, but simple randomization chosen. Ignoring Z in adjustment.")

.x.miss.err <- function() warning("Must have some covariates specified.")
.x.miss.warn <- function() warning("Covariates not specified, but adjustment desired. Changing method to ANOVA.")
.z.miss.err <- function() stop("No car_strata specified, but covariate-adaptive randomization desired. Please provide car_strata.")

.car.min.err <- function() stop("Minimization is not compatible with desired adjustment method. Please use ANHECOVA instead.")

.z.include.warn <- function() warning("Must include car_strata variables Z as covariates, modifying the covariates to include...")

# ERRORS FOR GLM

.form.warn <- function() warning("Formula specified; ignoring arguments for covariate variables and overriding adjustment method.")

.pu.warn <- function() warning("Prediction unbiasedness does not hold.")
.pu.z.warn <- function() warning("Prediction unbiasedness does not hold in all joint levels of Z.")
.pu.z.calibrate <- function() warning("We suggest you try out robin_calibrate.")

.homogeneous.min.error <- function() stop("Inapplicable settings. Try method = 'heterogeneous', include car_strata as covariates, and a g-family with canonical link.")

.miss.covariate.calibrate <- function() stop("Cannot find add_x covariate in original dataset.")

# ERRORS FOR TTE
.csl.nostrata.warn <- function() warning("If stratified test needs to be performed under simple randomization, please specify car_scheme as permuted-block or biased-coin.")
.cox.carscheme.error <- function() stop("No applicable robust score test for the specified car_scheme, try to use robtest_logrank() function.")
.tte.nomethod.error <- function() stop("No applicable adjustment methods for the specified adj_method.")
.sparse.car_strata.warn <- function() warning("Conditional variance in some car_strata equals 0, removing these car_strata.")

.lin.dep.error <- function() warning("Removing model variables that are linearly dependent with the stratification variables.")
.lin.dep.strat.error <- function() warning("The adjusted covariates are linearly dependent on car_strata, covariate stratified logrank test reduces to the stratified logrank test.")
