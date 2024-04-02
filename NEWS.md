# RobinCar 0.2.0

## Features

* Added a `NEWS.md` file to track changes to the package.
* Added two new functions: `robincar_linear2` and `robincar_glm2`. These are wrappers for `robincar_linear` and `robincar_glm`, respectively, but they give the user more full control over their covariate adjustment settings. The differences are listed below:
  * The difference between `robincar_linear` and `robincar_linear2` is that in `robincar_linear2`, if you want to include strata variables as covariates in the working model, you need to add those strata variables to the `covariate_cols` as well as the `car_strata_cols`.
  * `robincar_glm2` is exactly `robincar_glm` but only with the `formula` working model functionality.

## Bugfixes

* Fixed two issues with the `formula` argument in `robincar_glm`. The first is if you had overlapping names for your treatment/response/covariate/strata variables (e.g., "A" for treatment col, and "A2" for covariate col), the formula would parse incorrectly. The second is if you used parentheses (e.g., "A*(x1+x2)"), the formula would also parse incorrectly. This would have given breaking errors to the user, and now the behavior is fixed.

# RobinCar 0.1.2

* Updated variance estimator for covariate-adjusted hazard ratio so that it uses the un-adjusted HR estimate rather than the adjusted HR estimate. Both result in consistent variance estimators but using un-adjusted theta is the expected behavior for the user.

# RobinCar 0.1.1

* Support for negative binomial glm (with known or unknown dispersion parameter)
* Beta version of SuperLearner working model
* Updated variance calculation that fixes negative variance when all/most outcome observations are 0

# RobinCar 0.1.0

*Initial Release*
