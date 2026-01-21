# RobinCar 1.1.0

## Features

* Included new contrast options for `robincar_glm`: "log_ratio" and "log_odds_ratio". Now, when users ask for "ratio" or "odds_ratio" they get a warning that they should switch to the log scale, as this has better performance. To obtain risk ratios, or odds ratios, exponentiate the estimates.

# RobinCar 1.0.0

## Features

* Including a new function `robincar_mh`, which calculates a Mantel-Haenszel statistic. Please see [Xiaoyu Qiu, Yuhan Qian, Jaehwan Yi, Jinqiu Wang, Yu Du, Yanyao Yi, Ting Ye (2025)](https://arxiv.org/pdf/2408.12541) for methods,
and our new vignette.
* Added two new variance estimators for `robincar_glm`. They are asymptotically equivalent to the original estimator, but have different finite sample properties. See the description of the function for details.

## Bugfixes

* Fixed an issue with printing the results from a negative binomial working model with unknown dispersion parameter in `robincar_glm`.
* Fix when there are ties in time-to-event data. This impacts `robincar_logrank`, `robincar_coxscore`, and `robincar_covhr`.

# RobinCar 0.3.1

## Bugfixes

* When the following things happened simultaneously, the variance results were incorrect: (1) when the system locale was not set to "C", (2) the user did not specify their treatment variable as a factor, and (3) when inconsistent cases were used in naming treatment levels (e.g., "a", "B"). We have fixed this issue.

# RobinCar 0.3.0

## Features

* Changed the name description of RobinCar.

## Breaking Changes

* Simplified user experience for `robincar_glm`: removed `robincar_glm`,
and renamed `robincar_glm2` to `robincar_glm`. Same with `robincar_linear` and `robincar_linear2`. In effect,
the older versions of `robincar_linear` and `robincar_glm` have been deprecated.

## Bugfixes

* Previously, the factor level order passed by the user for the treatment variable was ignored. We have fixed this,
so all estimates will be presented in the order of treatment levels specified by user.

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
