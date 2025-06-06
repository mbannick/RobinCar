#' Robust (potentially stratified) logrank adjustment
#'
#' Perform a robust covariate-adjusted logrank test ("CL") that can be stratified ("CSL") if desired.
#'
#' Note: Since RobinCar version 0.4.0, the variance of the test statistic has changed
#' to better accommodate tied event times.
#'
#' \eqn{\hat{\sigma}_{\rm CL}^2} and \eqn{\hat{\sigma}_{\rm CSL}^2} for the
#' covariate-adjusted stratified log-rank test are given by Ye, Shao, and Yi (2024) after
#' equation (8) on page 700, with \eqn{\hat{\sigma}_{\rm SL}^2} replaced
#' by the following estimator, which is the standard denominator of the logrank test:
#'
#' \deqn{\frac1n \sum_{j} \sum_{i} v_j(t_i)}
#' \deqn{v_j(t_i) = \frac{Y_{0,j}(t)Y_{1,j}(t)d_j(t)\left[Y_j(t) - d_j(t)\right]}{Y_j(t)^2\left[Y_j(t) - 1 \right]}}
#' where \eqn{t_i} are strata-specific unique failure times, \eqn{d_j(t)} is the number of events at time \eqn{t} in strata \eqn{j},
#' \eqn{Y_j(t)} is the number at risk within strata \eqn{j} at time \eqn{t},
#' and \eqn{Y_{a,j}(t)} is the number at risk within strata \eqn{j} and treatment \eqn{a} at time \eqn{t}.
#'
#' Please see Ye, Shao, and Yi (2024)'s "Covariate-adjusted log-rank test: guaranteed efficiency
#' gain and universal applicability" in \emph{Biometrika} for more details about \eqn{\hat{\sigma}_{\rm CSL}^2}.
#'
#' @param adj_method Adjustment method, one of "CL", "CSL"
#' @param ... Additional arguments to `robincar_tte`
#'
#' @export
#'
#' @examples
#' library(magrittr)
#' library(dplyr)
#' library(forcats)
#' set.seed(0)
#' n=100
#' data.simu0=data_gen(n=n,
#'                     theta=0,
#'                     randomization="permuted_block",
#'                     p_trt=0.5,
#'                     case="case2") %>% mutate(strata1=sample(letters[1:3],n,replace=TRUE),
#'                                              strata2=sample(LETTERS[4:5],n,replace=TRUE))
#'
#' out <- robincar_logrank(df=data.simu0,
#'                         treat_col="I1",
#'                         p_trt=0.5,
#'                         ref_arm=0,
#'                         response_col="t",
#'                         event_col="delta",
#'                         covariate_cols=c("model_z1", "model_z2"),
#'                         car_scheme="simple",
#'                         adj_method=c("CL"))
#'
#' set.seed(0)
#' n=100
#' data.simu0=data_gen(n=n,
#'                     theta=0,
#'                     randomization="permuted_block",
#'                     p_trt=0.5,
#'                     case="case1")
#'
#' data.simu <- data.simu0 %>%
#'   tidyr::pivot_longer(cols=starts_with("car_strata"),
#'                       names_prefix="car_strata",
#'                       names_to="strt") %>%
#'   filter(value==1) %>% select(-value) %>%
#'   mutate(strt=forcats::as_factor(strt)) %>%
#'   select(t,strt) %>%
#'   left_join(data.simu0, .)
#'
#' out1 <- robincar_logrank(df=data.simu,
#'                          treat_col="I1",
#'                          p_trt=0.5,
#'                          ref_arm=0,
#'                          response_col="t",
#'                          event_col="delta",
#'                          car_strata_cols="strt",
#'                          covariate_cols=NULL,
#'                          car_scheme=c("permuted-block"),
#'                          adj_method=c("CSL")
#' )
#'
#' @returns A result object with the following attributes:
#'
#' \item{result}{A list: "statistic" is the adjusted logrank test statistic which can be used to obtain p-values; "U" and "se" are the numerator and denominator of the test statistic, respectively.}
#' \item{settings}{The covariate adjustment settings used.}
#' \item{original_df}{The dataset supplied by the user.}
#'
robincar_logrank <- function(adj_method, ...){

  .check.adj_method.logrank(adj_method)
  result <- robincar_tte(adj_method, ...)

  return(result)

}
