#' Robust cox score adjustment
#'
#' @param ... Arguments to robincar_tte, other than `adj_method`
#'
#' @export
#'
#' @returns A result object with the following attributes:
#'
#' \item{result}{A list: "statistic" is the robust Cox score test statistic which can be used to obtain p-values; "U" and "se" are the numerator and denominator of the test statistic, respectively.}
#' \item{settings}{The covariate adjustment settings used.}
#' \item{original_df}{The dataset supplied by the user.}
#'
robincar_coxscore <- function(...){

  result <- robincar_tte(
    adj_method="coxscore", ...
  )

  return(result)

}
