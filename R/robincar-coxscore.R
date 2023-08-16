#' Robust cox score adjustment
#'
#' @param ... Arguments to robincar_tte, other than `adj_method`
#'
#' @export
robincar_coxscore <- function(...){

  result <- robincar_tte(
    adj_method="coxscore", ...
  )

  return(result)

}
