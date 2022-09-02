#' Robust (potentially stratified) logrank adjustment
#'
#' @param inheritParams robincar_coxscore
#' @param adj_method Adjustment method, one of "CL", "CSL"
#'
#' @import dplyr
#' @import magrittr
#' @export
robincar_logrank <- function(adj_method, ...){

  .check.adj_method.logrank(adj_method)
  result <- robincar_tte(adj_method, ...)

  return(result)

}
