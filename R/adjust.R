#' Generic function for adjustment
#'
#' @param model Object of class Model
adjust <- function (model, ...) {
  UseMethod("adjust", model)
}
