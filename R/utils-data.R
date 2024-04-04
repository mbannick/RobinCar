
paste.dash <- function(...) paste(..., sep="-")
joint.car_strata <- function(cols) as.factor(do.call(paste.dash, as.list(cols)))
