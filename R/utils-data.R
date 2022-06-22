
paste.dash <- function(...) paste(..., sep="-")
joint.strata <- function(cols) as.factor(do.call(paste.dash, as.list(cols)))
