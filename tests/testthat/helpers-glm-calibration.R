Robin_calibrate <- function(robin.g_object = G_zx$robin.g_object, # robin.g_object includes robin.g_data, car.scheme and car.z
                            robin.add_x = c("x1"),
                            robin.add_x_to_include_z = TRUE, # if robin.g_object$car.scheme = SR, this arg will be ignored (force to FALSE)
                            robin.vcovHC = c("HC0", "HC1", "HC3")){

  # Robin_calibrate applies ANHECOVA (robin_linear(robin.adj_method = "ANHECOVA")) to robin.g_object with some additional options/arguments
  # customized calibrate on all "g_pred_"s and customized "robin.add_x"
  # default double calibration can use all "g_pred_"s and all "x" and "z" used in robin.g_object to generate "g_pred_"

  robin.c_data <- robin.g_object$robin.g_data %>%
    rename_with(.cols = contains("g_pred_"), .fn = ~ str_replace(.x, "g_pred_", "mu_"))

  robin.c_x <- c(robin.add_x, grep("mu_", x = colnames(robin.c_data), value=TRUE)) # !is.null(robin.add_x) --> double calibration

  if(robin.g_object$car.scheme == "minimization"){
    robin.add_x_to_include_z <- TRUE # give a warning if robin.add_x_to_include_z=F under minimization
  }

  # We can also use robin_linear with ANHECOVA
  robin.c_return <- Robin_g(robin.data = robin.c_data,
                            car.scheme = robin.g_object$car.scheme,
                            car.z = robin.g_object$car.z,
                            robin.x = robin.c_x,
                            robin.formula = NULL,
                            robin.x_to_include_z = robin.add_x_to_include_z, # default
                            robin.adj_method = "heterogeneous",
                            robin.g_family = gaussian(link = "identity"), # ANHECOVA
                            robin.g_accuracy = 7,
                            robin.vcovHC = c("HC0", "HC1", "HC3"))

  return(robin.c_return) # robin.c_return can also be robin objective, with similar structure as robin_g or robin_linear
}
