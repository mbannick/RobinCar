# SILENCE R CMD CHECK WORRIES
# regarding no visible binding for global variables

if(getRversion() >= "2.15.1"){
  # This is for using formula notation with glm
  utils::globalVariables(c("."))

  # This is for print.TTEResult and data.table dplyr issues
  utils::globalVariables(c("car_strata", "treat", "N", "name", "strata_col"))
}
