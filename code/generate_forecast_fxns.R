################################################################################
##  generate_forecast_fxns.R: Script to generate a forecast time series. 
##  The simulated data can then be used to test the forecast partitioning
##  approach. The forecast generation process includes error from 
##  initial conditions uncertainty, parameter error, driver uncertainty, and 
##  process error.
##
##  Author: Andrew Tredennick
##  Date created: March 13, 2018
################################################################################

##  Clear the workspace
rm(list = ls(all.names = TRUE))


####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse)
library(dplyr)



####
####  DEFINE SIMULATOR FUNCTION ------------------------------------------------
####
generate_forecast <- function(z_bar, z_sigma, 
                              theta_bar, theta_sigma, 
                              x_bar, x_sigma, proc_sigma, 
                              ntimes, niters){
  
}