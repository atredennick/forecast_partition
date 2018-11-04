################################################################################
##  predict_swe.R: R script to predict SWE at West Yellowstone SNOTEL using the
##  algorithm from Tercek et al. 2016, PLoS One.
##
##  ____________________________________________________________________________
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##                Mike Tercek
##  Date created: December 1, 2017
################################################################################

##  Clear everything...
rm(list = ls(all.names = TRUE))

##  Set working directory to source file location
root      <- "~/Repos/bison_forecast/"
cmip_dir  <- paste0(root, "data/tercek_data/bias_corrected_daily_t_and_p/")
calib_dir <- paste0(root, "data/tercek_data/")
out_dir   <- paste0(root, "data/swe_forecasts/")
setwd(paste0(root,"code/"))



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse)



####
####  DEFINE FUNCTIONS FOR PREDICTING SWE --------------------------------------
####
sim_snow_for_one_water_year <- function(tavg, precip, melt_factor = 0.35, melt_thresh_temperature = 0.5){
  todays_snow <- 0 # start on October 1
  snow_vector <- numeric(length(tavg))
  for(i in 1:length(tavg)){
    if(is.na(tavg[i]) == FALSE){
      todays_snow <- melt_one_day(todays_snow,tavg[i],melt_factor, melt_thresh_temperature)
      todays_snow <- accum_one_day(todays_snow,tavg[i],precip[i],melt_factor, melt_thresh_temperature)
      snow_vector[i] <- todays_snow 
    }
  }
  return(snow_vector)
}

melt_one_day <- function(start_swe, tavg, melt_factor, melt_thresh_temperature){
  if(tavg <= melt_thresh_temperature){
    return(start_swe)
  }else{
    swe_delta <-  (tavg - melt_thresh_temperature) * melt_factor
    end_swe <- max(0, start_swe - swe_delta)
    return(end_swe)
  }
}

accum_one_day <- function(start_swe, tavg, precip, precip_fraction, melt_thresh_temperature){
  if(tavg <= melt_thresh_temperature){
    rain <- 0.0
  }else if(tavg > (melt_thresh_temperature + 6)){
    rain <- 1.0
  }else{
    rain <- precip_fraction * (tavg - melt_thresh_temperature)
  }
  
  snow_increment <- max(0, (1.0 - rain) * precip)
  end_swe <- start_swe + snow_increment
  return(end_swe)
}



####
####  READ IN CMIP5 WEATHER PROJECTIONS AND PREDICT SWE ------------------------
####
weather_files <- as.data.frame(list.files(cmip_dir)) %>%
  dplyr::rename(fname = `list.files(cmip_dir)`) %>%
  dplyr::mutate(fnamesep = fname) %>%
  separate(fnamesep, into = c("source", "variable", "rcp", "model", "biassep"), sep = "_") %>%
  separate(biassep, into = c("bias", "filetype"), sep = "[.]")

swe_preds <- {}
for(do_model in unique(weather_files$model)){
  mod_files <- dplyr::filter(weather_files, model == do_model)
  for(do_rcp in unique(mod_files$rcp)){
    rcp_files <- dplyr::filter(mod_files, rcp == do_rcp)
    prfile <- rcp_files[which(rcp_files$variable == "pr"), "fname"]
    tmaxfile <- rcp_files[which(rcp_files$variable == "tasmax"), "fname"]
    tminfile <- rcp_files[which(rcp_files$variable == "tasmin"), "fname"]
    
    tmp_weather <- read_csv(paste0(cmip_dir,prfile))[,c("dates","West Yellowstone")] %>%
      dplyr::rename(pr = `West Yellowstone`) %>%
      cbind(read_csv(paste0(cmip_dir,tmaxfile))[,"West Yellowstone"]) %>%
      dplyr::rename(tmax = `West Yellowstone`) %>%
      cbind(read_csv(paste0(cmip_dir,tminfile))[,"West Yellowstone"]) %>%
      dplyr::rename(tmin = `West Yellowstone`) %>%
      dplyr::mutate(date = dates) %>%
      separate(dates, into = c("month", "day", "year"), sep = "/") %>%
      dplyr::filter(year %in% seq(2011,2017,by = 1)) %>%
      dplyr::mutate(tavg = (tmin+tmax)/2,
                    swe = sim_snow_for_one_water_year(tavg, pr)) %>%
      group_by(year) %>%
      summarise(accum_swe = sum(swe)) %>%
      dplyr::mutate(rcp = do_rcp, model = do_model)
    
    swe_preds <- rbind(swe_preds, tmp_weather)
  }
}



####
####  SAVE THE PREDICTIONS -----------------------------------------------------
####
dir.create(out_dir, showWarnings = FALSE)
saveRDS(object = swe_preds, file = paste0(out_dir, "swe_predictions.RDS"))








