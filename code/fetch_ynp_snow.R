################################################################################
##  fetch_ynp_snow.R: script to download and summarize SNOTEL data for
##  Yellowstone National Park.
##
## _____________________________________________________________________________
##  Author: Andrew Tredennick
##  Date created: September 5, 2017
################################################################################



rm(list = ls(all.names = T))

##  Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # only for RStudio


####
####  LOAD LIBRARIES ----
####
# library(devtools)
# install_github("khufkens/snotelr")
library(tidyverse)
library(dplyr)
library(stringr)
library(snotelr)
library(wux)



####
####  FETCH THE WEST YELLOWSTONE SNOTEL DATA (SITE #924) ----
####

### WEST YELLOWSTONE SNOTEL
snotel.info(path = "../data/.") 
download.snotel(site = 924) 
file.remove("snotel_metadata.csv")

##  Collate the data into water years (October - April)
water_months <- as.character(c(10,11,12,1,2,3,4))
addyear_months <- as.character(c(10,11,12))

ynp_snotel <- read.csv("../data/snotel_924.csv", skip = 58) %>%
  dplyr::rename(date = Date,
                snow_water_eq_mm    = Snow.Water.Equivalent..mm..Start.of.Day.Values,
                precip_accum_mm     = Precipitation.Accumulation..mm..Start.of.Day.Values,
                max_air_temp_degC   = Air.Temperature.Maximum..degC.,
                min_air_temp_degC   = Air.Temperature.Minimum..degC.,
                avg_air_temp_degC   = Air.Temperature.Average..degC.,
                precip_increment_mm = Precipitation.Increment..mm.) %>%
  separate(date, into = c("year", "month", "day"), sep = "-") %>%
  filter(month %in% water_months) %>%
  mutate(year = as.numeric(year),
         newyear = ifelse(month %in% addyear_months, year+1, year)) %>%
  group_by(newyear) %>%
  summarise(mean_snow_water_equiv_mm  = mean(snow_water_eq_mm, na.rm=TRUE),
            accum_snow_water_equiv_mm = sum(snow_water_eq_mm, na.rm=TRUE),
            max_snow_water_equiv_mm   = max(snow_water_eq_mm, na.rm=TRUE),
            sd_snow_water_equiv_mm    = sd(snow_water_eq_mm, na.rm=TRUE)) %>%
    dplyr::rename(year = newyear) %>%
    arrange(year)

write.csv(ynp_snotel, "../data/west_yellowstone_snotel_summary.csv")



####
####  DOWNLOAD AND PROCESS SNOW DEPTH DATA FROM GCMs ----
####
# dir.create("../data/CMIP5/", showWarnings = FALSE)
# CMIP5fromESGF(save.to = "/Users/atredenn/repos/ecocast_compare/data/CMIP5/",
#               models = c("CanESM2"),
#               variables = c("tas"),
#               experiments= c("historical"))


