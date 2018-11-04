################################################################################
##  format_weather_projections: R script to collate, reformat, and summarize
##  CMIP5 weather projections by model, run, and scenario.
##
##  ----------------------------------------------------------------------------
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: January 2, 2018
################################################################################

##  Clear the workspace
rm(list = ls(all.names = TRUE))

##  Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # only for RStudio



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse)
library(dplyr)
library(stringr)



####
####  READ IN DATA, ADD ROW COLUMN INFORMATION, SAVE ---------------------------
####
col_names <- c("year",
               "month",
               "day",
               as.character(as.data.frame(read.table("../data/CMIP_YNP/bcca5/COLS_pr.txt"))[,1])
               )

read_csv("../data/CMIP_YNP/bcca5/pr.csv", col_names = col_names) %>%
  gather(key = model, value = ppt, -year, -month, -day) %>%
  separate(model, into = c("model_name", "rep", "scenario"), sep = "[.]") %>%
  group_by(year, month, model_name, scenario) %>%
  summarise(total_ppt_mm = sum(ppt)) %>%
  dplyr::filter(month == 1) %>%
  dplyr::select(model_name, scenario, year, month, total_ppt_mm) %>%
  dplyr::arrange(model_name, scenario, year, month) %>%
  saveRDS(file = "../data/cmip_ynp_formatted_projections.RDS")



