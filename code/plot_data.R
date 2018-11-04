################################################################################
##  plot_data.R: R script to plot the time series of YNP bison counts and
##  the West Yellowstone SNOTEL soil water equivalent (annual mean).
##
##  ____________________________________________________________________________
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: December 1, 2017
################################################################################

##  Clear everything
rm(list = ls(all.names = TRUE))

##  Set working directory to source file location
root <- "~/Repos/bison_forecast/"
setwd(paste0(root,"code/"))



####
####  LOAD LIBRARIES ----
####
library(tidyverse) # Data science functions
library(dplyr)     # Data wrangling
library(ggthemes)  # Pleasing themes for ggplot2
library(cowplot)   # Combining ggplots
source("./utilities/plotting_theme.R") # Source my plotting theme



####
####  LOAD DATA ----------------------------------------------------------------
####
snow_ynp  <- read.csv("../data/west_yellowstone_snotel_summary.csv", row.names = 1) 
weather_dat <- read.csv("../data/PRISM_ppt_tmin_tmean_tmax_tdmean_vpdmin_vpdmax_provisional_4km_197001_201711_44.8090_-110.5728.csv", skip = 10)
bison_raw <- read.csv("../data/YNP_bison_population_size.csv")

##  Reformat weather data from PRISM
weather_dat <- weather_dat %>%
  dplyr::select(-tdmean..degrees.F., -vpdmin..hPa., -vpdmax..hPa.) %>%
  dplyr::rename(date = Date,
                ppt_in = ppt..inches.,
                tmin_F = tmin..degrees.F.,
                tmean_F = tmean..degrees.F.,
                tmax_F = tmax..degrees.F.) %>%
  separate(date, into = c("year", "month"), sep = "-")

precip_dat <- weather_dat %>%
  dplyr::select(year, month, ppt_in) %>%
  filter(month %in% c("01")) %>%
  mutate(year = as.integer(year))

##  Reformat bison data and combine with weather data
bison_dat <- bison_raw %>% 
  dplyr::select(-ends_with("source")) %>%     # drop the source column
  mutate(set = ifelse(year < 2011, "training", "validation")) %>% # make new column for data splits
  left_join(snow_ynp, by="year") %>% # merge in SNOTEL data
  left_join(precip_dat, by="year")



####
####  PLOT BISON AND SNOW DATA -------------------------------------------------
####
plot_data <- bison_dat %>%
  dplyr::select(year, set, count.mean, count.sd, wint.removal, ppt_in)

docolor  <- "#278DAF"
altcolor <- "#CF4C26"
docolor  <- "black"
altcolor <- "black"
bison_plot <- ggplot(plot_data, aes(x = year, y = count.mean, color = set))+
  geom_line(alpha = 0.6)+
  geom_point(size=1.5)+
  geom_errorbar(aes(ymin = count.mean-count.sd, ymax = count.mean+count.sd), width=0.5, size=0.5)+
  geom_col(aes(y = wint.removal), color = "grey55", fill = "grey55", width = 0.3)+
  scale_color_manual(values = c(docolor, altcolor))+
  scale_y_continuous(breaks = seq(0,6000,1000))+
  ylab("Number of bison")+
  xlab("Year")+
  theme_few()+
  guides(color = FALSE)

snow_plot <- ggplot(plot_data, aes(x = year, y = ppt_in, color = set))+
  geom_line(alpha = 0.6)+
  geom_point(size=1.5)+
  scale_color_manual(values = c(docolor, altcolor))+
  ylab("January Precipitation (in)")+
  xlab("Year")+
  theme_few()+
  guides(color = FALSE)

the_plots <- list(bison_plot, snow_plot)
suppressWarnings( # Ignore warning about 6 NA rows for errorbars where sd not reported
  plot_grid(plotlist = the_plots, labels = "AUTO", ncol = length(the_plots))
)
ggsave(filename = "../figures/bison_data_plots.png", height = 3, width = 8, units = "in", dpi = 120)



####
####  OLD
####
# bison_growth_data <- bison_dat %>%
#   dplyr::select(year, count.mean) %>%
#   mutate(id = 1) %>% # constant id to work with ave()
#   mutate(growth_rate = ave(count.mean, id, FUN=function(x) c(0, diff(log(x)))))

# bison_growth <- ggplot(bison_growth_data, aes(x = year, y = growth_rate))+
#   geom_line(color = docolor, alpha = 0.6)+
#   geom_point(size=1.5, color = docolor)+
#   ylab("Population growth rate (r)")+
#   xlab("Year")+
#   my_theme
