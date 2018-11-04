################################################################################
##  partition_forecast.R: R script to plot observations and posterior
##  predictions from the bison population model and then partition the forecast
##  variance following Dietze 2017, Ecological Applications
##  http://onlinelibrary.wiley.com/doi/10.1002/eap.1589/full
##
##  ____________________________________________________________________________
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: December, 4 2017
################################################################################

##  Clear everything...
rm(list = ls(all.names = TRUE))

##  Set working directory to source file location
root <- "~/Repos/bison_forecast/"
setwd(paste0(root,"code/"))



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse) # Data science functions
library(dplyr)     # Data wrangling
library(ggthemes)  # Pleasing themes for ggplot2
library(cowplot)   # Combining ggplots
library(rjags)     # Fitting Bayesian models with JAGS
library(coda)      # MCMC summaries
# library(devtools) # For installing packages from GitHub
# install_github("atredennick/ecoforecastR") # get latest version
library(ecoforecastR) # MCMC manipulation (by M. Dietze)
source("./utilities/plotting_theme.R")



####
####  LOAD DATA AND MCMCs ------------------------------------------------------
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

##  Split into training and validation sets
training_dat   <- filter(bison_dat, set == "training")
validation_dat <- filter(bison_dat, set == "validation")

##  Set up SWE knowns (2011-2017), relative to scaling of observations
ppt_mean     <- mean(training_dat$ppt_in)
ppt_sd       <- sd(training_dat$ppt_in)
forecast_ppt <- precip_dat %>%
  filter(year %in% validation_dat$year) %>%
  pull(ppt_in)
scl_fut_ppt  <- (forecast_ppt - ppt_mean) / ppt_sd

##  Load MCMC
ppt_est_posts <- as.mcmc.list(list(readRDS("../results/ppt_est_posteriors_chain1.RDS"),
                                   readRDS("../results/ppt_est_posteriors_chain2.RDS"),
                                   readRDS("../results/ppt_est_posteriors_chain3.RDS")))
ppt_avg_posts <- as.mcmc.list(list(readRDS("../results/ppt_avg_posteriors_chain1.RDS"),
                                   readRDS("../results/ppt_avg_posteriors_chain2.RDS"),
                                   readRDS("../results/ppt_avg_posteriors_chain3.RDS")))



####
####  REFORMAT KNOWN SWE MCMC RESULTS ------------------------------------------
####
## Split output
out          <- list(params=NULL, predict=NULL)
mfit         <- as.matrix(ppt_est_posts,chains=TRUE)
pred.cols    <- union(grep("z[",colnames(mfit),fixed=TRUE),grep("mu[",colnames(mfit),fixed=TRUE))
chain.col    <- which(colnames(mfit)=="CHAIN")
out$predict  <- mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
out$params   <- mat2mcmc.list(mfit[,-pred.cols])
fitted_model <- out

## Collate predictions
predictions        <- rbind(fitted_model$predict[[1]],
                            fitted_model$predict[[2]],
                            fitted_model$predict[[3]])
median_predictions <- apply(predictions, MARGIN = 2, FUN = "median")
upper_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.975)})
lower_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.025)})
prediction_df      <- data.frame(year = bison_dat$year,
                                 set = bison_dat$set,
                                 observation = bison_dat$count.mean,
                                 upper_observation = bison_dat$count.mean+bison_dat$count.sd,
                                 lower_observation = bison_dat$count.mean-bison_dat$count.sd,
                                 median_prediction = median_predictions,
                                 upper_prediction = upper_predictions,
                                 lower_prediction = lower_predictions)



####
####  PLOT DATA AND POSTERIOR PREDICTIONS --------------------------------------
####
# obs_color <- "#CF4C26"
# pred_color  <- "#278DAF"
pred_color <- "black"
obs_color  <- "black"
calibration_plot <- ggplot(prediction_df, aes(x=year))+
  geom_ribbon(aes(ymax=upper_prediction, ymin=lower_prediction),
              fill=pred_color, 
              color=NA, 
              alpha=0.2)+
  geom_line(aes(y=median_prediction), color=pred_color, size = 0.2)+
  geom_errorbar(aes(ymin=lower_observation, ymax=upper_observation), 
                width=0.5, 
                color=obs_color, 
                size=0.2)+
  geom_point(aes(y=observation, shape = set, fill = set), 
             color=obs_color, 
             size=1)+
  geom_vline(aes(xintercept=2010), linetype=2,color="grey55")+
  geom_col(data = bison_dat, aes(x = year, y = wint.removal),
           color = "grey55", 
           fill = "grey55",
           width = 0.3)+
  scale_y_continuous(breaks = seq(0,10000,2000))+
  scale_shape_manual(values = c(19,21), 
                     name = NULL, 
                     labels = c("Training data", "Validation data"))+
  scale_fill_manual(values = c("black","white"), 
                    name = NULL, 
                    labels = c("Training data", "Validation data"))+
  ylab("Number of bison")+
  xlab("Year")+
  theme_few()+
  guides(shape = guide_legend(override.aes = list(size=2, 
                                                  fill = c("black","white"))))+
  theme(legend.position = c(0.21,0.85))



####
####  SET UP GCM PROJECTION MATRIX ---------------------------------------------
####
# Set up column names for GCM projection file
col_names <- c("year",
               "month",
               "day",
               as.character(as.data.frame(read.table("../data/CMIP_YNP/bcca5/COLS_pr.txt"))[,1])
)

# Read in GCM projections and format as matrix
gcm_precip <- read_csv("../data/CMIP_YNP/bcca5/pr.csv", col_names = col_names) %>%
  gather(key = model, value = ppt, -year, -month, -day) %>%
  separate(model, into = c("model_name", "rep", "scenario"), sep = "[.]") %>%
  group_by(year, month, model_name, scenario) %>%
  summarise(total_ppt_mm = sum(ppt),
            total_ppt_in = total_ppt_mm*0.0393701) %>%
  ungroup() %>%
  dplyr::filter(month == 1) %>%
  dplyr::filter(year %in% c(2011,2012,2013,2014,2015,2016,2017)) %>%
  dplyr::select(model_name, scenario, year, month, total_ppt_mm, total_ppt_in) %>%
  dplyr::arrange(model_name, scenario, year, month) %>%
  dplyr::mutate(stdzd_precip = (total_ppt_in-ppt_mean) / ppt_sd) %>%
  dplyr::mutate(model_rcp = paste(model_name, scenario,"::")) %>%
  dplyr::select(model_rcp, year, stdzd_precip) %>%
  spread(model_rcp, stdzd_precip)



####
####  PARTITION FORECAST UNCERTAINTY -------------------------------------------
####
##  Function for the ecological process (Gompertz population growth)
iterate_process <- function(Nnow, xnow, r, b, b1, sd_proc, E) {
  xnow[xnow>5] <- 5
  # Log integration of extractions
  e <- log( abs( 1 - (E / Nnow ) ) )
  
  # Determinstic process; log scale
  mu <- log(Nnow) + e + r + b*(log(Nnow) + e) + b1*xnow
  
  # Stochastic process; log scale
  zlog <- rnorm(length(mu), mu, sd_proc)
  
  # Back-transform to arithmetic scale
  N <- exp(zlog)
}


##  Initial condition uncertainty: make forecasts from all MCMC iterations of
##    the final year, but use mean parameter values and no process error.
forecast_steps <- 7
num_iters      <- 5000
E              <- validation_dat$wint.removal
z              <- sample(predictions[,nrow(training_dat)], num_iters, replace = TRUE)
param_summary  <- summary(fitted_model$params)$quantile
r              <- param_summary[6,3]
b              <- param_summary[1,3]
b1             <- param_summary[2,3]
sd_proc        <- param_summary[7,3]
x              <- scl_fut_ppt
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  z <- iterate_process(Nnow = z, xnow = x[t], r, b, b1, sd_proc = 0, E = E[t])
  forecasts[,t] <- z
}
varI <- apply(forecasts,2,var)


##  Initial conditions and parameter uncertainty
forecast_steps <- 7
num_iters      <- 100000
E              <- validation_dat$wint.removal
z              <- sample(predictions[,nrow(training_dat)], num_iters, replace = TRUE)
params         <- as.matrix(fitted_model$params)
sample_params  <- sample.int(nrow(params), size = num_iters, replace = TRUE)
r              <- params[sample_params,"r"]
b              <- params[sample_params,"b"]
b1             <- params[sample_params,"b1"]
sd_proc        <- param_summary[7,3]
x              <- scl_fut_ppt
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  z <- iterate_process(Nnow = z, xnow = x[t], r, b, b1, sd_proc = 0, E = E[t])
  forecasts[,t] <- z
}
varIP <- apply(forecasts,2,var)


##  Initial conditions, parameter, and driver uncertainty
z              <- sample(predictions[,nrow(bison_dat)], num_iters, replace = TRUE)
r              <- params[sample_params,"r"]
b              <- params[sample_params,"b"]
b1             <- params[sample_params,"b1"]
sd_proc        <- param_summary[7,3]
xsamps         <- sample(x = ncol(gcm_precip[2:ncol(gcm_precip)]), size = num_iters, replace = TRUE)
x              <- as.matrix(gcm_precip[2:ncol(gcm_precip)])
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  z <- iterate_process(Nnow = z, xnow = as.numeric(x[t,xsamps]), r, b, b1, sd_proc = 0, E = E[t])
  forecasts[,t] <- z
}
varIPD <- apply(forecasts,2,var)

##  Initial conditions, parameter, driver, and process uncertainty
z              <- sample(predictions[,nrow(bison_dat)], num_iters, replace = TRUE)
r              <- params[sample_params,"r"]
b              <- params[sample_params,"b"]
b1             <- params[sample_params,"b1"]
sd_proc        <- param_summary[7,3]
x              <- as.matrix(gcm_precip[2:ncol(gcm_precip)])
forecasts      <- matrix(data = NA, nrow = num_iters, ncol = forecast_steps)

for(t in 1:forecast_steps){
  z <- iterate_process(Nnow = z, xnow = as.numeric(x[t,xsamps]), r, b, b1, sd_proc = sd_proc, E = E[t])
  forecasts[,t] <- z
}
varIPDE <- apply(forecasts,2,var)


V.pred.sim     <- rbind(varIPDE,varIPD,varIP,varI)
V.pred.sim.rel <- apply(V.pred.sim,2,function(x) {x/max(x)})



####
####  PLOT THE FORECASTING UNCERTAINTY PARTITION -------------------------------
####
var_rel_preds <- as.data.frame(t(V.pred.sim.rel*100))
var_rel_preds$x <- 1:nrow(var_rel_preds)
my_cols <- c("black", "grey55", "grey70","grey90")
variance_plot <- ggplot(data=var_rel_preds, aes(x=x))+
  geom_ribbon(aes(ymin=0, ymax=varIPDE), fill=my_cols[4])+
  geom_ribbon(aes(ymin=0, ymax=varIPD), fill=my_cols[3])+
  geom_ribbon(aes(ymin=0, ymax=varIP), fill=my_cols[2])+
  geom_ribbon(aes(ymin=0, ymax=varI), fill=my_cols[1])+
  ylab("Percent of uncertainty")+
  xlab("Forecast steps")+
  scale_x_continuous(breaks=seq(1,forecast_steps,by=1), 
                     labels=paste(seq(1,forecast_steps,by=1), "yrs"))+
  scale_y_continuous(labels=paste0(seq(0,100,25),"%"))+
  theme_few()

##  For presentations

my_cols <- c("#EB6C26", "#D5DD7F", "#408CAD","grey")
tmpvar <- var_rel_preds
colnames(tmpvar) <- c("AvarIPDE", "BvarIPD", "CvarIP", "DvarI", "x")
var2 <- tmpvar %>%
  gather(simtype, variance, -x)

ggplot(var2, aes(x=x, fill = simtype))+
  geom_ribbon(aes(ymin=0, ymax=variance), color = "black")+
  ylab("Percentage of total variance (%)")+
  xlab("Forecast steps")+
  scale_fill_manual(values = my_cols, name = NULL, 
                    labels = c("Process error", "Driver uncertainty", "Parameter uncertainty", "Initial conditions"))+
  scale_x_continuous(breaks=seq(2,forecast_steps,by=2), 
                     labels=paste(seq(2,forecast_steps,by=2), "yrs"),
                     expand = c(0, 0))+
  scale_y_continuous(labels=paste0(seq(0,100,25),"%"),
                     expand = c(0, 0))+
  theme_few()



####
####  COMBINE PLOTS AND SAVE ---------------------------------------------------
####
plot_grid(calibration_plot, variance_plot, nrow = 2, labels = "AUTO")
ggsave(filename = "../figures/bison_combined.png",
       width = 4,
       height = 6,
       units = "in",
       dpi =200)



####
####  PLOT THE WEATHER PROJECTIONS ---------------------------------------------
####
clim_proj <- as.data.frame(x) %>%
  mutate(year = 1:7) %>%
  gather(key = model, value = precip, -year)

library(viridis)
ggplot(filter(clim_proj, precip < 2.6), aes(x = year+2010, y = precip, color = model))+
  geom_line()+
  scale_color_viridis(discrete = TRUE)+
  guides(color = FALSE)+
  xlab("Year")+
  ylab("Standardized Precipitation")+
  theme_few()
