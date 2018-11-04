################################################################################
##  validate_models.R: R script to compare forecasts to observations.
##
##  ____________________________________________________________________________
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: December 4, 2017
################################################################################

##  Clear everything...
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
library(rjags)     # Fitting Bayesian models with JAGS
library(coda)      # MCMC summaries
# library(devtools) # For installing packages from GitHub
# install_github("atredennick/ecoforecastR") # get latest version
library(ecoforecastR) # MCMC manipulation (by M. Dietze)
source("./utilities/plotting_theme.R")



####
####  LOAD DATA AND MCMCs ------------------------------------------------------
####
##  Data
snow_ynp  <- read.csv("../data/west_yellowstone_snotel_summary.csv", row.names = 1) 
bison_raw <- read.csv("../data/YNP_bison_population_size.csv")
bison_dat <- bison_raw %>% 
  dplyr::select(-source) %>%     # drop the source column
  mutate(set = ifelse(year < 2011, "training", "validation")) %>% # make new column for data splits
  left_join(snow_ynp, by="year") # merge in SNOTEL data

##  MCMC
swe_est_posts <- as.mcmc.list(list(readRDS("../results/swe_est_posteriors_chain1.RDS"),
                                   readRDS("../results/swe_est_posteriors_chain2.RDS"),
                                   readRDS("../results/swe_est_posteriors_chain3.RDS")))
swe_avg_posts <- as.mcmc.list(list(readRDS("../results/swe_avg_posteriors_chain1.RDS"),
                                   readRDS("../results/swe_avg_posteriors_chain2.RDS"),
                                   readRDS("../results/swe_avg_posteriors_chain3.RDS")))



####
####  REFORMAT MCMC RESULTS ----------------------------------------------------
####
## Split output
out          <- list(params=NULL, predict=NULL)
mfit         <- as.matrix(swe_est_posts,chains=TRUE)
pred.cols    <- union(grep("z[",colnames(mfit),fixed=TRUE),grep("mu[",colnames(mfit),fixed=TRUE))
chain.col    <- which(colnames(mfit)=="CHAIN")
out$predict  <- mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
out$params   <- mat2mcmc.list(mfit[,-pred.cols])
fitted_model <- out

## Collate predictions
est_swe_preds <- rbind(fitted_model$predict[[1]],
                       fitted_model$predict[[2]],
                       fitted_model$predict[[3]])

## Split output
out          <- list(params=NULL, predict=NULL)
mfit         <- as.matrix(swe_avg_posts,chains=TRUE)
pred.cols    <- union(grep("z[",colnames(mfit),fixed=TRUE),grep("mu[",colnames(mfit),fixed=TRUE))
chain.col    <- which(colnames(mfit)=="CHAIN")
out$predict  <- mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
out$params   <- mat2mcmc.list(mfit[,-pred.cols])
fitted_model <- out

## Collate predictions
avg_swe_preds <- rbind(fitted_model$predict[[1]],
                       fitted_model$predict[[2]],
                       fitted_model$predict[[3]])



####
####  COMPARE PREDICTIONS TO OBSERVATIONS --------------------------------------
####
avg_swe_df <- as.data.frame(avg_swe_preds) 
colnames(avg_swe_df) <- bison_dat$year
avg_swe_df <- avg_swe_df %>%
  mutate(iteration = 1:nrow(avg_swe_df)) %>%
  gather(key = year, value = estimate, -iteration) %>%
  mutate(year = as.numeric(year)) %>%
  inner_join(filter(bison_dat, set == "validation"), by = "year") %>%
  dplyr::select(iteration, year, estimate, count.mean) %>%
  mutate(sq_error = (estimate - count.mean)^2,
         model = "avg_swe") %>%
  group_by(model, year) %>%
  summarise(rmse = sqrt(mean(sq_error)))

est_swe_df <- as.data.frame(est_swe_preds) 
colnames(est_swe_df) <- bison_dat$year
est_swe_df <- est_swe_df %>%
  mutate(iteration = 1:nrow(est_swe_df)) %>%
  gather(key = year, value = estimate, -iteration) %>%
  mutate(year = as.numeric(year)) %>%
  inner_join(filter(bison_dat, set == "validation"), by = "year") %>%
  dplyr::select(iteration, year, estimate, count.mean) %>%
  mutate(sq_error = (estimate - count.mean)^2,
         model = "est_swe") %>%
  group_by(model, year) %>%
  summarise(rmse = sqrt(mean(sq_error)))

rmse_df <- rbind(est_swe_df, avg_swe_df)



####
####  MAKE THE RMSE PLOT AND SAVE ----------------------------------------------
####
mycol <- c("#278DAF", "#CF4C26")
ggplot(rmse_df, aes(x = year, y = rmse, color = model))+
  geom_line()+
  geom_point(size = 3, color = "#EFEFEF")+
  geom_point(size = 2)+
  scale_x_continuous(breaks = rmse_df$year)+
  scale_color_manual(values = mycol, name = NULL, labels = c("Mean SWE","Known SWE"))+
  xlab("Forecast Year")+
  ylab("RMSE (number of bison)")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color="white"),
        panel.background   = element_rect(fill = "#EFEFEF"),
        axis.text          = element_text(size=10, color="grey35", family = "Arial Narrow"),
        axis.title         = element_text(size=12, family = "Arial Narrow", face = "bold"),
        panel.border       = element_blank(),
        axis.line.x        = element_line(color="black"),
        axis.line.y        = element_line(color="black"),
        strip.background   = element_blank(),
        strip.text         = element_text(size=10, color="grey15", family = "Arial Narrow"),
        legend.position    = c(0.2,0.8),
        legend.key         = element_rect(fill = NA),
        legend.background  = element_rect(fill = NA),
        legend.text        = element_text(size=10, color="grey15", family = "Arial Narrow"))
ggsave(filename = "../figures/rmse_by_mod.png", width = 4, height = 4, units = "in", dpi = 120)
