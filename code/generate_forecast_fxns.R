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
library(ggthemes)



####
####  DEFINE SIMULATOR FUNCTION ------------------------------------------------
####
generate_forecast <- function(z_bar, 
                              z_sigma, 
                              a_bar,
                              a_sigma,
                              theta_bar = NULL, 
                              theta_sigma = NULL, 
                              x_bar = NULL, 
                              x_sigma = NULL, 
                              proc_sigma, 
                              n_times, 
                              n_iters,
                              seed = NULL){
  if(is.null(seed) == FALSE){
    set.seed(seed)
  }
  
  if(is.null(x_bar) == FALSE){
    if(is.matrix(x_bar) == FALSE){
      stop("x_bar must be a matrix even if only providing one covariate.
           Coerce to a single column matrix, with the number of rows equal
           to the number of time steps to be forecast.")
    }
    
    if((ncol(x_bar)-1) != length(theta_bar)){
      stop("The number of covariates [ncol(x_bar)] and the number of
           coefficients [length(theta_bar)] do not match.")
    }
  }
 
  
  # Generate vector of initial conditions
  Z <- rnorm(n_iters, z_bar, z_sigma)
  
  # Generate vector of autoregressive terms
  A <- rnorm(n_iters, a_bar, a_sigma)
  
  # Generate array of time-varying covariates
  if(is.null(x_bar) == FALSE){
    X <- array(NA, dim = c((ncol(x_bar)+1), n_iters, n_times))
    X[1,,] <- 1 # make first rows 1 for the intercept
    for(dot in 1:n_times){
      for(dox in 1:ncol(x_bar)){
        X[dox+1,,dot] <- rnorm(n_iters, x_bar[dot,dox], x_sigma[dot,dox])
      }
    }
  }
  
  
  # Generate matrix of time-invariant parameters
  if(is.null(theta_bar) == FALSE){
    Theta <- matrix(NA, length(theta_bar), n_iters)
    for(dotheta in 1:length(theta_bar)){
      Theta[, dotheta] <- rnorm(n_iters, theta_bar[dotheta], theta_sigma[dotheta])
    }
  }
  
  # Generate matrix of forecasts
  forecasts <- matrix(NA, n_times+1, n_iters)
  forecasts[1,] <- Z
  
  if(is.null(theta_bar)){
    for(iiter in 1:n_times){
      Z <- A*Z
      forecasts[iiter+1, ] <- rnorm(n_iters, Z, proc_sigma) 
    }
  }else{
    for(iiter in 1:n_times){
      mu <- A*Z[iiter] + X[iiter, ]*Theta
      forecasts[iiter, ] <- rnorm(n_iters, mu, proc_sigma) 
    }
  }
  
  return(forecasts)
}



####
####  GENERATE FORECASTS -------------------------------------------------------
####
z_bar <- -20
z_sigma <- 1.5
a_bar <- 0.8
a_sigma <- 0.1
proc_sigma <- 1.6
n_times <- 10
n_iters <- 100
seed <- 126

outcasts_all <- generate_forecast(z_bar,
                                  z_sigma, 
                                  a_bar,
                                  a_sigma,
                                  theta_bar = NULL, 
                                  theta_sigma = NULL, 
                                  x_bar = NULL, 
                                  x_sigma = NULL, 
                                  proc_sigma, 
                                  n_times, 
                                  n_iters,
                                  seed)

outcasts_noproc <- generate_forecast(z_bar,
                                  z_sigma, 
                                  a_bar,
                                  a_sigma,
                                  theta_bar = NULL, 
                                  theta_sigma = NULL, 
                                  x_bar = NULL, 
                                  x_sigma = NULL, 
                                  proc_sigma = 0, 
                                  n_times, 
                                  n_iters,
                                  seed)

outcasts_noparam <- generate_forecast(z_bar,
                                  z_sigma, 
                                  a_bar,
                                  a_sigma = 0,
                                  theta_bar = NULL, 
                                  theta_sigma = NULL, 
                                  x_bar = NULL, 
                                  x_sigma = NULL, 
                                  proc_sigma, 
                                  n_times, 
                                  n_iters,
                                  seed)

outcasts_noinit <- generate_forecast(z_bar,
                                      z_sigma = 0, 
                                      a_bar,
                                      a_sigma,
                                      theta_bar = NULL, 
                                      theta_sigma = NULL, 
                                      x_bar = NULL, 
                                      x_sigma = NULL, 
                                      proc_sigma, 
                                      n_times, 
                                      n_iters,
                                      seed)

outcasts_combined <- as.data.frame(outcasts_all) %>%
  mutate(year = 0:n_times) %>%
  gather(key,value,-year) %>%
  mutate(simulation = "A) All sources") %>%
  rbind(as.data.frame(outcasts_noproc) %>%
          mutate(year = 0:n_times) %>%
          gather(key,value,-year) %>%
          mutate(simulation = "D) No process error")) %>%
  rbind(as.data.frame(outcasts_noinit) %>%
          mutate(year = 0:n_times) %>%
          gather(key,value,-year) %>%
          mutate(simulation = "B) No initial conditions error")) %>%
  rbind(as.data.frame(outcasts_noparam) %>%
          mutate(year = 0:n_times) %>%
          gather(key,value,-year) %>%
          mutate(simulation = "C) No parameter error"))

ggplot(outcasts_combined, aes(x = year, y = value, group = key))+
  geom_line(color = "darkgrey", alpha = 0.8)+
  facet_wrap(~simulation, ncol = 4)+
  scale_x_continuous(breaks = c(0,2,4,6,8,10))+
  ylab("State of interest")+
  xlab("Forecast year")+
  theme_few()+
  theme(axis.text.y = element_blank())
# ggsave(filename = "../figures/forecast_uncertainty_example.pdf", width = 8.5, height = 3, units = "in")



####
####  TEST CASE OF ANOVA WITH 2 SOURCES OF ERROR ONLY --------------------------
####
z_bar <- -20
z_sigma <- 1.5
a_bar <- 0.8
a_sigma <- 0.1
proc_sigma <- 0
n_times <- 10
n_iters <- 100
seed <- 126

outcasts_all <- generate_forecast(z_bar,
                                  z_sigma, 
                                  a_bar,
                                  a_sigma,
                                  theta_bar = NULL, 
                                  theta_sigma = NULL, 
                                  x_bar = NULL, 
                                  x_sigma = NULL, 
                                  proc_sigma, 
                                  n_times, 
                                  n_iters,
                                  seed)

outcasts_noinit <- generate_forecast(z_bar,
                                     z_sigma = 0, 
                                     a_bar,
                                     a_sigma,
                                     theta_bar = NULL, 
                                     theta_sigma = NULL, 
                                     x_bar = NULL, 
                                     x_sigma = NULL, 
                                     proc_sigma, 
                                     n_times, 
                                     n_iters,
                                     seed)

outcasts_noparam <- generate_forecast(z_bar,
                                      z_sigma, 
                                      a_bar,
                                      a_sigma = 0,
                                      theta_bar = NULL, 
                                      theta_sigma = NULL, 
                                      x_bar = NULL, 
                                      x_sigma = NULL, 
                                      proc_sigma, 
                                      n_times, 
                                      n_iters,
                                      seed)

outcasts_combined <- as.data.frame(outcasts_all) %>%
  mutate(year = 0:n_times) %>%
  gather(key,value,-year) %>%
  mutate(simulation = "all") %>%
  rbind(as.data.frame(outcasts_noinit) %>%
          mutate(year = 0:n_times) %>%
          gather(key,value,-year) %>%
          mutate(simulation = "noinit")) %>%
  rbind(as.data.frame(outcasts_noparam) %>%
          mutate(year = 0:n_times) %>%
          gather(key,value,-year) %>%
          mutate(simulation = "noparam"))

ggplot(outcasts_combined, aes(x = year, y = value, group = key))+
  geom_line(color = "darkgrey", alpha = 0.8)+
  facet_wrap(~simulation, ncol = 4)+
  scale_x_continuous(breaks = c(0,2,4,6,8,10))+
  ylab("State of interest")+
  xlab("Forecast year")+
  theme_few()+
  theme(axis.text.y = element_blank())

##  Calculate variance
forecast_var <- outcasts_combined %>%
  group_by(year, simulation) %>%
  summarise(variance = var(value)) %>%
  spread(simulation, variance) %>%
  mutate(interaction = noinit + noparam) %>%
  gather(simulation, variance, -year)

ggplot(forecast_var, aes(x = year, y = variance))+
  geom_ribbon(aes(ymax = forecast_var[forecast_var$simulation=="all","variance"],
                  ymin = forecast_var[forecast_var$simulation=="interaction","variance"]),
              fill = "grey88")+
  geom_line(aes(color = simulation), size = 1.3)+
  annotate("text", x = 8, y = 30, label = "Interaction effect")+
  xlab("Forecast year")+
  ylab("Forecast variance")+
  scale_color_colorblind(name = NULL, labels = c("All sources","IC + Param. Error","Param Error","IC Error"))+
  scale_x_continuous(breaks = c(0,2,4,6,8,10))+
  theme_few()+
  theme(legend.position = c(0,1),
        legend.justification=c(-0.1, 1.1))
ggsave(filename = "../figures/example_interaction_effect.pdf", width = 4.5, height = 4, units = "in")

interaction_var <- forecast_var %>%
  spread(simulation, variance) %>%
  mutate(interaction_effect = all - interaction)

ggplot(interaction_var, aes(x = year, y = interaction_effect))+
  geom_line()
  
