################################################################################
##  analyze_posteriors.R: R script to plot the posterior distributions of model
##  parameters and to run MCMC diagnostics
##
##  ____________________________________________________________________________
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: October 19, 2016
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
source("./utilities/plotting_theme.R")
# library(devtools) # For installing packages from GitHub
# install_github("atredennick/ecoforecastR") # get latest version
library(ecoforecastR) # MCMC manipulation (by M. Dietze)



####
####  LOAD FITTED MODELS -------------------------------------------------------
####
swe_est_posts <- as.mcmc.list(list(readRDS("../results/ppt_est_posteriors_chain1.RDS"),
                                   readRDS("../results/ppt_est_posteriors_chain2.RDS"),
                                   readRDS("../results/ppt_est_posteriors_chain3.RDS")))
swe_avg_posts <- as.mcmc.list(list(readRDS("../results/ppt_avg_posteriors_chain1.RDS"),
                                   readRDS("../results/ppt_avg_posteriors_chain2.RDS"),
                                   readRDS("../results/ppt_avg_posteriors_chain3.RDS")))



####
#### PLOT POSTERIOR DISTRIBUTIONS ----------------------------------------------
####
## Redefine priors for plotting
r_mu_prior <- log(1.11) # lambda = 1.11 in Hobbs et al. 2015
r_sd_prior <- sd(log(rnorm(100000,1.11,0.024))) # sd_lambda = 0.024 in Hobbs et al. 2015

## Split output for swe_est model
out          <- list(params=NULL, predict=NULL)
mfit         <- as.matrix(swe_est_posts, chains=TRUE)
pred.cols    <- union(grep("z[",colnames(mfit),fixed=TRUE),grep("mu[",colnames(mfit),fixed=TRUE))
chain.col    <- which(colnames(mfit)=="CHAIN")
out$predict  <- mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
out$params   <- mat2mcmc.list(mfit[,-pred.cols])
fitted_model <- out

post_params <- as.data.frame(as.matrix(fitted_model$params))
max_iters   <- nrow(post_params)
post_params <- post_params %>%
  mutate(iteration = 1:max_iters) %>%
  dplyr::select(iteration,b,b1,r,sigma_proc) %>%
  gather(key = parameter, value = estimate, -iteration) %>%
  mutate(prior = c(rnorm(max_iters,0,1000), # b prior
                   rnorm(max_iters,0,1000), # b1 prior
                   rnorm(max_iters,r_mu_prior,r_sd_prior), # r prior
                   runif(max_iters,0,10))) # sd prior

docolor   <- "#278DAF"
prior_col <- "#CF4C26"
docolor   <- "grey90"
prior_col <- "black"
# docolor <- "#307129"
# prior_col <- "#EB6625"
post_params$parameter <- as.factor(post_params$parameter)
levels(post_params$parameter) <- c(expression(beta[0]), expression(beta[1]), "r", expression(sigma[p]))

ggplot(post_params, aes(x = estimate, y = ..density..))+
  geom_histogram(fill = docolor, color = "black", bins = 30)+
  # geom_line(data = filter(post_params, parameter == "r"),
  #           aes(x = prior), 
  #           stat = "density", 
  #           color = "white",
  #           size = 1.2)+
  geom_line(data = filter(post_params, parameter == "r"),
            aes(x = prior), 
            stat = "density", 
            color = prior_col,
            linetype = 2,
            size = 0.8)+
  facet_wrap(~parameter, scales = "free", ncol = 4, labeller = label_parsed)+
  ylab("Posterior density")+
  xlab("Parameter estimate")+
  theme_few()
ggsave(filename = "../figures/bison_post_params.png", 
       height = 3, 
       width = 10, 
       units = "in", 
       dpi = 120)


##  Just the climate effect, for presentations
docolor <- "#307129" #green to match curren beamer style
ggplot(filter(post_params, parameter == "beta[1]"), aes(x = estimate, y = ..density..))+
  geom_histogram(fill = docolor, color = "white", bins = 30)+
  geom_vline(aes(xintercept = 0), color = "grey35", linetype = 2)+
  ylab("Posterior density")+
  xlab("Effect of Jan. Precip.")+
  theme_few()


####
####  MCMC DIAGNOSTICS ---------------------------------------------------------
####
# traceplot(swe_avg_posts)   ## looks good!
# traceplot(swe_est_posts)   ## looks good!
# gelman.diag(swe_est_posts) ## all < 1.1
# heidel.diag(swe_est_posts) ## all passed
# gelman.diag(swe_avg_posts) ## all < 1.1
# heidel.diag(swe_avg_posts) ## all passed


