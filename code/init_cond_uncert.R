################################################################################
##  init_cond_uncert.R: Script to generate figures that demonstrate how
##  uncertainty can be propogated, or not, in forecasts using MCMC output.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
################################################################################

rm(list = ls(all.names = TRUE))


####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse)
library(ggthemes)



####
####  GENERATE FAKE DATA AND MAKE BOXPLOTS -------------------------------------
####
zt <- data.frame(value = rnorm(1000,1000,200),
                 iteration = 1:1000,
                 set = "A")

ggplot(zt, aes(x = set, y=value))+
  geom_jitter(width = 0.15, size = 0.2)+
  geom_boxplot(width = 0.4, alpha = 0.7, outlier.color = "white", outlier.alpha = 0, size = 0.8)+
  geom_point(data = data.frame(set = "A", mu = median(zt$value)), aes(y = mu), size = 3)+
  ylab(expression(z[(t)]))+
  xlab(NULL)+
  theme_classic()+
  theme(axis.text = element_blank())
ggsave(filename = "../figures/example_boxplot.pdf", width = 2, height = 4, units = "in", scale = 0.5)
