##  Script to test equivalence of post-hoc simulation from posterior distribution
##  assuming some parameters are set, and refitting the model with constant
##  parameters, set to values from initial fit.

rm(list = ls(all.names = TRUE))



# Load packages -----------------------------------------------------------

library(tidyverse)
library(rjags)

# Simulate data -----------------------------------------------------------

get_ar1 <- function(nsims, z1, x, a, b, c, eps) {
  z <- numeric(nsims+1)
  z[1] <- z1
  for(t in 2:(nsims+1)){
    z[t] <- a + b*z[t-1] + c*x[t] + rnorm(1,0,eps)
  }
  return(z)
}

years <- 30
x <- arima.sim(list("ar",1,0,2), n = years+1) # simulate environment
the_data <- get_ar1(years, 4, x, 10, 0.8, -1, 2)
plot(the_data, type = "l")



