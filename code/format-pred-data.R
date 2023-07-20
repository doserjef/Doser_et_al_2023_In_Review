# format-pred-data.R: script to format the grassland and maximum 
#                     temperature data for predicting their effects
#                     across the continental US.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(stars)
library(sf)
library(viridis)

# Load prediction grid ----------------------------------------------------
load("data/pred-coordinates.rda")

# Load LULC values --------------------------------------------------------
# TODO: hardcoded
years <- 2010:2019
J <- nrow(coords.0)
n.years <- length(years)
grass <- matrix(NA, J, n.years)
for (t in 1:n.years) {
  load(paste("data/lulc-data/pred-lulc-vals-", years[t], sep = ''))
  grass[, t] <- grass.cov
}
# Load Climate data -------------------------------------------------------
load("data/pred-tmax.rda")

# Final processing of climate data ----------------------------------------
# Average values
grass.0 <- apply(grass, 1, mean)
tmax.0 <- apply(tmax.vals.mat, 1, mean)

# Save to file
save(grass.0, tmax.0, coords.0, file = "data/full-pred-data.rda")
