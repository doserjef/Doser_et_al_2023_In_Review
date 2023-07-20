# spOcc-format-data.R: format data in the required data list format for
#                      fitting all models in spOccupancy.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(stars)
library(sf)
library(viridis)

# Load formatted BBS data -------------------------------------------------
load("data/bbs-data-y-det-covs.rda")

# Convert coordinates to Albers Equal Area. These are the centroids
# of each route. 
coords.sf <- st_as_sf(data.frame(coords),
		      coords = c("Longitude", "Latitude"),
		      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Albers equal area across contiguous US.
coords.sf.albers <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")
# Get coordinates in Albers Equal Area (in kilometers)
coords.albers <- st_coordinates(coords.sf.albers)

# Load LULC values --------------------------------------------------------
years <- 2010:2019
J <- nrow(coords.sf)
n.years <- length(years)
grass <- matrix(NA, J, n.years)
for (t in 1:n.years) {
  load(paste("data/lulc-data/lulc-vals-", years[t], sep = ''))
  grass[, t] <- grass.cov
}
# Load Climate data -------------------------------------------------------
load("data/tmax.rda")

# Final processing of climate data ----------------------------------------
# Average values
grass.mean <- apply(grass, 1, mean)
tmax.mean <- apply(tmax.vals.mat, 1, mean)

# Format for spOccupancy --------------------------------------------------
# Format occurrence covariates
occ.covs <- data.frame(grass = grass.mean, 
		       tmax = tmax.mean)

data.list <- list(y = y, 
		  occ.covs = occ.covs, 
		  det.covs = det.covs, 
		  coords = coords.albers)

# Get breeding range data for spOccupancy ---------------------------------
load('data/range.data.rda')
data.list$range.ind <- range.ind
# Only use sites that are within at least one species range
range.indx <- apply(data.list$range.ind, 2, sum)
bad.indx <- which(range.indx == 0)
data.list$y <- data.list$y[, -bad.indx, ]
data.list$occ.covs <- data.list$occ.covs[-bad.indx, ]
data.list$det.covs$day <- data.list$det.covs$day[-bad.indx]
data.list$det.covs$obs <- data.list$det.covs$obs[-bad.indx]
data.list$det.covs$stop <- data.list$det.covs$stop[-bad.indx, ]
data.list$coords <- data.list$coords[-bad.indx, ]
data.list$range.ind <- data.list$range.ind[, -bad.indx]

save(data.list, file = "data/full-data-spOccupancy.rda")
