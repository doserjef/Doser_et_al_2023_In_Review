# lulc-data-prep.R: script to prepare the USGS EROS LULC data to calculate
#                   grassland area across the BBS routes. Note the raw 
#                   LULC data files are not provided on GitHub.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(stars)
library(sf)

# Get filename to process -------------------------------------------------
# Grab the name from the commandline when running the script. Alternatively, 
# could comment this out and just write the file name manually. 
# Name should be full path name relative to the working directory
file.name <- commandArgs(trailingOnly = TRUE)
if(length(file.name) == 0) base::stop('Need to give the file name to process')

# Load formatted BBS data -------------------------------------------------
load("data/bbs-data-y-det-covs.rda")
# Loads y (detection-nondetection data array), det.covs (detection coordinates), 
# and coords (the spatial coordinates of the beginning of each route). 

coords.sf <- st_as_sf(data.frame(coords),
		      coords = c("Longitude", "Latitude"),
		      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Albers equal area across contiguous US.
coords.sf.albers <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")
# Get coordinates in Albers Equal Area
coords.albers <- st_coordinates(coords.sf.albers)

# Calculating covariates within a 1km buffer radius around start of route 
buffer.radius <- 1000

# Land Cover Classes ------------------
# 1: Water
# 2: Developed
# 3: Mechanically Disturbed National Forests
# 4: Mechanically Disturbed Other Public Lands
# 5: Mechanically Disturbed Private
# 6: Mining
# 7: Barren
# 8: Deciduous Forest
# 9: Evergreen Forest
# 10: Mixed Forest
# 11: Grassland
# 12: Shrubland
# 13: Cropland
# 14: Hay/Pasture Land
# 15: Herbaceous Wetland
# 16: Woody Wetland
# 17: Perennial Ice/Snow

# Load rasters ------------------------------------------------------------
lulc.curr <- read_stars(paste("~/DSFBGWZ22/data/eros-lulc/", file.name, sep = ''))
# Convert coordinates to the coordinates of the raster. 
coords.ACEA <- coords.sf.albers %>%
  st_transform(crs = st_crs(lulc.curr))
# Buffer the coordinates
coords.ACEA.buffered <- st_buffer(coords.ACEA, dist = buffer.radius)
sum.cover <- function(a, val) {
  return(sum(a == val))
}
prop.grass <- function(a) {
  mean(a %in% c(11), na.rm = TRUE)
}
# Create matrix to store all values for grassland variables
J <- nrow(coords.albers)
# Loop through all the sites
vals <- split(1:J, ceiling(seq_along(1:J)/100))
grass.cov <- rep(0, J)
lulc.covs <- matrix(NA, nrow(coords.ACEA), 2)
for (j in 1:length(vals)) {
  print(j)
  coords.curr <- coords.ACEA.buffered[vals[[j]], ]
  tmp <- aggregate(lulc.curr, by = coords.curr, FUN = prop.grass)
  grass.cov[vals[[j]]] <- tmp[[1]]
}

# Save results ------------------------------------------------------------
# Extract year from the file name
curr.year <- parse_number(file.name)
save(grass.cov, file = paste("data/lulc-data/lulc-vals-", curr.year, sep = ''))
