# tmax-pred-data-prep.R: script to extract maximum temperature data across
#                        the prediction grid.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(stars)
library(prism)

# Load formatted BBS data -------------------------------------------------
load("data/pred-coordinates.rda")
# Loads the prediction coordinates (coords.0)

coords.sf <- st_as_sf(data.frame(coords.0),
		      coords = c("X", "Y"),
		      crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")

# Set PRISM directory -----------------------------------------------------
prism_set_dl_dir("~/Dropbox/DSFBGWZ22/data/prism-climate")

# Extract data from 2019 --------------------------------------------------
tmax.vals.list <- list()
years <- 2010:2019
# Create matrix to store all values for 3 forest classes and 1 developed class
J <- nrow(coords.sf)
# Loop through all the sites
vals <- split(1:J, ceiling(seq_along(1:J)/100))
tmax.vals.mat <- matrix(NA, J, length(years))
for (i in 1:length(years)) {
  curr.year <- years[i]
  print(paste("Currently on year ", years[i], sep = ''))
  get_prism_monthlys(type = "tmax", year = curr.year, mon = 5:7, keepZip = FALSE)
  # Get file name to data of interest
  tmax.curr.path.5 <- prism_archive_subset("tmax", "monthly", years = curr.year, mon = 5)
  tmax.curr.path.6 <- prism_archive_subset("tmax", "monthly", years = curr.year, mon = 6)
  tmax.curr.path.7 <- prism_archive_subset("tmax", "monthly", years = curr.year, mon = 7)
  # Get absolute file path
  tmax.curr.abs.5 <- pd_to_file(tmax.curr.path.5)
  tmax.curr.abs.6 <- pd_to_file(tmax.curr.path.6)
  tmax.curr.abs.7 <- pd_to_file(tmax.curr.path.7)
  # Download the raster
  tmax.curr.5 <- read_stars(tmax.curr.abs.5)
  tmax.curr.6 <- read_stars(tmax.curr.abs.6)
  tmax.curr.7 <- read_stars(tmax.curr.abs.7)
  # Take the average value across all three months
  tmax.all <- st_apply(c(tmax.curr.5, tmax.curr.6, tmax.curr.7, along = 3), 
        	       c('x', 'y'), mean)
  if (i == 1) {
    # Get coordinates with buffer for extracting the climate data
    coords.proj <- coords.sf %>%
      st_transform(crs = st_crs(tmax.curr.5))
  }
  for (j in 1:length(vals)) {
    print(j)
    # Get maximum temperature at the start of each route.
    tmp <- st_extract(tmax.all, at = coords.proj[vals[[j]], ])
    # Save the mean values in the list
    tmax.vals.mat[vals[[j]], i] <- tmp[[1]]
  }
}

save(tmax.vals.mat, file = "data/pred-tmax.rda")
