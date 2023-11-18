# bbs-data-prep.R: this script takes the raw BBS data and preps it for analysis. Note
#                  that this script does not generate the covariate data. This script
#                  is provided to show how raw BBS data were extracted from data 
#                  downloaded from the USGS BBS website, but the raw BBS data files
#                  to run the script are not available on GitHub
# Author: Jeffrey W. Doser

rm(list = ls())
library(tidyverse)
library(lubridate)
library(sp)
library(raster)
library(sf)
library(stars)
# For linking bird codes with species info.
library(wildlifeR)

# Read in BBS Data --------------------------------------------------------
# This code extracts the BBS data directly from BBS downloaded data. It is 
# memory intensive, and so it is commented out and I read in the saved object
# bbs.dat below that contains the data for the species of interest across the 
# entire continental US from 1972-present. 
# These are the 10-stop summary data
# bbs.dat <- list.files(path = "data/BBS/States/", full.names = TRUE) %>%
#   lapply(read.csv) %>%
#   bind_rows()
# # Get associated route data
 route.dat <- read.csv("data/BBS/routes.csv")
# # Note that route id is nested within state. 
# # Join BBS data with route data
# bbs.dat <- left_join(bbs.dat, route.dat, by = c('Route', 'CountryNum', 'StateNum'))
# # Grab BBS data from 1970 and later, as well as routes with their starting points 
# # east of the 100th meridian. 
# bbs.dat <- bbs.dat %>%
#   filter(Year >= 2010, Year < 2020, RPID == 101)
# # Select columns of interest
# bbs.dat <- bbs.dat %>%
#   dplyr::select(RouteDataID, RouteName, Latitude, Longitude, Year, AOU, starts_with("Count")) %>%
#   # Data set only include US records, so can remove CountryNum
#   dplyr::select(-CountryNum)
# # Fill in implicit zeros. For this situation, gives a value for each species in 
# # each existing combination of RouteDataID, Latitude, Longitude, and Year
# # RouteDataID is a unique identifier for each combo of CountryNum, StateNum, Route, RPID, and Year
# bbs.dat <- bbs.dat %>%
#   complete(AOU, nesting(RouteDataID, RouteName, Latitude, Longitude, Year))
# # Replace NAs with 0s for all columns at once.
# bbs.dat <- bbs.dat %>%
#   replace(is.na(.), 0)
# # Extract community of eastern forest birds.
# aou.info <- AOU_species_codes %>%
#   mutate(AOU = spp.num)
# bbs.dat <- left_join(bbs.dat, aou.info, by = c("AOU"))
# comm.group.dat <- read.csv("data/bird-species-table-bateman.csv")
# my.sp.code <- comm.group.dat %>%
#   filter(Group %in% c("Grasslands"))
# bbs.dat <- bbs.dat %>%
#   mutate(alpha.code = as.character(alpha.code)) %>%	
#   filter(alpha.code %in% my.sp.code$Code)
# # Save the BBS data to avoid having to run the above memory intensive code a bunch.
# save(bbs.dat, file = "data/bbs-raw-data.rda")

# Load in the above BBS data ----------------------------------------------
load("data/bbs-raw-data.rda")
# Only use 2019 data.
bbs.dat <- bbs.dat %>%
  filter(Year == 2019)
# Get associated weather data
weather.dat <- read.csv("data/BBS/weather.csv")

# Join data with Weather data
# Get date in proper format
weather.dat <- weather.dat %>% 
  unite('date', sep = '-', Year, Month, Day, remove = FALSE)
weather.dat$date <- as.Date(weather.dat$date, tz = "America/New_York")
# Get julian date of each survey
weather.dat$julian <- as.numeric(format(weather.dat$date, '%j'))
weather.covs <- weather.dat %>%
  filter(RouteDataID %in% unique(bbs.dat$RouteDataID))

# Get data in format for spOccupancy --------------------------------------
bbs.dat <- left_join(bbs.dat, weather.covs, by = c('RouteDataID', 'Year'))
# Sort data by species, Year, Longitude, Latitude
bbs.dat <- bbs.dat %>%
  arrange(AOU, Year, Longitude, Latitude)
# bbs.dat %>%
#   group_by(RouteDataID) %>%
#   summarize(n.coords = n_distinct(Longitude)) %>%
#   arrange(desc(n.coords))
coords <- unique(bbs.dat[, c('Longitude', 'Latitude')])
# Number of routes
J <- nrow(coords)
# Create a site index for easy linking. 
coords$site.index <- 1:J
# Join the site index with the full data
bbs.dat <- bbs.dat %>%
  left_join(coords, by = c('Longitude', 'Latitude'))

# Detection Covariates ----------------
# Number of species
N <- n_distinct(bbs.dat$AOU)
sp.codes <- unique(bbs.dat$AOU)
sp.names <- AOU_species_codes$alpha.code[which(AOU_species_codes$spp.num %in% sp.codes)]
# Number of stop replicates (here reduced to 5, not using the full 50 stop data)
K <- 5
# Three dimensional array of detection-nondetection data.
y <- array(NA, dim = c(N, J, K))
# Detection covariates
day <- rep(NA, J)
tod <- rep(NA, J)
obs <- rep(NA, J)
# Unique site indices for linking
site.ordered.indx <- sort(unique(bbs.dat$site.index))
# Calculate
for (j in 1:J) {
  tmp <- bbs.dat %>%
    filter(site.index == site.ordered.indx[j])
  day[j] <- first(tmp$julian)
  obs[j] <- first(tmp$ObsN)
}
for (j in 1:J) {
  print(paste("Currently on site ", j, " out of ", J, sep = ''))
  # Get detection-nondetection data for each species. 
  # group_by is necessary to only grab data from one observer for the rare 
  # occassions where there are multiple observers for a given route
  # in a given year. 
  y[, j, ] <- bbs.dat %>% 
    filter(site.index == site.ordered.indx[j]) %>% 
    arrange(AOU) %>% 
    dplyr::select(starts_with('Count'), -CountryNum) %>%
    as.matrix()
}

y <- ifelse(y > 0, 1, 0)

# Remove species with less than 30 observation or species with less than 50 spatial locations
keep.sp <- which((apply(y, 1, sum, na.rm = TRUE) > 30) & 
                 (apply(apply(y, c(1, 2), function(a) sum(a > 0, na.rm = TRUE)), 1, 
			 function(a) sum(a > 0)) > 50))
y <- y[keep.sp, , ]
sp.names <- as.character(sp.names[keep.sp])

# Reorder species to aid in mixing of MCMC chains.
start.sp <- c('WEME', 'EAME', 'UPSA', 'GRSP', 'LBCU')
# Other species codes
indices <- rep(NA, 5)
for (i in 1:5) {
  indices[i] <- which(sp.names == start.sp[i])
}
indices.other <- 1:nrow(y)
indices.other <- indices.other[-indices]
# Create the ordered y data frame
y <- y[c(indices, indices.other), , ]
# Updated species codes
sp.names <- sp.names[c(indices, indices.other)]
dimnames(y)[[1]] <- sp.names

# Remove nocturnal species
nocturnal.sp <- which(sp.names %in% c('BUOW', 'SEOW'))
y <- y[-nocturnal.sp, , ]
sp.names <- sp.names[-nocturnal.sp]

det.covs <- list(day = day, 
		 obs = obs, 
                 stop = matrix(1:dim(y)[3], dim(y)[2], dim(y)[3], byrow = TRUE))
coords <- as.matrix(coords[c('Longitude', 'Latitude')])

# Save the data -----------------------------------------------------------
save(coords, y, det.covs, sp.codes = sp.names, file = "data/bbs-data-y-det-covs.rda")
