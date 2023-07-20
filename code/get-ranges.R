# get-ranges.R: code to extract species ranges from BirdLife International.
#               Note this script will not run since the BirdLife International
#               data are only available via a data sharing agreement.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(coda)
library(sf)
library(lwgeom)

# Read in data for grassland bird community -------------------------------
load("data/bbs-data-y-det-covs.rda")
sp.names <- dimnames(y)[[1]]
comm.group.dat <- read.csv("data/bird-species-table-bateman.csv")
scientific.names <- comm.group.dat %>%
  filter(Code %in% sp.names) %>%
  select(Code, name = Scientific.Name, common.name = Common.Name) %>%
  mutate(name = tolower(name)) %>%
  unique()

# Manually adjust species with different classifications
# Sedge Wren
scientific.names$name[which(scientific.names$Code == 'SEWR')] <- 'cistothorus stellaris'

# Read in extracted bird-life data ----------------------------------------
# Object read in is called bird.life
load("data/bird-life-data.rda")

# Only extract data for the eastern US
coords.sf <- st_as_sf(as.data.frame(coords),
		      coords = c("Longitude", "Latitude"),
		      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
coords.sf <- coords.sf %>%
  st_transform("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")

st_crs(coords.sf) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# Map of US states
usa <- st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(st_crs(coords.sf))
usa <- st_union(st_make_valid(usa))

# Extract necessary information -------------------------------------------
# sf object for the area within the species range.
within.sf.list <- list()
# sf object for the area outside the species range. 
outside.sf.list <- list()
bird.life.sp <- unique(bird.life$sci_name)
for (i in 1:length(bird.life.sp)) {
  print(i)
  tmp.indx <- which(scientific.names$name == tolower(bird.life.sp[i]))
  if (length(tmp.indx) == 0) next
  indx <- which(sp.names == scientific.names$Code[tmp.indx])
  print(scientific.names$common.name[tmp.indx])
  # Presence == 1 means extant range, seasonal %in% c(1, 2) means
  # resident or breeding season ranges.
  curr.full <- bird.life %>%
    filter(sci_name == bird.life.sp[i], 
           presence == 1, seasonal %in% c(1, 2))
  # Need to manually do some filtering for one species
  # Grasshopper sparrow, just extract breeding range.
  if (scientific.names$Code[tmp.indx] == 'GRSP') {
    curr.full <- bird.life %>%
      filter(sci_name == bird.life.sp[i], 
	     presence == 1, seasonal == 2)    
  }
  curr.sp.range <- curr.full %>%
    st_transform(st_crs(coords.sf))
  curr.sp.range <- st_union(st_make_valid(curr.sp.range))
  filtered.range <- st_intersection(curr.sp.range, usa)
  # Buffering the range by 200km
  filtered.range <- st_buffer(filtered.range, 200)
  filtered.range <- st_intersection(filtered.range, usa)
  within.sf.list[[indx]] <- filtered.range
}

N <- nrow(y)
J <- ncol(y)
range.ind <- matrix(0, N, J)
for (i in 1:N) {
  indx <- unlist(c(st_contains(within.sf.list[[i]], coords.sf)))
  range.ind[i, indx] <- 1
}

save(within.sf.list, range.ind, file = "data/range.data.rda")
