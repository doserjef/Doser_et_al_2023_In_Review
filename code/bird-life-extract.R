# bird-life-extract.R: this script extracts data from Bird Life International
#                      for determination of species ranges. Note the BirdLife
#                      data are not available on GitHub, so this script will not
#                      run.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(coda)
library(sf)
library(lwgeom)

# Read in species name data for eastern forest bird community -------------
load("data/full-data-spOccupancy.rda")
sp.names <- dimnames(data.list$y)[[1]]
comm.group.dat <- read.csv("data/bird-species-table-bateman.csv")
scientific.names <- comm.group.dat %>%
  filter(Code %in% sp.names) %>%
  select(Code, name = Scientific.Name, common.name = Common.Name) %>%
  mutate(name = tolower(name)) %>%
  unique()

# Manually adjust species with different classifications
# Baird's Sparrow
scientific.names$name[which(scientific.names$Code == 'BAIS')] <- 'passerculus bairdii'
# Henslow's Sparrow
scientific.names$name[which(scientific.names$Code == 'HESP')] <- 'passerculus henslowii'
# Ammospiza nelsoni (Nelson's Sparrow)
scientific.names$name[which(scientific.names$Code == 'NESP')] <- 'ammospiza nelsoni'
# Le Conte's Sparrow
scientific.names$name[which(scientific.names$Code == 'LCSP')] <- 'ammospiza leconteii'
# Sedge Wren
scientific.names$name[which(scientific.names$Code == 'SEWR')] <- 'cistothorus stellaris'
# White-tailed Hawk
scientific.names$name[which(scientific.names$Code == 'WTHA')] <- 'geranoaetus albicaudatus'

# Read in bird-life data --------------------------------------------------
# What layers are in the GDB? 
st_layers(dsn = '~/Dropbox/DSFBGWZ22/data/BOTW/BOTW.gdb')
# Takes about 30 min or so.
# bird.life <- st_read(dsn = 'data/BOTW/BOTW.gdb', query = "SELECT * FROM All_Species limit 2000" )
bird.life <- st_read(dsn = '~/Dropbox/DSFBGWZ22/data/BOTW/BOTW.gdb', 
		     query = "SELECT * 
		              FROM All_Species
                              WHERE sci_name IN ('passerculus bairdii', 'dolichonyx oryzivorus',
						 'peucaea botterii', 'peucaea cassinii',         
                                                 'calcarius ornatus', 'spizella pallida',         
                                                 'spiza americana', 'tyrannus tyrannus',        
                                                 'sturnella magna', 'buteo regalis',            
                                                 'ammodramus savannarum', 'passerculus henslowii',    
                                                 'eremophila alpestris', 'calamospiza melanocorys',  
                                                 'ammospiza leconteii', 'lanius ludovicianus',      
                                                 'numenius americanus', 'rhynchophanes mccownii',   
                                                 'charadrius montanus', 'ammospiza nelsoni',        
                                                 'passerculus sandwichensis', 'tyrannus forficatus',      
                                                 'cistothorus stellaris', 'anthus spragueii',         
                                                 'buteo swainsoni', 'bartramia longicauda',     
                                                 'pooecetes gramineus', 'tyrannus verticalis',      
                                                 'sturnella neglecta', 'geranoaetus albicaudatus')") 
# Save data ---------------------------------------------------------------
save(bird.life, file = 'data/bird-life-data.rda')


