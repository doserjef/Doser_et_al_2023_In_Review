# main-univariate-svc-interaction.R: script to fit a univariate SVC model
#                                    with an interaction between maximum
#                                    temperature and grassland cover for
#                                    the predictive performance assessment.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

# Get info from command line ----------------------------------------------
# This code is to extract the current species name from the command line
# to easily run the script for different species
args <- commandArgs(trailingOnly = TRUE)
# Current species
curr.sp <- args[1]
# Alternatively, if not running the script from the command line:
# curr.sp <- 'LCSP'
if(length(args) == 0) base::stop('Need to tell spOccupancy the current species')
print(curr.sp)

# Read in data ------------------------------------------------------------
load('data/full-data-spOccupancy.rda')
sp.names <- dimnames(data.list$y)[[1]]
sp.indx <- which(sp.names == curr.sp)
site.indx <- which(data.list$range.ind[sp.indx, ] == 1)
data.list$y <- data.list$y[sp.indx, site.indx, ]
data.list$occ.covs <- data.list$occ.covs[site.indx, ]
data.list$det.covs$day <- data.list$det.covs$day[site.indx]
data.list$det.covs$obs <- data.list$det.covs$obs[site.indx]
data.list$det.covs$stop <- data.list$det.covs$stop[site.indx, ]
data.list$coords <- data.list$coords[site.indx, ]

# Run the model -----------------------------------------------------------
# Change as desired
n.batch <- 8000
n.burn <- 100000
n.thin <- 100
n.chains <- 3

# Set priors and initial values.
tuning.list <- list(phi = 0.3, sigma.sq = 0.5)
dist.bbs <- dist(data.list$coords)
min.dist <- quantile(dist.bbs, 0.05)
max.dist <- max(dist.bbs)
mean.dist <- mean(dist.bbs)
inits <- list(phi = 3 / mean.dist, sigma.sq = 1, sigma.sq.p = 1)

priors <- list(phi.unif = list(a = 3 / max.dist, b = 3 / min.dist))

out <- svcPGOcc(occ.formula = ~ scale(grass) * scale(tmax),
                det.formula = ~ scale(day) + I(scale(day)^2) +
                                factor(stop),
               data = data.list,
	       inits = inits,
               priors = priors,
               n.batch = n.batch,
	       svc.cols = c(1, 3),
               batch.length = 25,
               accept.rate = 0.43,
               cov.model = "exponential",
               tuning = tuning.list,
               n.omp.threads = 1, # Change as desired
               verbose = TRUE,
               NNGP = TRUE,
               n.neighbors = 10,
               n.report = 10,
               n.burn = n.burn,
               n.thin = n.thin,
               n.chains = n.chains)

# Save to hard drive ------------------------------------------------------
save(out, file = paste('results/univariate-svc-interaction/', curr.sp, 
		       '-univariate-svc-interaction-results.rda', sep = ''))
