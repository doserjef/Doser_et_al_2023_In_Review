# main-joint.R: script to fit the joint model with no spatially-varying
#               coefficients.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

# Read in data ------------------------------------------------------------
load('data/full-data-spOccupancy.rda')
sp.names <- dimnames(data.list$y)[[1]]

# Run the model -----------------------------------------------------------
n.batch <- 20000
n.burn <- 200000
n.thin <- 150
n.chains <- 1

# Set priors and initial values.
tuning.list <- list(phi = 0.3)
dist.bbs <- dist(data.list$coords)
min.dist <- quantile(dist.bbs, 0.05)
max.dist <- max(dist.bbs)
mean.dist <- mean(dist.bbs)
inits <- list(phi = 3 / mean.dist, fix = TRUE)

priors <- list(phi.unif = list(a = 3 / max.dist, b = 3 / min.dist))

out <- svcMsPGOcc(occ.formula = ~ scale(grass) + scale(tmax),
                 det.formula = ~ scale(day) + I(scale(day)^2) +
                                 factor(stop),
               data = data.list,
	       inits = inits,
               priors = priors,
               n.batch = n.batch,
               svc.cols = 1,
               batch.length = 25,
               accept.rate = 0.43,
               cov.model = "exponential",
               std.by.sp = TRUE,
               tuning = tuning.list,
               n.omp.threads = 5, # Change as desired
               verbose = TRUE,
               NNGP = TRUE,
	       n.factors = 6,
               n.neighbors = 10,
               n.report = 1,
               n.burn = n.burn,
               n.thin = n.thin,
               n.chains = n.chains)


# Save results ------------------------------------------------------------
save(out, file = 'results/joint-constant/joint-constant-results.rda')
