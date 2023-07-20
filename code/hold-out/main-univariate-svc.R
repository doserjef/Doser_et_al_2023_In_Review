# main-univariate-svc.R: script to fit a univariate SVC model for use
#                        in the predictive performance assessment.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)
library(pROC)

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
load('data/hold-out-fit-data.rda')
sp.names <- dimnames(data.list.fit$y)[[1]]
sp.indx <- which(sp.names == curr.sp)
site.indx <- which(data.list.fit$range.ind[sp.indx, ] == 1)
data.list.fit$y <- data.list.fit$y[sp.indx, site.indx, ]
data.list.fit$occ.covs <- data.list.fit$occ.covs[site.indx, ]
data.list.fit$det.covs$day <- data.list.fit$det.covs$day[site.indx]
data.list.fit$det.covs$obs <- data.list.fit$det.covs$obs[site.indx]
data.list.fit$det.covs$stop <- data.list.fit$det.covs$stop[site.indx, ]
data.list.fit$coords <- data.list.fit$coords[site.indx, ]


# Run the model -----------------------------------------------------------
n.batch <- 8000
n.burn <- 100000
n.thin <- 100
n.chains <- 1

# Set priors and initial values.
tuning.list <- list(phi = 0.3, sigma.sq = 0.5)
dist.bbs <- dist(data.list.fit$coords)
min.dist <- min(dist.bbs)
low.dist <- quantile(dist.bbs, 0.05)
max.dist <- max(dist.bbs)
mean.dist <- mean(dist.bbs)
inits <- list(phi = 3 / mean.dist, sigma.sq = 1, sigma.sq.p = 1)

priors <- list(phi.unif = list(a = 3 / max.dist, b = 3 / low.dist))

out <- svcPGOcc(occ.formula = ~ scale(grass) + scale(tmax),
                det.formula = ~ scale(day) + I(scale(day)^2) +
                                factor(stop),
               data = data.list.fit,
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

# Predict at hold out locations -------------------------------------------
# NOTE: you're only predicting at sites within a given species range.
load('data/hold-out-pred-data.rda')
coords.0 <- data.list.pred$coords
X.0 <- as.matrix(data.frame(intercept = 1,
		  grass = (data.list.pred$occ.covs$grass - mean(data.list.fit$occ.covs$grass)) /
			  sd(data.list.fit$occ.covs$grass),
		  tmax = (data.list.pred$occ.covs$tmax - mean(data.list.fit$occ.covs$tmax)) /
			  sd(data.list.fit$occ.covs$tmax)))
curr.indx <- which(data.list.pred$range.ind[sp.indx, ] == 1)
X.0 <- X.0[curr.indx, ]
coords.0 <- coords.0[curr.indx, ]
out.pred <- predict(out, X.0, coords.0)
# Calculate AUC -----------------------
n.samples <- out$n.post * out$n.chains
auc.vals <- rep(NA, n.samples)
y.sum.curr <- apply(data.list.pred$y[sp.indx, curr.indx, ], 1, function(a) ifelse(sum(a, na.rm = TRUE) > 0, 1, 0))
auc.vals <- rep(NA, n.samples)
for (j in 1:n.samples) {
  print(j)
  auc.vals[j] <- auc(response = y.sum.curr, predictor = c(out.pred$z.0.samples[j, ]))
} # j (iteration)

# Save AUC values ---------------------------------------------------------
save(auc.vals, file = paste('results/univariate-svc/', curr.sp, '-auc-hold-out-samples.rda',
			    sep = ''))
