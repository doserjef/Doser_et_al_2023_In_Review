# main-joint.R: script to fit the joint model with no spatially-varying 
#               coefficients for the predictive performance assessment.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)
library(pROC)

# Read in data ------------------------------------------------------------
load('data/hold-out-fit-data.rda')
sp.names <- dimnames(data.list.fit$y)[[1]]

# Run the model -----------------------------------------------------------
# Change as desired
n.batch <- 8000
n.burn <- 100000
n.thin <- 100
n.chains <- 1

# Set priors and initial values.
tuning.list <- list(phi = 0.3)
dist.bbs <- dist(data.list.fit$coords)
min.dist <- min(dist.bbs)
low.dist <- quantile(dist.bbs, 0.05)
max.dist <- max(dist.bbs)
mean.dist <- mean(dist.bbs)
inits <- list(phi = 3 / mean.dist, sigma.sq.p = 1, fix = TRUE)

priors <- list(phi.unif = list(a = 3 / max.dist, b = 3 / low.dist))

out <- svcMsPGOcc(occ.formula = ~ scale(grass) + scale(tmax),
                 det.formula = ~ scale(day) + I(scale(day)^2) +
                                 factor(stop),
               data = data.list.fit,
	       inits = inits,
               priors = priors,
               n.batch = n.batch,
               batch.length = 25,
               accept.rate = 0.43,
               cov.model = "exponential",
               tuning = tuning.list,
	       std.by.sp = TRUE,
               n.omp.threads = 5, # Change as desired
               svc.cols = 1,
               verbose = TRUE,
               NNGP = TRUE,
	       n.factors = 6,
               n.neighbors = 10,
               n.report = 1,
               n.burn = n.burn,
               n.thin = n.thin,
               n.chains = n.chains)

# Predict at hold out locations -------------------------------------------
# NOTE: you should only assess AUC at locations within a given species range, 
#       since you're assuming those data are fixed.
load('data/hold-out-pred-data.rda')
coords.0 <- data.list.pred$coords
X.0 <- as.matrix(data.frame(intercept = 1,
		  grass = (data.list.pred$occ.covs$grass - mean(data.list.fit$occ.covs$grass)) /
			  sd(data.list.fit$occ.covs$grass),
		  tmax = (data.list.pred$occ.covs$tmax - mean(data.list.fit$occ.covs$tmax)) /
			  sd(data.list.fit$occ.covs$tmax)))
out.pred <- predict(out, X.0, coords.0)
# Calculate AUC -----------------------
n.samples <- out$n.post * out$n.chains
y.sum <- apply(data.list.pred$y, c(1, 2), function(a) ifelse(sum(a, na.rm = TRUE) > 0, 1, 0))
N <- nrow(y.sum)
auc.vals <- matrix(NA, n.samples, N)
for (i in 1:N) {
  print(paste0("Currently on species ", i, " out of ", N))
  curr.indx <- which(data.list.pred$range.ind[i, ] == 1)
  for (j in 1:n.samples) {
    auc.vals[j, i] <- auc(response = y.sum[i, curr.indx], 
			  predictor = out.pred$z.0.samples[j, i, curr.indx])
  } # j (iteration)
} # i (species)

save(auc.vals, file = 'results/joint-constant/auc-hold-out-samples.rda')

