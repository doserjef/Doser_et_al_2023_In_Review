# pred-joint-svc-interaction.R: script to predict species-specific 
#                               occupancy probability and the spatially-varying
#                               coefficient across the continental US. Note this
#                               script requires the full object from the resulting
#                               model fit, which is not available on GitHub since
#                               it is quite large.
# Author: Jeffrey W. Doser 
rm(list = ls())
library(spOccupancy)

# Load the data used to fit the model -------------------------------------
load('data/full-data-spOccupancy.rda')
sp.names <- dimnames(data.list$y)[[1]]

# Load the model fit object -----------------------------------------------
load('results/joint-svc-interaction/joint-svc-interaction-results.rda')

# Load the prediction data ------------------------------------------------
load('data/full-pred-data.rda')
# Standardize/scale variables by those used to fit the model
# Note that the spOccupancy predict function for svcMsPGOcc models will 
# automatically standardize the covariate values appropriately for each species
# based on the model used to fit the data if the model was fit with std.by.sp = TRUE.
J.0 <- nrow(coords.0)
X.0 <- matrix(1, J.0, ncol(out$X))
colnames(X.0) <- c('(Intercept)', 'scale(grass)', 'scale(tmax)', 'grassTmaxInt')
X.0[, 'scale(grass)'] <- (grass.0 - mean(data.list$occ.covs$grass)) / sd(data.list$occ.covs$grass)
X.0[, 'scale(tmax)'] <- (tmax.0 - mean(data.list$occ.covs$tmax)) / sd(data.list$occ.covs$tmax)
X.0[, 'grassTmaxInt'] <- X.0[, 'scale(grass)'] * X.0[, 'scale(tmax)'] 
vals <- split(1:J.0, ceiling(seq_along(1:J.0) / 500))
N <- nrow(data.list$y)
n.samples <- out$n.post * out$n.chains
psi.quants <- array(NA, dim = c(3, N, J.0))
q <- out$q
p.svc <- length(out$svc.cols) 
w.samples <- array(NA, dim = c(n.samples, q, J.0, p.svc))
for (j in 1:length(vals)) {
  print(paste("Currently on set ", j, " out of ", length(vals), sep = ''))
  curr.indx <- vals[[j]]
  out.pred <- predict(out, X.0[curr.indx, ], coords.0[curr.indx, ],
        	      n.omp.threads = 10, verbose = TRUE)
  w.samples[, , curr.indx, ] <- out.pred$w.0.samples
  psi.quants[, , curr.indx] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))
}

# Save samples and quantiles ----------------------------------------------
save(w.samples, psi.quants, file = 'results/joint-svc-interaction/pred-results.rda')
