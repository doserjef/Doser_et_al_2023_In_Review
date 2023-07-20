# extract-results.R: script to extract the main results from the top
#                    performing model (the joint SVC model with 
#                    an interaction) for subsequent analyses. This is
#                    necessary since the resulting model files from 
#                    spOccupancy are quite large and too big for GitHub.
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
J.0 <- nrow(coords.0)
X.0 <- matrix(1, J.0, ncol(out$X))
X.0 <- matrix(1, J.0, ncol(out$X))
colnames(X.0) <- c('(Intercept)', 'scale(grass)', 'scale(tmax)', 'grassTmaxInt')
X.0[, 'scale(grass)'] <- (grass.0 - mean(data.list$occ.covs$grass)) / sd(data.list$occ.covs$grass)
X.0[, 'scale(tmax)'] <- (tmax.0 - mean(data.list$occ.covs$tmax)) / sd(data.list$occ.covs$tmax)
X.0[, 'grassTmaxInt'] <- X.0[, 'scale(grass)'] * X.0[, 'scale(tmax)']

# Read in prediction results ----------------------------------------------
load("results/joint-svc-interaction/pred-results.rda")

# Generate different forms of the SVCs for plotting -----------------------
n.samples <- out$n.post * out$n.chains
J <- nrow(coords.0)
N <- nrow(out$y)
p.occ <- ncol(out$X)
q <- out$q
p.svc <- length(out$svc.cols)
svc.cols <- out$svc.cols
# Save non-spatial species-specific regression coefficients before continuing
beta.samples <- out$beta.samples
save(beta.samples, file = 'results/joint-svc-interaction/beta-samples.rda')
beta.samples <- array(out$beta.samples, dim = c(n.samples, N, p.occ))
lambda.samples <- array(out$lambda.samples, dim = c(n.samples, N, q, p.svc))

# Scale and center grass by species-specific means and sds
# NOTE: hardcoding here
grass.list <- list()
for (i in 1:N) {
  grass.list[[i]] <- (X.0[, 2] - out$species.means[i, 2]) / out$species.sds[i, 2]
}

# Spatial component of maximum temperature effect
beta.star.tmax.samples <- array(NA, dim = c(n.samples, N, J))
beta.star.int.samples <- array(NA, dim = c(n.samples, N, J))
# NOTE: lots of hardcoding here. 2
for (j in 1:n.samples) {
  print(j)
  for (i in 1:N) {
    tmp.1 <- lambda.samples[j, i, , 2] %*% w.samples[j, , , 2]
    tmp.2 <- t(apply(as.matrix(beta.samples[j, i, 4]), 1,   
                                    function(a) a * grass.list[[i]]))
    tmp.3 <- beta.samples[j, i, 3]
    beta.star.tmax.samples[j, i, ] <- tmp.1 + tmp.2 + tmp.3

    tmp.1 <- lambda.samples[j, i, , 1] %*% w.samples[j, , , 1]
    tmp.3 <- beta.samples[j, i, 1]
    beta.star.int.samples[j, i, ] <- tmp.1 + tmp.3
  }
}

# Calculate the posterior quantiles for each species
beta.star.tmax.quants <- apply(beta.star.tmax.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))
beta.star.tmax.prob.pos <- apply(beta.star.tmax.samples, c(2, 3), function(a) mean(a > 0))
beta.star.int.quants <- apply(beta.star.int.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))

# Save quantiles and probability of positive effect to object
save(beta.star.tmax.quants, beta.star.int.quants, 
     beta.star.tmax.prob.pos, psi.quants, file = 'results/joint-svc-interaction/pred-quantiles.rda')
