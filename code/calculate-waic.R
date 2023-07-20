# calculate-waic.R: this script calculates the WAIC for each individual
#                   species and each candidate model and saves the objects
#                   in subsequent files. Full model objects are too large
#                   for GitHub, and so this script will not run successfully, 
#                   but is provided to show how WAIC can be calculated with 
#                   spOccupancy::waicOcc()
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

# Load full model results -------------------------------------------------
load("results/joint-constant/joint-constant-results.rda")
sp.names <- out$sp.names
waic.joint.constant <- waicOcc(out, by.sp = TRUE)
rownames(waic.joint.constant) <- sp.names
save(waic.joint.constant, file = 'results/joint-constant/waic-joint-constant.rda')
load("results/joint-interaction/joint-interaction-results.rda")
sp.names <- out$sp.names
waic.joint.interaction <- waicOcc(out, by.sp = TRUE)
rownames(waic.joint.interaction) <- sp.names
save(waic.joint.interaction, file = 'results/joint-interaction/waic-joint-interaction.rda')
load("results/joint-svc/joint-svc-results.rda")
sp.names <- out$sp.names
waic.joint.svc <- waicOcc(out, by.sp = TRUE)
rownames(waic.joint.svc) <- sp.names
save(waic.joint.svc, file = 'results/joint-svc/waic-joint-svc.rda')
load("results/joint-svc-interaction/joint-svc-interaction-results.rda")
sp.names <- out$sp.names
waic.joint.svc.interaction <- waicOcc(out, by.sp = TRUE)
rownames(waic.joint.svc.interaction) <- sp.names
save(waic.joint.svc.interaction, file = 'results/joint-svc-interaction/waic-joint-svc-interaction.rda')


waic.uni.constant <- array(NA, dim = dim(waic.joint.constant))
waic.uni.interaction <- array(NA, dim = dim(waic.joint.constant))
waic.uni.svc <- array(NA, dim = dim(waic.joint.constant))
waic.uni.svc.interaction <- array(NA, dim = dim(waic.joint.constant))
for (i in 1:length(sp.names)) {
  print(paste0("Currently on species ", i, " out of ", length(sp.names)))
  load(paste0("results/univariate-constant/", sp.names[i], "-univariate-constant-results.rda"))
  waic.uni.constant[i, ] <- waicOcc(out)
  load(paste0("results/univariate-interaction/", sp.names[i], "-univariate-interaction-results.rda"))
  waic.uni.interaction[i, ] <- waicOcc(out)
  load(paste0("results/univariate-svc/", sp.names[i], "-univariate-svc-results.rda"))
  waic.uni.svc[i, ] <- waicOcc(out) 
  load(paste0("results/univariate-svc-interaction/", sp.names[i], "-univariate-svc-interaction-results.rda"))
  waic.uni.svc.interaction[i, ] <- waicOcc(out)
}
save(waic.uni.constant, file = 'results/univariate-constant/waic-univariate-constant-results.rda')
save(waic.uni.interaction, file = 'results/univariate-interaction/waic-univariate-interaction-results.rda')
save(waic.uni.svc, file = 'results/univariate-svc/waic-univariate-svc-results.rda')
save(waic.uni.svc.interaction, file = 'results/univariate-svc-interaction/waic-univariate-svc-interaction-results.rda')
