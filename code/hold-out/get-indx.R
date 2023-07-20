# get-indx.R: script to extract the hold-out locations used for model fitting
#             and prediction in the predictive performance assessment using AUC.
# Author: Jeffrey W. Doser
rm(list = ls())
set.seed(10101)
load("data/full-data-spOccupancy.rda")

J <- nrow(data.list$coords)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
data.list.fit <- data.list
data.list.pred <- data.list

data.list.fit$y <- data.list$y[, -pred.indx, ]
data.list.fit$occ.covs <- data.list$occ.covs[-pred.indx, ]
data.list.fit$det.covs$day <- data.list$det.covs$day[-pred.indx]
data.list.fit$det.covs$obs <- data.list$det.covs$obs[-pred.indx]
data.list.fit$det.covs$stop <- data.list$det.covs$stop[-pred.indx, ]
data.list.fit$coords <- data.list$coords[-pred.indx, ]
data.list.fit$range.ind <- data.list$range.ind[, -pred.indx]

data.list.pred$y <- data.list$y[, pred.indx, ]
data.list.pred$occ.covs <- data.list$occ.covs[pred.indx, ]
data.list.pred$det.covs$day <- data.list$det.covs$day[pred.indx]
data.list.pred$det.covs$obs <- data.list$det.covs$obs[pred.indx]
data.list.pred$det.covs$stop <- data.list$det.covs$stop[pred.indx, ]
data.list.pred$coords <- data.list$coords[pred.indx, ]
data.list.pred$range.ind <- data.list$range.ind[, pred.indx]

apply(data.list.fit$range.ind, 1, sum)
apply(data.list.pred$range.ind, 1, sum)
apply(data.list.pred$range.ind, 1, sum) / apply(data.list$range.ind, 1, sum)


save(data.list.fit, file = 'data/hold-out-fit-data.rda')
save(data.list.pred, file = 'data/hold-out-pred-data.rda')
