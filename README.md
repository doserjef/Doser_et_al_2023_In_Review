# Modeling complex species-environment relationships through spatially-varying coefficient occupancy models

### In Review

### Jeffrey W. Doser, Andrew O. Finley, Sarah P. Saunders, Marc K&eacute;ry, Aaron S. Weed, Elise F. Zipkin

### spOccupancy Package [Website](https://www.jeffdoser.com/files/spoccupancy-web/) and [Repository](https://github.com/doserjef/spOccupancy/)

### Please contact the first author for questions about the code or data used in the manuscript: Jeffrey W. Doser (doserjef@msu.edu)

---------------------------------

## Abstract

Occupancy models are frequently used by ecologists to quantify spatial variation in species distributions while accounting for observational biases (e.g., false-negative errors) in collection of detection-nondetection data. However, the common assumption that a single set of regression coefficients can adequately explain species-environment relationships is often biologically unrealistic across large spatial domains. Here we develop computationally-efficient single-species (i.e., univariate) and multi-species (i.e., multivariate) spatially-varying coefficient (SVC) occupancy models to account for spatially-varying species-environment relationships. We employ Nearest Neighbor Gaussian Processes and \pg data augmentation in a hierarchical Bayesian framework to yield computationally efficient Gibbs samplers, which we implement in the spOccupancy R package. For multi-species models, we use spatial factor dimension reduction to efficiently model datasets with large numbers of species (e.g., > 10) while simultaneously accounting for species correlations. The hierarchical Bayesian framework readily enables generation of posterior predictive maps of the SVCs, with fully propagated uncertainty. We apply our SVC models to quantify spatial variability in the relationship between maximum breeding season temperature and occurrence probability of 21 grassland bird species across the US.  Jointly modeling species in a multivariate model generally outperformed univariate models, which all revealed substantial spatial variability in species occurrence relationships with maximum temperatures. Our models are particularly relevant for quantifying species-environment relationships using detection-nondetection data from large-scale monitoring programs, which are becoming increasingly prevalent for answering macroscale ecological questions regarding wildlife responses to global change.

## Repository Directory

All model results were generated using `spOccupancy` v0.7.0. This version is available for download from [GitHub](https://github.com/doserjef/spOccupancy) and will be available on CRAN in mid-August 2023.  

### [code](./code)

Contains all code used in the BBS case study and simulation study

+ `bbs-data-prep.R`: preps raw BBS data for analysis.
+ `bird-life-extract.R`: extracts data from BirdLife International for determination of species ranges.
+ `calculate-waic.R`: script to calculate WAIC for all species from all candidate models.
+ `extract-results.R`: script to extract posterior quantiles and other information used to generate figures and summarize results for the top performing model.
+ `format-pred-data.R`: script to format prediction data for use in predicting species-specific occupancy and covariate effects across the continental US
+ `get-prediction-coords.R`: script to extract the prediction grid across the continental US.
+ `get-ranges.R`: script to extract the species ranges from BirdLife for restricting the area across which models are run for each individual species.
+ `hold-out/`: directory containing all scripts to run out-of-sample predictive performance assessment for the eight candidate models.
+ `lulc-data-pre.R`: script to calculate grassland area at BBS routes using USGS EROS data.
+ `lulc-data-pre.R`: script to calculate grassland area across the prediction grid using USGS EROS data.
+ `main-*.R`: scripts to fit the different candidate models, where the name of the candidate model replaces `*`. 
+ `pred-joint-svc-interaction.R`: script to predict occupancy probability and spatially-varying coefficients across the continental US using the top performing model (the multi-species model with an SVC and interaction).
+ `sims/main-sim.R`: script to run the simulation study proof of concept.
+ `sims/summary-sim.R`: script to summarize the results from the simulation study proof of concept.
+ `spOcc-format-data.R`: script to format data into spOccupancy format.
+ `summary.R`: summarizes results from all models fit in the case study and produces all case study figures in the manuscript.
+ `tmax-data-prep.R`: script to extract maximum temperature data from PRISM at BBS routes.
+ `tmax-pred-data-prep.R`: script to extract maximum temperature data from PRISM across the prediction grid.

### [data](./data)

Contains data used for the BBS case study.

+ `hold-out-fit-data.rda`: subset of data in spOccupancy format for fitting models in the hold-out assessment. 
+ `hold-out-pred-data.rda`: subset of data in spOccupancy format for prediction in the hold-out assessment.
+ `bird-species-table-bateman.csv`: CSV file from Bateman et al. (2020) containing information on species classifications that was used to select the grassland bird community.
+ `pred-coordinates.rda`: coordinates that form a grid across the US for prediction.
+ `full-pred-data.rda`: formatted prediction data for prediction across the US.
+ `full-data-spOccupancy.rda`: BBS data and covariates formatted for use in spOccupancy model fitting functions.
