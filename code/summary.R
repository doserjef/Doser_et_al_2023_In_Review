# summary.R: script to summarize case study results and generate all figures
#            in the manuscript and supplemental material
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(spOccupancy)
library(viridis)
library(pals)
library(ggpubr)
library(stars)

load("data/full-data-spOccupancy.rda")

# Coordinates and maps for plotting ---------------------------------------
coords.sf <- st_as_sf(as.data.frame(data.list$coords),
		      coords = c("X", "Y"),
		      crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")

# Map of US states
usa <- st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(st_crs(coords.sf))
usa <- st_union(st_make_valid(usa))

# Summarize hold out assessments with AUC ---------------------------------
sp.names <- dimnames(data.list$y)[[1]]
N <- length(sp.names)

# Read in auc values for all the candidate models
load("results/joint-constant/auc-hold-out-samples.rda")
auc.joint.constant <- auc.vals
colnames(auc.joint.constant) <- sp.names
load("results/joint-interaction/auc-hold-out-samples.rda")
auc.joint.interaction <- auc.vals
colnames(auc.joint.interaction) <- sp.names
load("results/joint-svc/auc-hold-out-samples.rda")
auc.joint.svc <- auc.vals
colnames(auc.joint.svc) <- sp.names
load("results/joint-svc-interaction/auc-hold-out-samples.rda")
auc.joint.svc.interaction <- auc.vals
colnames(auc.joint.svc.interaction) <- sp.names
auc.uni.constant <- array(NA, dim = dim(auc.joint.constant))
auc.uni.interaction <- array(NA, dim = dim(auc.joint.constant))
auc.uni.svc <- array(NA, dim = dim(auc.joint.constant))
auc.uni.svc.interaction <- array(NA, dim = dim(auc.joint.constant))
for (i in 1:N) {
  load(paste0("results/univariate-constant/", sp.names[i], "-auc-hold-out-samples.rda"))
  auc.uni.constant[, i] <- auc.vals
  load(paste0("results/univariate-interaction/", sp.names[i], "-auc-hold-out-samples.rda"))
  auc.uni.interaction[, i] <- auc.vals
  load(paste0("results/univariate-svc/", sp.names[i], "-auc-hold-out-samples.rda"))
  auc.uni.svc[, i] <- auc.vals
  load(paste0("results/univariate-svc-interaction/", sp.names[i], "-auc-hold-out-samples.rda"))
  auc.uni.svc.interaction[, i] <- auc.vals
}
colnames(auc.uni.constant) <- sp.names
colnames(auc.uni.interaction) <- sp.names
colnames(auc.uni.svc) <- sp.names
colnames(auc.uni.svc.interaction) <- sp.names


# Create a data frame with the means for each model
auc.df <- data.frame(u.const = apply(auc.uni.constant, 2, mean),
		     u.interaction = apply(auc.uni.interaction, 2, mean),
		     u.svc = apply(auc.uni.svc, 2, mean),
		     u.svc.interaction = apply(auc.uni.svc.interaction, 2, mean), 
		     joint.const = apply(auc.joint.constant, 2, mean),
		     joint.interaction = apply(auc.joint.interaction, 2, mean),
		     joint.svc = apply(auc.joint.svc, 2, mean),
		     joint.svc.interaction = apply(auc.joint.svc.interaction, 2, mean))
apply(auc.df, 1, which.max)

auc.max.by.sp <- apply(auc.df, 1, max)
auc.plot.df <- apply(auc.df, 1, function(a) a - max(a))

# Generate a plot ---------------------
n.cand <- 8
auc.long.df <- data.frame(auc = c(auc.plot.df), 
			  species = rep(sp.names, each = n.cand), 
                          model = factor(c('SS', 'SS-INT', 
				           'SS-SVC', 'SS-INT-SVC', 
				           'MS', 'MS-INT', 
				           'MS-SVC', 'MS-INT-SVC'), 
			                 ordered = TRUE, 
					 levels = rev(c('MS', 'MS-INT', 
				           'SS', 'SS-INT', 
				           'SS-SVC', 'SS-INT-SVC', 
				           'MS-SVC', 'MS-INT-SVC'))))
auc.long.df$val <- c(t(auc.df))
for (i in 1:N) {
  auc.long.df[which(auc.long.df$species == sp.names[i]), 'val'] <- ifelse(auc.long.df[which(auc.long.df$species == sp.names[i]), 'val'] == auc.max.by.sp[i], auc.long.df[which(auc.long.df$species == sp.names[i]), 'val'], NA)
}
# Generate Figure S2
ggplot(auc.long.df, aes(x = model, y = species, fill = auc)) + 
  geom_tile(color = 'black') +
  scale_fill_gradientn(expression(paste(Delta, " AUC")), colors = rev(brewer.reds(1000)),
                       guide = guide_colourbar(title.position="top", reverse = FALSE),
                       na.value = NA) +
  theme_bw(base_size = 18) + 
  geom_text(aes(label = round(val, digits = 2)), family = 'LM Roman 10', size = 4) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  labs(x = 'Model', y = 'Species', fill = 'Trend') + 
  theme(axis.ticks.y = element_blank(), 
	axis.ticks.x = element_blank(), 
	legend.title = element_text(size = 14),
        text = element_text(family="LM Roman 10"), 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file = 'figures/Figure-S2.png', width = 7, height = 8, units = 'in', bg = 'white')

# Summarize model fit with WAIC -------------------------------------------

# Read in waic values for all models
load("results/joint-constant/waic-joint-constant.rda")
load("results/joint-interaction/waic-joint-interaction.rda")
load("results/joint-svc/waic-joint-svc.rda")
load("results/joint-svc-interaction/waic-joint-svc-interaction.rda")
load("results/univariate-constant/waic-univariate-constant-results.rda")
load("results/univariate-interaction/waic-univariate-interaction-results.rda")
load("results/univariate-svc/waic-univariate-svc-results.rda")
load("results/univariate-svc-interaction/waic-univariate-svc-interaction-results.rda")

# Create a data frame with the means for each model
waic.df <- data.frame(u.const = waic.uni.constant[, 3],
		     u.interaction = waic.uni.interaction[, 3],
		     u.svc = waic.uni.svc[, 3],
		     u.svc.interaction = waic.uni.svc.interaction[, 3], 
		     joint.const = waic.joint.constant[, 3],
		     joint.interaction = waic.joint.interaction[, 3],
		     joint.svc = waic.joint.svc[, 3],
		     joint.svc.interaction = waic.joint.svc.interaction[, 3])
apply(waic.df, 1, which.min)

# Generate plot of results ------------
waic.min.by.sp <- apply(waic.df, 1, min)
waic.plot.df <- apply(waic.df, 1, function(a) (a - min(a)))
colnames(waic.plot.df) <- sp.names
# waic.plot.df <- apply(waic.df, 1, function(a) (a - min(a)) / min(a))
waic.sum.all <- apply(waic.df, 2, sum)
waic.plot.df <- cbind(waic.plot.df, waic.sum.all)


# Number of candidate models
n.cand <- 8
waic.long.df <- data.frame(waic = c(waic.plot.df), 
			  species = factor(rep(c(sp.names, 'COMM'), each = n.cand), 
					   levels = c('COMM', sort(sp.names)), order = TRUE),
                          model = factor(c('SS', 'SS-INT', 
				           'SS-SVC', 'SS-INT-SVC', 
				           'MS', 'MS-INT', 
				           'MS-SVC', 'MS-INT-SVC'), 
			                 ordered = TRUE, 
					 levels = rev(c('MS', 'MS-INT', 
				           'SS', 'SS-INT', 
				           'SS-SVC', 'SS-INT-SVC', 
				           'MS-SVC', 'MS-INT-SVC'))))
waic.long.df$val <- NA
for (i in 1:N) {
  waic.long.df[which(waic.long.df$species == sp.names[i]), 'val'] <- ifelse(waic.long.df[which(waic.long.df$species == sp.names[i]), 'waic'] < 2, '*', NA)
  waic.long.df[which(waic.long.df$species == sp.names[i]), 'waic'] <- waic.long.df[which(waic.long.df$species == sp.names[i]), 'waic'] / waic.min.by.sp[i]
}
# Do it separately for the community sum
waic.long.df[(nrow(waic.long.df) - n.cand + 1):nrow(waic.long.df), 'waic'] <- (waic.long.df[(nrow(waic.long.df) - n.cand + 1):nrow(waic.long.df), 'waic'] - min(waic.sum.all))
waic.long.df[(nrow(waic.long.df) - n.cand + 1):nrow(waic.long.df), 'val'] <- ifelse(waic.long.df[(nrow(waic.long.df) - n.cand + 1):nrow(waic.long.df), 'waic'] < 2, '*', NA)
waic.long.df[(nrow(waic.long.df) - n.cand + 1):nrow(waic.long.df), 'waic'] <- (waic.long.df[(nrow(waic.long.df) - n.cand + 1):nrow(waic.long.df), 'waic'] / min(waic.sum.all))

# Figure 2
ggplot(waic.long.df, aes(x = model, y = species, fill = waic)) + 
  geom_tile(color = 'black') +
  scale_fill_gradientn("Proportional\nchange in WAIC", colors = brewer.reds(1000),
                       guide = guide_colourbar(title.position="top", reverse = FALSE),
                       na.value = NA) +
  theme_bw(base_size = 18) + 
  geom_text(aes(label = val), family = 'LM Roman 10', size = 4) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  labs(x = 'Model', y = 'Species', fill = 'Trend') + 
  theme(axis.ticks.y = element_blank(), 
	axis.ticks.x = element_blank(), 
	legend.title = element_text(size = 14),
        text = element_text(family="LM Roman 10"), 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file = 'figures/Figure-2.png', width = 7, height = 8, units = 'in', bg = 'white')

# Summarize model results from best performing model ----------------------
# The full posterior distributions for all model parameters are too big for GitHub, 
# so instead this loads quantiles for model objects that are shown in plots.
# This loads the following objects: 
#      beta.tmax.int.quants (the interaction adjustment from the overall mean tmax effect), 
#      w.star.tmax.quants (the spatial component of the tmax effect), and 
#      beta.tmax.quants (the non-spatial tmax effect). 
load("results/joint-svc-interaction/pred-quantiles.rda")
# Load prediction grid
load('data/full-pred-data.rda')
J.0 <- nrow(coords.0)
# Species range data from BirdLife (commented out as ranges must be obtained via
# data sharing agreement with BirdLife
load('data/range.data.rda')

# Get coordinates set up and info on which sites are within a species range
coords.0.sf <- st_as_sf(as.data.frame(coords.0),
		      coords = c("X", "Y"),
		      crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")
range.pred.full <- array(0, dim = c(N, J.0))
for (i in 1:N) {
  indx <- unlist(c(st_contains(within.sf.list[[i]], coords.0.sf)))
  range.pred.full[i, indx] <- 1
}
save(range.pred.full, file = 'data/prediction-range-binary.rda')
# Load in the binary value that indicates whether a given pixel is within BirdLife
# range or not. Loads an object called range.pred.full
load('data/prediction-range-binary.rda')

# Summary of effects across all species
# Cutoffs used to generate mapes related to the probability of a positive effect:
#      0-.2, .2-.4, .4-.6, .6-.8, .8-1
prop.lowest <- apply(beta.star.tmax.prob.pos, 1, function(a) mean(a <= .2))
prop.low <- apply(beta.star.tmax.prob.pos, 1, function(a) mean(a > .2 & a <= .4))
prop.med <- apply(beta.star.tmax.prob.pos, 1, function(a) mean(a >.4 & a <=.6))
prop.high <- apply(beta.star.tmax.prob.pos, 1, function(a) mean(a > .6 & a <= .8))
prop.highest <- apply(beta.star.tmax.prob.pos, 1, function(a) mean(a > .8))

# Generate data frame for plotting
plot.df <- data.frame(prop = c(prop.lowest, prop.low, prop.med, prop.high,
			       prop.highest),
		      type = factor(rep(c('Strong Negative', 'Moderate Negative', 'No Support', 
                                          'Moderate Positive', 'Strong Positive'), each = N), 
                                           levels = c('Strong Negative', 'Moderate Negative', 
                                                      'No Support', 'Moderate Positive', 
                                                      'Strong Positive')),
		      species = rep(sp.names, times = 5), 
                      mean.val = rep(apply(beta.star.tmax.quants[2, , ], 1, mean), 
		                     times = 5))
plot.df <- plot.df %>%
  arrange(desc(mean.val))
plot.df$species <- factor(plot.df$species, levels = unique(plot.df$species), order = TRUE)

# Generate Figure 3
ggplot(plot.df, aes(x = species, y = prop, fill = type)) +
  geom_bar(stat = 'identity', width = 1, color = 'grey') +
  theme_bw(base_size = 16) +
  # scale_fill_viridis_d() +
  scale_fill_brewer(palette = 'RdBu') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_blank(),
	text = element_text(family  = 'LM Roman 10'),
 	legend.position = 'bottom') +
  labs(x = 'Species', y = 'Proportion of Sites',
       fill = '')
ggsave(file = 'figures/Figure-3.png', units = 'in', device = 'png', 
       width = 10, height = 5, bg = 'white')

# Generate summary figures for each species -------------------------------
for (i in 1:N) {
  print(i)
  curr.sp <- i
  plot.df <- data.frame(x = coords.0[, 1], 
  		      y = coords.0[, 2], 
  		      tmax.med = beta.star.tmax.quants[2, curr.sp, ],
  		      tmax.prob.pos = beta.star.tmax.prob.pos[curr.sp, ],
  		      tmax.ci.width = beta.star.tmax.quants[3, curr.sp, ] - 
  			              beta.star.tmax.quants[1, curr.sp, ], 
                        psi.med = psi.quants[2, curr.sp, ],
                        psi.ci.width = psi.quants[3, curr.sp, ] - 
                                       psi.quants[1, curr.sp, ])
  plot.df[range.pred.full[curr.sp, ] == 0, -c(1, 2)] <- NA
  
  plot.df <- st_as_stars(plot.df, dims = c('x', 'y'))
  
  tmax.plot <- ggplot() +
    geom_stars(data = plot.df, aes(x = x, y = y, fill = tmax.med), interpolate = TRUE) +	
    geom_sf(data = usa, alpha = 0, col = 'grey') +
    geom_sf(data = within.sf.list[[curr.sp]], alpha = 0, col = 'black') + 
    scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
    	               na.value = NA) + 
    theme_bw(base_size = 18) +
    labs(x = "Longitude", y = "Latitude", fill = "", 
         title = '(A) Median TMAX effect') +
    theme(legend.position = c(0.92, 0.25), 
          legend.background = element_rect(fill = NA), 
  	text = element_text(family = 'LM Roman 10'),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.title = element_text(size = 12),
  	plot.title = element_text(size = 18),
          legend.text = element_text(size = 12))
  
  tmax.prob.plot <- ggplot() +
    geom_stars(data = plot.df, aes(x = x, y = y, fill = tmax.prob.pos), interpolate = TRUE) +	
    geom_sf(data = usa, alpha = 0, col = 'grey') +
    geom_sf(data = within.sf.list[[curr.sp]], alpha = 0, col = 'black') + 
    scale_fill_steps2(midpoint = 0.5, low = '#B2182B', mid = 'white', high = '#2166AC', 
    	               na.value = NA, limits = c(0, 1), n.breaks = 6) + 
    theme_bw(base_size = 18) +
    labs(x = "Longitude", y = "Latitude", fill = "", 
         title = '(B) P(TMAX effect > 0)') +
    theme(legend.position = c(0.92, 0.25), 
          legend.background = element_rect(fill = NA), 
  	text = element_text(family = 'LM Roman 10'),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.title = element_text(size = 12),
  	plot.title = element_text(size = 18),
          legend.text = element_text(size = 12))
  
  psi.plot <- ggplot() +
    geom_stars(data = plot.df, aes(x = x, y = y, fill = psi.med), interpolate = TRUE) +	
    geom_sf(data = usa, alpha = 0, col = 'grey') +
    geom_sf(data = within.sf.list[[curr.sp]], alpha = 0, col = 'black') + 
    scale_fill_gradientn("", colors = ocean.tempo(5), limits = c(0, 1),
                             guide = guide_colourbar(title.position="top", reverse = FALSE), 
  			   na.value = NA) +
    theme_bw(base_size = 18) +
    labs(x = "Longitude", y = "Latitude", fill = "", 
         title = '(C) Median occurrence probability') +
    theme(legend.position = c(0.92, 0.25), 
          legend.background = element_rect(fill = NA), 
  	text = element_text(family = 'LM Roman 10'),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.title = element_text(size = 12),
  	plot.title = element_text(size = 18),
          legend.text = element_text(size = 12))
  psi.ci.plot <- ggplot() +
    geom_stars(data = plot.df, aes(x = x, y = y, fill = psi.ci.width), interpolate = TRUE) +	
    geom_sf(data = usa, alpha = 0, col = 'grey') +
    geom_sf(data = within.sf.list[[curr.sp]], alpha = 0, col = 'black') + 
    scale_fill_gradientn("", colors = ocean.tempo(5), limits = c(0, 1),
                             guide = guide_colourbar(title.position="top", reverse = FALSE), 
  			   na.value = NA) +
    theme_bw(base_size = 18) +
    labs(x = "Longitude", y = "Latitude", fill = "", 
         title = '(D) Occurrence 95% CI Width') +
    theme(legend.position = c(0.92, 0.25), 
          legend.background = element_rect(fill = NA), 
  	text = element_text(family = 'LM Roman 10'),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.title = element_text(size = 12),
  	plot.title = element_text(size = 18),
          legend.text = element_text(size = 12))
  # Generate Supplemental Figures S4-S22, as well as Figure 4 and 5 in the main text
  plot(ggarrange(tmax.plot, tmax.prob.plot, psi.plot, psi.ci.plot))
  ggsave(file = paste0('figures/species-figs/', sp.names[i], '-plot.png'), device = 'png', units = 'in', 
         width = 13, height = 9.5, bg = 'white')
}

# Generate map of BBS route locations and covariates ----------------------
# Point locations
point.plot <- ggplot(coords.sf) +
    geom_sf(size = 1) +
    geom_sf(data = usa, fill = NA, color=alpha("black", 0.75)) +
    theme_bw(base_size = 14) +
    labs(x = 'Longitude', y = 'Latitude', 
	 title = '(A) BBS Locations') +
  theme(text = element_text(family = 'LM Roman 10'),
	plot.title = element_text(size = 13))
# Combined with maps of covariates to form Figure 1 later in the script
plot.df <- data.frame(x = coords.0[, 1], 
		      y = coords.0[, 2], 
		      tmax = tmax.0, 
                      grass = grass.0)

plot.df <- st_as_stars(plot.df, dims = c('x', 'y'))
tmax.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = tmax), interpolate = TRUE) +	
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  scale_fill_viridis_c(na.value = NA) +
  theme_bw(base_size = 14) +
  labs(x = "Longitude", y = "Latitude", fill = "", 
       title = '(B) Maximum Temperature (degrees Celsius)') +
  theme(legend.position = c(0.92, 0.25), 
        legend.background = element_rect(fill = NA), 
	legend.key.size = unit(.4, 'cm'),
	text = element_text(family = 'LM Roman 10'),
	plot.title = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

grass.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = grass), interpolate = TRUE) +	
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  scale_fill_viridis_c(na.value = NA) +
  theme_bw(base_size = 14) +
  labs(x = "Longitude", y = "Latitude", fill = "", 
       title = '(C) Proportion of Grassland Area') +
  theme(legend.position = c(0.92, 0.25), 
        legend.background = element_rect(fill = NA), 
	legend.key.size = unit(.4, 'cm'),
	text = element_text(family = 'LM Roman 10'),
	plot.title = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
# Figure 1
ggarrange(point.plot, tmax.plot, grass.plot, ncol = 1)
ggsave(file = 'figures/Figure-1.png', device = 'png', height = 12, width = 6, units = 'in', 
       bg = 'white')

# Assess interaction between max temperature and grassland effect ---------
load('results/joint-svc-interaction/beta-samples.rda')

beta.int.samples <- beta.samples[, (ncol(beta.samples) - N + 1):ncol(beta.samples)]

beta.int.df <- data.frame(med = apply(beta.int.samples, 2, median),
			  low = apply(beta.int.samples, 2, quantile, 0.025),
			  high = apply(beta.int.samples, 2, quantile, 0.975),
			  prob.neg = apply(beta.int.samples, 2, function(a) mean (a < 0)), 
                          sp = sp.names)

plot.order <- sp.names[order(beta.int.df$prob.neg)]

# Generate Figure S3
beta.int.df %>%
  mutate(sp = factor(sp, levels = plot.order, ordered = TRUE)) %>%
    ggplot(aes(x = med, y = sp, fill = prob.neg)) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_segment(aes(x = low, y = sp, xend = high, yend = sp), 
		 lineend = 'butt', linewidth = 1, col = 'lightgray') +
    geom_point(size = 4, pch = 21) +
    scale_fill_gradient2(midpoint = 0.5, high = '#B2182B', mid = 'white', low = '#2166AC',
    	               na.value = NA) +
    theme_classic(base_size = 17) +
    labs(x = 'Interaction Effect Size',
	 y = 'Species', fill = 'P(effect < 0)') + 
    theme(text = element_text(family = 'LM Roman 10'))
  ggsave(file = 'figures/Figure-S3.png', width = 7, height = 7, bg = 'white')

