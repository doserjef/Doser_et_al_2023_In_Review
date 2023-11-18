# summary-sim.R: this script summarizes results from the simulation study
#                in which we perform a proof of concept for the single-species SVC
#                occupancy model across varying types of spatial variability in the 
#                covariate effect. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(viridis)
library(devtools)
library(ggthemes)
library(spOccupancy)

# SVC occupancy model results ---------------------------------------------
# Load results ------------------------------------------------------------
load("results/sim-results.rda")
# Number of sites
J <- dim(psi.mean.samples)[1]
# Number of simulations for each scenario/model combo
n.sims <- dim(psi.mean.samples)[2]
# Number of scenarios
n.scenarios <- dim(psi.mean.samples)[3]
# Number of models
n.models <- dim(psi.mean.samples)[4]

# Simulated spatial variance parameters
sigma.sq.vals <- c(0.1, 0.5, 1, 2)
phi.vals <- c(3 / 0.1, 3 / 0.3, 3 / 0.8, 3 / 3, 3 / 100)
scenario.vals <- expand.grid(sigma.sq = sigma.sq.vals, phi = phi.vals)

# Assessment of model deviance --------------------------------------------
# Average deviance score across all simulations
deviance.avg <- apply(deviance.vals, c(2, 3), mean)
deviance.low <- apply(deviance.vals, c(2, 3), quantile, 0.025)
deviance.high <- apply(deviance.vals, c(2, 3), quantile, 0.975)
plot.deviance <- data.frame(phi = factor(rep(c(scenario.vals$phi, 0), times = n.models)),
			    sigma.sq = factor(rep(c(scenario.vals$sigma.sq, 0), times = n.models)),
			    # deviance = c(deviance.avg) - deviance.avg[, 1],
			    deviance = c(deviance.avg),
			    type = factor(rep(c('SVC', 'SVI', 'OCC'), each = n.scenarios), 
			                  levels = c('OCC', 'SVI', 'SVC'), ordered = TRUE))
# levels(plot.deviance$phi) <- c('High Range', 'Medium Range', 'Low Range')
# levels(plot.deviance$sigma.sq) <- c('Low Variance', 'Medium Low Variance', 'Medium High Variance', 'High Variance')

# Deviance results from simulation study.
plot.deviance %>%
  arrange(desc(phi), sigma.sq, type)

# Assessment of WAIC ------------------------------------------------------
# Average WAIC score across all simulations
# waic.diff <- array(NA, dim = dim(waic.vals))
# for (i in 1:nrow(waic.vals)) {
#   for (j in 1:ncol(waic.vals)) {
#     waic.diff[i, j, ] <- waic.vals[i, j, ] - waic.vals[i, j, 1]
#   }
# }
waic.avg <- apply(waic.vals, c(2, 3), mean)
waic.low <- apply(waic.vals, c(2, 3), quantile, 0.025)
waic.high <- apply(waic.vals, c(2, 3), quantile, 0.975)
plot.waic <- data.frame(phi = factor(rep(c(scenario.vals$phi, 0), times = n.models)),
			    sigma.sq = factor(rep(c(scenario.vals$sigma.sq, 0), times = n.models)),
			    waic = c(waic.avg),
			    type = factor(rep(c('SVC', 'SVI', 'OCC'), each = n.scenarios), 
			                  levels = c('OCC', 'SVI', 'SVC'), ordered = TRUE))
# levels(plot.waic$phi) <- c('High Range', 'Medium Range', 'Low Range')
# levels(plot.waic$sigma.sq) <- c('Low Variance', 'Medium Low Variance', 'Medium High Variance', 'High Variance')
# 
# WAIC results from simulation study.
plot.waic %>%
  arrange(desc(phi), sigma.sq, type)

# Create overall figure to summarize simulations --------------------------
# Get rid of the no SVC case for the figure
plot.deviance.small <- plot.deviance %>%
  filter(phi != 0)
plot.waic.small <- plot.waic %>%
  filter(phi != 0)
plot.df <- left_join(plot.deviance.small, plot.waic.small, by = c('phi', 'sigma.sq', 'type'))
plot.long.df <- data.frame(val = c(plot.df$deviance, plot.df$waic),
			   phi = rep(plot.df$phi, 2),
			   sigma.sq = rep(plot.df$sigma.sq, 2),
			   type = rep(plot.df$type, 2),
			   criterion = rep(c('Deviance', 'WAIC'), each = nrow(plot.df)))
plot.long.df <- plot.long.df %>%
  mutate(phi = factor(rep(scenario.vals$phi, n.models * 2),
		      levels = c(30, 10, 3.75, 1, 0.03), ordered = TRUE, 
		#       labels = c('ESR = 10%', 'ESR = 30%', 'ESR = 80%', 
		# 		 'ESR = 300%', 'ESR = 10,000%')),
		      labels = c(expression(paste("ESR", " = 10%")), 
				 expression(paste("ESR", " = 30%")), 
				 expression(paste("ESR", " = 80%")), 
				 expression(paste("ESR", " = 300%")), 
				 expression(paste("ESR", " = 10,000%")))), 
		#       labels = c(expression(paste(phi, " = 3 / .1")), 
		# 		 expression(paste(phi, " = 3 / .3")), 
		# 		 expression(paste(phi, " = 3 / .8")), 
		# 		 expression(paste(phi, " = 3 / 3")), 
		# 		 expression(paste(phi, " = 3 / 100")))), 
	 sigma.sq = factor(rep(scenario.vals$sigma.sq, n.models * 2), 
			   labels = c(expression(paste(sigma, " "^2, " = 0.1")), 
				      expression(paste(sigma, " "^2, " = 0.5")),
				      expression(paste(sigma, " "^2, " = 1.0")),
				      expression(paste(sigma, " "^2, " = 2.0"))))) 

ggplot(data = plot.long.df, aes(x = type, y = val, col = criterion, group = criterion)) + 
  geom_point() + 
  geom_line() +
  scale_color_manual(values = c('#B2182B', '#2166AC')) + 
  facet_grid(phi ~ sigma.sq, 
             labeller = label_parsed) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        text = element_text(family = 'LM Roman 10'), 
	legend.position = 'bottom') +
  labs(x = 'Model', y = 'Value', col = 'Criterion')
ggsave(file = 'figures/Figure-1.png', device = 'png', units = 'in', width= 8, height = 7)

# Figure of spatial range and spatial variance ----------------------------
# This code creates Figure S1, which gives an example of how the spatial
# variance and spatial decay influence the resulting spatially-varying
# coefficient under the different conditions. 
# Spatial locations
J.x <- 50
J.y <- 50
J <- J.x * J.y
# Replicates
# n.rep <- sample(5, J, replace = TRUE)
n.rep <- rep(5, J)
# Occurrence coefficients -------------
beta <- c(0, 0)
# Detection coefficients --------------
alpha <- c(-0.2, 0.4)
# Spatial parameters ------------------
sp <- TRUE
svc.cols <- c(1, 2)
cov.model <- 'exponential'
# Values for spatially-varying intercept
sigma.sq.int <- 1
phi.int <- 3 / 0.4
plot.data <- list()
for (i in 1:nrow(scenario.vals)) {
  set.seed(100)
  print(i)
  dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
                sp = sp, svc.cols = svc.cols, cov.model = cov.model,
                sigma.sq = c(sigma.sq.int, scenario.vals$sigma.sq[i]),
                phi = c(phi.int, scenario.vals$phi[i]))
  
  plot.data[[i]] <- data.frame(x = dat$coords[, 1], 
  		               y = dat$coords[, 2], 
  		               w = dat$w[, 2])
}
plot.data.df <- do.call(rbind, plot.data)
plot.data.df <- plot.data.df %>%
  mutate(phi = factor(rep(scenario.vals$phi, each = J), 
		      levels = c(30, 10, 3.75, 1, 0.03), ordered = TRUE, 
		      labels = c(expression(paste(phi, " = 3 / .1")), 
				 expression(paste(phi, " = 3 / .3")), 
				 expression(paste(phi, " = 3 / .8")), 
				 expression(paste(phi, " = 3 / 3")), 
				 expression(paste(phi, " = 3 / 100")))), 
	 sigma.sq = factor(rep(scenario.vals$sigma.sq, each = J), 
			   labels = c(expression(paste(sigma, " "^2, " = 0.1")), 
				      expression(paste(sigma, " "^2, " = 0.5")),
				      expression(paste(sigma, " "^2, " = 1.0")),
				      expression(paste(sigma, " "^2, " = 2.0"))))) 
# Generate Figure 1
my.plot <- ggplot(data = plot.data.df, aes(x = x, y = y, fill = w)) + 
    scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
      	               na.value = NA) + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    geom_raster() + 
    facet_grid(phi ~ sigma.sq, 
	       labeller = label_parsed) + 
    theme_light(base_size = 14) + 
    labs(x = 'Easting', y = 'Northing', fill = 'Effect') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(), 
	  text = element_text(family = 'LM Roman 10'),
          axis.ticks.y = element_blank(), 
          strip.text.x = element_text(color = 'black'), 
          strip.text.y = element_text(color = 'black')) 
my.plot
ggsave(file = 'figures/Figure-S1.png', device = 'png', units = 'in', width= 8, height = 6)
