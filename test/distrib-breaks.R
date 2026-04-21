library(ggplot2)
library(MASS)


## Have a proper look  with the mean / variance parameters, if I go exp() on x but not y ?
# impact on correlation
# multivariate normal law distribution
# TODO - truncated ??
n.obs <- 10000L
peak <- 9
beta.onset <- 10/71 # 71 min if plateau at 10 SDI
mu = c(peak/beta.onset, peak) 
sd.x <- 33 # study value is 33, this one is for showing the correlation strategy only
sd.y <- 1
rho <- beta.onset*sd.x/sd.y # positive correlation condition
# TODO CHECK negative correlation shift
rho <- -beta.onset*sd.y/sd.x # correlation < 0 condition (orthogonal line, higher ratio variance)
rho <- -1 # arbitrary value
## condition for specifying correlation so ratio y/x has mean beta.target
if(rho > 1 | rho < -1){
  warning(paste0("Variance parameters forces rho outside of range [-1;1] (rho = ", round(rho,2), ")"))
  if(rho > 1) {
    rho <- 1
  }else if (rho < -1){
    rho <- -1
  }
  warning(paste("rho set to", rho, "instead of", round(beta.onset*sd.x/sd.y, 2)))
}

cor.mat <- matrix(c(sd.x^2, rho*sd.x*sd.y, rho*sd.x*sd.y, sd.y^2), ncol=2)
gen.norm <- as.data.frame(mvrnorm(n.obs, mu, cor.mat))
names(gen.norm) <- c("x", "y")

## breakpoint 1 density plot
ggplot(gen.norm, aes(x=x, y=y)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", colour="white") +
  
  # Adding reference slope
  geom_abline(slope = 1.4/10, colour = "black") +
  geom_abline(slope = rho*sd.y/sd.x, intercept = mu[2]-mu[1]*rho*sd.y/sd.x,
              linetype="dashed", color = "red") +
  
  # custom limits and ticks
  scale_x_continuous(breaks = seq(0, 120, by=20), limits = c(0, 120)) + # custom x ticks (SPS)
  scale_y_continuous(breaks = seq(0, 12, by=2), limits = c(6, 12)) + # custom y ticks
  
  # labels
  labs(title = paste0("Breakpoints' density plot (correlation = ", round(rho, 2),")"),
       x = "Time since drug intake (minutes)", y = "SDI scale")

# see https://en.wikipedia.org/wiki/Multivariate_normal_distribution
# "loci are squeezed toward the following line : 
# y(x) = sgn(rho)*sd.y/sd.x*(x-mu.x)+mu.y
# bcs expression with sgn(rho) replaced by rho is BLUP of y given x"
# CONCLUSIONS:
# Considering the expression, the correlation strategy could be successfully applied iff:
# sd.y/sd.x > beta.onset (easier to check and interpret) OR
# (mu.y/mu.x) * (sd.x/sd.y) \in [-1;1]
# which typically is not the case given observed simulation (sd.x=33, sd.y=??)
# Could possibly try and save the day by considering orthogonal line passing through
# the point with negative correlation


## onset slope density plot
ggplot(gen.norm, aes(x=y/x)) +
  geom_density() +
  
  # add normal law density on top
  stat_function(fun = function(u) dnorm(u, mean(gen.norm$y/gen.norm$x), sd(gen.norm$y/gen.norm$x)),
                color = "red", linetype = "dotted", size = 1) +
  xlim(.075,.2125) +
  
  # labels
  labs(title = paste("Onset slope density plot - (x,y) ~ MVN with correlation =", round(rho,2)),
       subtitle = paste("mean =", round(mean(gen.norm$y/gen.norm$x), 4),
                        "; sd =", round(sd(gen.norm$y/gen.norm$x), 4)),
       x = "Onset slope value", y = "Density value")


## TODO - Simulated multiple datasets to assess average value of ratio 
# (un)biased whenever we met (or not) the conditions on sd ratio

################################################################################

## Needs a trajData object in sim.data

## density plots (need enough observations to be beautiful - and meaningful)
# slope up VS down
ggplot(sim.data$sim.gen.model, aes(x=beta.1, y=beta.3)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", colour="white") +
  geom_point(alpha=0)
cor.test(sim.data$sim.gen.model$beta.1,sim.data$sim.gen.model$beta.3)

# library(ggExtra) 3 marginal density plots
# ggMarginal(p, type = "density", fill = "transparent") # not so good - p<-previous_graph

# break 1/2/3 (x, y) coordinates - one at a time
ggplot(sim.data$sim.gen.model, aes(x=break.x1, y=break.y1)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", colour="white") 


# all breakpoints (x, y) densities
ggplot(sim.data$sim.gen.model, aes(x=break.x1, y=break.y1), alpha=.5) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", colour="gray10") +
  stat_density_2d(aes(x=break.x2, y=break.y2, fill = after_stat(level), alpha=.5), 
                  geom = "polygon", colour="gray10") +
  stat_density_2d(aes(x=break.x3, y=break.y3, fill = after_stat(level), alpha=.5), 
                  geom = "polygon", colour="gray10") +
  
  # color customization
  scale_fill_distiller(palette = "Purples", direction = 1) +
  
  # custom limits and ticks
  scale_x_continuous(breaks = seq(0, 500, by=60), limits = c(0, 500)) + # custom x ticks for SPS
  scale_y_continuous(breaks = seq(0, 10, by=2), limits = c(0, 10)) + # custom y ticks
  
  # labels
  labs(title = "Breakpoints' density plot (independent coordinates)",
       x = "Time since drug intake (minutes)", y = "SDI scale")

# gradient scale is poorly adjusted due to heterogeneous variances


## OBSERVATIONS
# Here we can assess well the escalation in variation due to modelling
# the distances and summing independent r.v.