library(ggplot2)
library(MASS)


## Have a proper look  with the mean / variance parameters, if I go exp() on x but not y ?
# impact on correlation
# multivariate normal law distribution
# TODO - truncated ??
n.obs <- 10000L
peak <- 9
mu = c(71*peak/10, peak) # 71 min if plateau at 10 SDI (ie beta.onset = 10/71 = 0.141)
sd.x <- 20
sd.y <- 3
rho <- mu[2]/mu[1]*sd.x/sd.y
# rho <- beta.onset * sd.x/sd.y
# rho <- .75 # arbitrary value
## condition for specifying correlation so ratio y/x has mean beta.target
if(rho > 1 | rho < -1){
  stop(paste0("Variance parameters forces rho outside of range [-1;1] (rho = ", round(rho,2), ")"))
}

cor.mat <- matrix(c(sd.x^2, rho*sd.x*sd.y, rho*sd.x*sd.y, sd.y^2), ncol=2)
gen.norm <- as.data.frame(mvrnorm(n.obs, mu, cor.mat))
names(gen.norm) <- c("x", "y")

## breakpoint 1 density plot
ggplot(gen.norm, aes(x=x, y=y)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", colour="white") +
  
  # Adding reference slope
  geom_abline(slope = 1.4/10) +
  geom_abline(slope = rho*sd.y/sd.x, intercept = mu[2]-mu[1]*rho*sd.y/sd.x,
              linetype="dashed", color = "red") +
  
  # custom limits and ticks
  scale_x_continuous(breaks = seq(0, 100, by=20), limits = c(0, 100)) + # custom x ticks (SPS)
  scale_y_continuous(breaks = seq(0, 12, by=2), limits = c(0, 12)) + # custom y ticks
  
  # labels
  labs(title = "Breakpoints' density plot (independent coordinates)",
       x = "Time since drug intake (minutes)", y = "SDI scale")

# see https://en.wikipedia.org/wiki/Multivariate_normal_distribution
# "loci are squeezed toward the following line : 
# y(x) = sgn(rho)*sd.y/sd.x*(x-mu.x)+mu.y
# bcs expression with sgn(rho) replaced by rho is BLUP of y given x"
# CONCLUSIONS:
# Considering the expression, the correlation strategy could be successfully applied iff:
# sd.y/sd.x > beta.onset OR
# (mu.y/mu.x) * (sd.x/sd.y) \in [-1;1]
# which typically is not the case given observed simulation (sd.x=33, sd.y=??)


## onset slope density plot
ggplot(gen.norm, aes(x=y/x)) +
  geom_density() +
  
  # add normal law density on top
  stat_function(fun = function(u) dnorm(u, mean(gen.norm$y/gen.norm$x), sd(gen.norm$y/gen.norm$x)),
                color = "red", linetype = "dotted", size = 1) +
  xlim(.075,.2125) +
  
  # labels
  labs(title = paste("Onset slope density plot - (x,y) ~ MVN with correlation =", round(rho,2)),
       x = "Onset slope value", y = "Density value")


## TODO - Relationship between correlation and beta values BLUP of y given x
# y(x) = rho*(sd.y/sd.x)*(x-mu.x) + mu.y