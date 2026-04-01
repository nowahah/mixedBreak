### simData10-test.R ----

source("R/DataLoader.R")

source("R/simData10.R")

library(ggplot2)

n.obs <- 6L # nb of patients
breakpoints <- data.frame(
  pattern = c(1, 0, 1, 0, NA),
  bp.x = c(0, 90, 180, 250, 360),    # psi, time coordinate of breakpoints
  bp.y = c(0, 9.5, 9.5, 2, 2),   # height of breakpoints
  bp.x.sd = c(0, 100, 100, 100, 100),  # noise levels
  bp.y.sd = c(0, 1, 0, 1, 0)
)

# different time specifications
time.reg <- list("value" = 20, "sd" = 0) # default regular "perfect" measurement strategy
time.noise <- list("value" = 20, "sd" = 2) # regular "perfect" measurement strategy with noise
time.spec <- list("value" = c(0, 10, 20, 30, 40, 60, 120), "sd" = 1) # specific measurement strategy

outlier.prob <- 1 / 2 # individual probability of having exactly one outlier observation
na.prob <- 1 / 10 # measurement probability to be missing
n.trail <- 0L # nb of trailing observations


sim.data <- simData10(
  n.obs = n.obs, score.sd = score.sd, times = time.noise,
  breakpoints = breakpoints,
  outlier.prob = outlier.prob,
  na.prob = na.prob,
  n.trail = n.trail
)

sim.dataset <- sim.data$sim.dataset
sim.gen.model <- sim.data$sim.gen.model

# TODO - need to rewrite plot, print, etc. according to new data structure
plot(sim.dataset) # default trajectory plot
print(sim.dataset) # data generation model

# print(sim.dataset, default = T) # default method for data.frame
print(sim.dataset, n.line = 10L, dgm = F) # dataset of observations
(sim.measurements <- subset(sim.dataset))
# Higher variance during plateau (especially if peak <≈ 10) -> not realistic (need lower variance during plateau)
