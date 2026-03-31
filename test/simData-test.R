### simData-test.R ----

source("R/DataLoader.R")

source("R/simData.R")

library(ggplot2)

n.obs <- 4L
breakpoints <- data.frame(
  bp = c(90, 180),
  bp.sd = c(10, 10)
)
b.onset <- list("value" = 2 * 1.4, "sd" = 2 * 0.5, "min" = 0.5) # onset coefficient per 20 min
b.return <- list("value" = 2 * -0.7, "sd" = 2 * 0.3, "max" = -0.2) # return to normal coefficient per 20 min
score.sd <- 1 # noise level on the measurements
# different time specifications
time.reg <- list("value" = 20, "sd" = 0) # default regular "perfect" measurement strategy
time.noise <- list("value" = 20, "sd" = 2) # regular "perfect" measurement strategy with noise
time.spec <- list("value" = c(0, 10, 20, 30, 40, 60, 120), "sd" = 1) # specific measurement strategy

outlier.prob <- 1 / 2 # individual probability of having exactly one outlier observation
na.prob <- 1 / 10 # measurement probability to be missing
n.trail <- 3L # nb of trailing observations


sim.dataset <- simData(
  n.obs = n.obs, score.sd = score.sd, times = time.noise,
  b.onset = b.onset,
  b.return = b.return,
  plateau = NULL,
  breakpoints = breakpoints,
  ending.times = NULL,
  outlier.prob = outlier.prob, na.prob = na.prob,
  n.trail = n.trail
)
plot(sim.dataset) # default trajectory plot
print(sim.dataset) # data generation model

# print(sim.dataset, default = T) # default method for data.frame
print(sim.dataset, n.line = 10L, dgm = F) # dataset of observations
(sim.measurements <- subset(sim.dataset))
# Higher variance during plateau (especially if peak <≈ 10) -> not realistic (need lower variance during plateau)
