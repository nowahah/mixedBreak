### simData-test.R ----

source("R/DataLoader.R")

source("R/simData.R")

library(ggplot2)

n.obs <- 4L
breakpoints <- data.frame(
  bp = c(90, 180),
  bp.sd = c(10, 10)
)
b.onset <- list("value" = 2*1.4, "sd" = 2*0.5, "min" = 0.5) # onset coefficient per 20 min
b.return <- list("value" = 2*-0.7, "sd" = 2*0.3, "max" = -0.2) # return to normal coefficient per 20 min
score.sd <- 1 # noise level on the measurements
# different time specifications
time.reg <- list("value" = 20, "sd" = 0) # default regular "perfect" measurement strategy
time.noise <- list("value" = 20, "sd" = 2) # regular "perfect" measurement strategy with noise
time.spec <- list("value" = c(0,10,20,30,40,60,120), "sd" = 1)   # specific measurement strategy

outlier.prob <- 1/5
na.prob <- 0.10
n.trail <- 0L


sim.dataset <- simData(n.obs, breakpoints, b.onset, b.return, score.sd,
                       times = time.noise,
                       outlier.prob = outlier.prob, na.prob = na.prob,
                       n.trail = n.trail)
plot(sim.dataset)  # default trajectory plot
print(sim.dataset) # data generation model

# print(sim.dataset, default = T) # default method for data.frame
print(sim.dataset, n.line = 10L, dgm = F) # dataset of observations
(sim.measurements <- subset(sim.dataset))

# Higher variance during plateau (especially if peak <≈ 10) -> not realistic (need lower variance during plateau)
