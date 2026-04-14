### simData10-test.R ----

# source("R/DataLoader.R")

# source("R/simData10.R")
source("R/trueTraj.R")
source("R/noiseTraj.R")
source("R/plot.R")
source("R/print.R")


# different time specifications are possible
time.reg <- list("value" = 20, "sd" = 0) # default regular "perfect" measurement strategy
time.noise <- list("value" = 20, "sd" = 2) # regular "perfect" measurement strategy with noise
time.spec <- list("value" = c(0, 10, 20, 30, 40, 60, 120), "sd" = 1) # specific measurement strategy

outlier.prob <- 1 / 2     # ind probability of having exactly one outlier 
na.prob <- 1 / 10         # measurement probability to be missing
n.trail <- -1L            # nb of trailing observations
break.min.dist <- list(x = 45, y = 0) # min between-breakpoints distance

n.obs <- 6L # nb of patients
score.sd <- .5
break.5 <- data.frame(
  pattern = c(1, 0, 1, 0, NA),
  bp.x = c(0, 90, 180, 250, 360),      # psi, time coordinate of breakpoints
  bp.y = c(0, 9.5, 9.5, 1, 1),         # height of breakpoints
  bp.x.sd = c(0, 100, 100, 100, 100),  # noise levels
  bp.y.sd = c(0, 1, 0, 1, 1)
)

break.4 <- data.frame(
  pattern = c(1, 0, 1, NA),
  bp.x = c(0, 90, 180, 300),       # psi, time coordinate of breakpoints
  bp.y = c(0, 9.5, 9.5, 1),        # height of breakpoints
  bp.x.sd = c(0, 50, 50, 50),      # noise levels
  bp.y.sd = c(0, 1, 0, 1)
)

break.3 <- data.frame(
  pattern = c(1, 1, NA),
  bp.x = c(0, 180, 300),    # psi, time coordinate of breakpoints
  bp.y = c(0, 9.5, 1),      # height of breakpoints
  bp.x.sd = c(0, 50, 50),   # noise levels
  bp.y.sd = c(0, 1, 1)
)

## TEST
# # simData10.R
# sim.data <- simData10(
#   n.obs = n.obs, score.sd = score.sd, times = time.noise,
#   breakpoints = break.4,
#   break.min.dist = break.min.dist,
#   outlier.prob = outlier.prob,
#   na.prob = na.prob,
#   n.trail = n.trail
# ) 

# splitted trueTraj.R & noiseTraj.R
times <- time.noise
breakpoints <- break.4
true.traj <- trueTraj(times[["value"]], breakpoints[1:3])
plot(true.traj)

sim.data <- noiseTraj(true.traj, n.obs, score.sd, times.sd = times[["sd"]],
                      breakpoints.sd = breakpoints[4:5], break.x.dist = 20,
                      outlier.prob = outlier.prob, na.prob = na.prob, n.trail = 1L)
# sim.dataset <- sim.data$sim.dataset
# sim.gen.model <- sim.data$sim.gen.model

## Default methods for trajData objects
# plot
plot(sim.data) # default trajectory plot
plot(sim.data, cluster = 1:4)
plot(sim.data, cluster = 1:4, breakpoints = F)
plot(sim.data, cluster = 1:4, lines = F)
plot(sim.data, cluster = 1:4, alpha = .5, true.color = "orange2")

# print
print(sim.data) # data and data generation model
print(sim.data, cluster = c(2,5))
print(sim.data, cluster = c(2,5), n.lines = 3L) # choose nb of lines of data to print
print(sim.data, cluster = c(2,5), dgm = F) # choose not to print data generation model


# ## slopes distribution density plot (need enough observations to be beautiful - and meaningful)
# ggplot(sim.data$sim.gen.model, aes(x=beta.1, y=beta.3)) +
#   stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", colour="white")
