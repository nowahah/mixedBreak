## Import useful scripts / packages ====
library(lmbreak)
library(segmented)

source("R/trueTraj.R")
source("R/noiseTraj.R")


## Simulation parameters ====

# FIXED PARAMETERS
n.sim <- 500L # 500L nb of repetition for each scenario
times <- list("value" = 20, "sd" = 2) # time of measurement in simulated experimentss
outlier.prob <- 1 / 2 # individual probability of having exactly one outlier observation
na.prob <- 1 / 20 # measurement probability to be missing
n.trail <- -1L # nb of trailing observations (default is -1)
break.min.dist <- list(x = 45, y = 0) # min between-breakpoints distance
score.sd <- list("slope" = 1, "plateau" = .25)


# VARYING PARAMETERS
n.obs <- c(10, 15, 20, 25, 50, 100) # nb of observations in the dataset


## Seed management ====
# nb of scenarios: product of length of parameter.range
n.scenario <- length(n.obs)

set.seed(20041418)
nsimAll <- n.sim * n.scenario
allseeds <- sample.int(n = .Machine$integer.max, size = nsimAll, replace=FALSE) 


## Allocating memory ====



# ## Running the simulation ====
# 
# ## FIRST LOOP on sample size
# start.time <- Sys.time()
# sim.nb = 1
# for (n.patients in n.obs){
#   
#   
#   ## Scenario set up
#   break.3 <- data.frame(
#     pattern = c(1, 0, NA),
#     bp.x = cumsum(c(0, 77, 78)),  # psi, time duration of each phase (from reproduced study)
#     bp.y = c(0, 9.3, 9.3),        # height of breakpoints
#     bp.x.sd = c(0, 38, 33),      # noise levels
#     bp.y.sd = c(0, 1, 0)
#   )
#   true.traj <- trueTraj(times[["value"]], break.3[1:3])
# 
#   ## REP LOOP on repetitions
#   for (j in 1:n.sim){
#     
#     # 0. update progress bar information
#     if (sim.nb %% 10 == 0){ print(paste(sim.nb, "/", nsimAll)) }
#     
#     # 1. set seed
#     set.seed(allseeds[sim.nb])
#     
#     # 2. simulate data
#     sim.data <- noiseTraj(
#       true.traj, n.patients, 
#       score.sd, times.sd = times[["sd"]],
#       breakpoints.sd = break.3[4:5], 
#       
#       break.min.dist = break.min.dist,
#       outlier.prob = outlier.prob, 
#       na.prob = na.prob, 
#       n.trail = n.trail
#     )
#     
#     # 3. evaluate models
#     try(T)
#     
#     # 4. store results / errors
#     
#     
#     # increment
#     sim.nb = sim.nb + 1
#   }
#   # END REP LOOP
#   
#   # store scenario in the dataset
#   
#   
# }
# # END
# end.time <- Sys.time()
# print("Simulation study lasted for:")
# end.time-start.time
