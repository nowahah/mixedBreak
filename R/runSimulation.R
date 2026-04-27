## Import useful scripts / packages ====
library(lmbreak)
library(segmented)

source("R/trueTraj.R")
source("R/noiseTraj.R")
source("R/normalize.R")


## Simulation parameters ====

# FIXED PARAMETERS
n.sim <- 2L # 500L nb of repetition for each scenario
times <- list("value" = 20, "sd" = 2) # time of measurement in simulated experimentss
outlier.prob <- 1 / 2 # individual probability of having exactly one outlier observation
na.prob <- 1 / 20 # measurement probability to be missing
n.trail <- -1L # nb of trailing observations (default is -1)
break.min.dist <- list(x = 45, y = 0) # min between-breakpoints distance
score.sd <- list("slope" = 1, "plateau" = .25)


# VARYING PARAMETERS
n.obs <- c(10, 15, 20, 25, 50, 100) # nb of observations in the datasets


## Seed management ====
# nb of scenarios: product of length of parameter.range
n.scenario <- length(n.obs)

set.seed(20041418)
nsimAll <- n.sim * n.scenario
allseeds <- sample.int(n = .Machine$integer.max, size = nsimAll, replace=FALSE) 


## Allocating memory ====

# scenarios dataset
scenarios <- data.frame(
  scenarID = 1:n.scenario,
  break.x1 = NA_real_,
  break.x1.sd = NA_real_,
  break.y1 = NA_real_,
  break.y1.sd = NA_real_,
  break.x2 = NA_real_,
  break.x2.sd = NA_real_,
  score.sd = NA # only for the slopes ? if plateau variance does not change
)

# true values datasets
true.parameters <- data.frame(
  scenarID = rep(1:n.scenario, times = n.sim*n.obs),
  simID = rep(1:nsimAll, times = rep(n.obs, each=n.sim)),
  patientID = sequence(rep(n.obs, each=n.sim)),
  break.x1 = NA_real_,
  break.x1.sd = NA_real_,
  break.y1 = NA_real_
)

# TODO estimates dataset
estimates.seglme <- data.frame(
  simID = rep(1:nsimAll, times = rep(n.obs, each=n.sim)),
  patientID = sequence(rep(n.obs, each=n.sim)),
  break.avg = NA_real_,
  break.x1 = NA_real_,
  break.x1.sd = NA_real_,
  break.y1 = NA_real_,
  breaks.CI95.low = NA_real_,
  breaks.CI95.up = NA_real_,
  intercept = NA_real_,
  error = NA
)
# TODO - update fields
estimates.segreg <- data.frame(
  simID = rep(1:nsimAll, times = rep(n.obs, each=n.sim)),
  patientID = sequence(rep(n.obs, each=n.sim)),
  break.x1 = NA_real_,
  break.x1.sd = NA_real_,
  break.y1 = NA_real_,
  error = NA_character_
)


## Running the simulation ====

## FIRST LOOP on sample size
start.time <- Sys.time()
sim.nb = 0 # iterator count
for (n.patients in n.obs){

  ## Scenario set up
  break.3 <- data.frame(
    pattern = c(1, 0, NA),
    bp.x = cumsum(c(0, 77, 78)),  # psi, time duration of each phase (from reproduced study)
    bp.y = c(0, 9.3, 9.3),        # height of breakpoints
    bp.x.sd = c(0, 38, 33),      # noise levels
    bp.y.sd = c(0, 1, 0)
  )
  true.traj <- trueTraj(times[["value"]], break.3[1:3])

  ## REP LOOP on repetitions
  for (jj in 1:n.sim){
    
    ## increment
    sim.nb = sim.nb + 1
    
    # # debug
    # if (sim.nb!=11) next
    
    ## 0. update 'progress bar' information
    if (sim.nb %% 1 == 0){ print(paste(sim.nb, "/", nsimAll)) }

    ## 1. set seed
    set.seed(allseeds[sim.nb])

    ## 2. simulate data and store the truth
    sim.data <- noiseTraj(
      true.traj, n.patients,
      score.sd, times.sd = times[["sd"]],
      breakpoints.sd = break.3[4:5],

      break.min.dist = break.min.dist,
      outlier.prob = outlier.prob,
      na.prob = na.prob,
      n.trail = n.trail
    )
    
    true.parameters[true.parameters$simID == sim.nb, c(3,4,6)] <- 
      sim.data$sim.gen.model[,c("ID", "break.x1", "break.y1")]
    true.parameters[true.parameters$simID == sim.nb, 5] <- NA
      

    ## 3. evaluate models
    #  3.1. segmented.lme
    tryCatch(
      {
        mod.lme <- lme(
          score ~ 0+time, 
          random = ~time|ID,
          data = sim.data$sim.dataset,
          na.action = na.omit,
          control = lmeControl(msMaxIter = 1000, msMaxEval = 1000)
        )
        
        mod.seg.lme <- segmented.lme(
          mod.lme, ~seg(0+time, est=c(1,0)),
          random = list(ID = pdDiag(~1+G0))
        )
        
      },
      error = function(e){ 
        print("Error during estimation of 'segmented.lme' model.")
        # browser()
        estimates.seglme[sim.nb, 'error'] <<- e$message # error assignement is incorrect
      }
    )
    
    #  3.2. segmented (with apply)
    tryCatch(
      {
        # mod.seg <- segreg(
        #   score ~  0 + seg(time, npsi=1, est=c(1,0)), 
        #   data = sim.data$sim.dataset
        # )
        
        mod.seg <- sim.data$sim.dataset %>% 
          group_by(ID) %>% 
          mutate(model = list(segreg(score ~  0 + seg(time, npsi = 1, est = c(1,0)), 
                                     data=pick(everything())))) %>%
          select(ID, model) %>%
          distinct()
        
      },
      error = function(e){ 
        print("Error during estimation of 'segreg' model.")
        # browser()
        estimates.segreg[sim.nb, 'error'] <<- paste(c(e$body, e$parent$message), 
                                                    collapse=" ==> ")
      }
    )
    

    ## 4. store results / errors
    if(!is.na(estimates.seglme[sim.nb, 'error'])){
      # store results
    }
    
  }
  # END REP LOOP

  # store scenario in the dataset


}
# END
ending.seed <- .Random.seed
end.time <- Sys.time()
print("Simulation study lasted for:")
print(end.time-start.time)
print("Average time/iteration:")
print((end.time-start.time)/nsimAll)
print("Estimated simulation study duration (in hours):")
print((end.time-start.time)/nsimAll*9500/60)
# Estimated 17.37 h of simulation