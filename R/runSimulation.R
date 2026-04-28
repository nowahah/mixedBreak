## Import useful scripts / packages ====
library(lmbreak)
library(segmented)
library(dplyr)

source("R/trueTraj.R")
source("R/noiseTraj.R")
source("R/normalize.R")


## Simulation parameters ====

# FIXED PARAMETERS
n.sim <- 559L # 559L nb of repetition for each scenario
times <- list("value" = 20, "sd" = 2) # time of measurement in simulated experiments
outlier.prob <- 1 / 2 # individual probability of having exactly one outlier observation
na.prob <- 1 / 20 # measurement probability to be missing
n.trail <- -1L # nb of trailing observations (default is -1)
break.min.dist <- list(x = 45, y = 0) # min between-breakpoints distance
plateau.sd <- .25


# VARYING PARAMETERS
# value: range of parameters
# ref.ind: index of case-based value in vector of value
n.obs <- list(value = c(10, 15, 20, 25, 50, 100), ref.ind = 4) # nb of observations in the datasets
score.sd.slopes <- list(value = c(0.5, 1, 2), ref.ind = 2) # level of noise during slopes
break01.sd <- list(value = c(5, 15, 33), ref.ind = 2) # level of noise on break1 time coordinate
break12.mean <- list(value = c(32, 78, 143), ref.ind = 2) # avg distance between break1 and end
break12.sd <- list(value = c(5, 15, 33), ref.ind = 2) # noise on distance between break1 and end


# SIMULATION GRID (varying parameters only - could add fixed ones also)
one.at.a.time <- function(...){
  factors <- list(...)
  n.scenario <- 1
  for (ii in seq_along(factors)){
    n.scenario = n.scenario + length(factors[[ii]]$value) - 1
  }
  
  # initialize scenario data
  scenario.data <- matrix(0, ncol = length(factors)+1, nrow = n.scenario)
  scenario.data <- data.frame(scenario.data)
  colnames(scenario.data) <- c("fac.var", names(factors))
  
  # complete scenario data
  for (ii in seq_along(factors)){
    fixed <- rep(NA, length(factors))
    fixed[-ii] <- unlist(lapply(factors, function(x) x$value[x$ref.ind]))[-ii]
    for (jj in seq_along(factors[[ii]]$value)){
      fixed[ii] <- factors[[ii]]$value[jj]
      fac.var <- if_else(jj==factors[[ii]]$ref.ind, "REFERENCE", names(factors)[[ii]])
      scenario.data[
        c(0, cumsum(unlist(lapply(factors, function(x) length(x$value)))))[ii] + jj,
        ] <- c(fac.var, fixed)
    }
  }
  
  # distinct() to not repeat reference scenario multiple times
  return(scenario.data %>% mutate(across(-fac.var, as.numeric)) %>% distinct())
}

scenario.data <- one.at.a.time(
  n.obs = n.obs, 
  score.sd.slopes = score.sd.slopes, 
  break01.sd = break01.sd, 
  break12.sd = break12.sd, 
  break12.mean = break12.mean
)
rm(n.obs, score.sd.slopes, break01.sd, break12.sd, break12.mean)


## Seed management ====
# nb of scenarios: product of length of parameter.range
n.scenario <- nrow(scenario.data)

set.seed(20041418)
nsimAll <- n.sim * n.scenario
allseeds <- sample.int(n = .Machine$integer.max, size = nsimAll, replace=FALSE) 


## Allocating memory ====

# true values datasets
true.parameters <- data.frame(
  scenarID = rep(1:n.scenario, times = n.sim*as.numeric(scenario.data[,"n.obs"])),
  simID = rep(1:nsimAll, times = rep(as.numeric(scenario.data[,"n.obs"]), each = n.sim)),
  patientID = sequence(rep(as.numeric(scenario.data[,"n.obs"]), each = n.sim)),
  break.x1 = NA_real_,
  break.y1 = NA_real_
)

# TODO estimates dataset
estimates.seglme <- data.frame(
  simID = true.parameters$simID,
  patientID = true.parameters$patientID,
  break.avg = NA_real_,
  break.x1 = NA_real_,
  break.x1.sd = NA_real_,
  break.x1.random = NA_real_,
  break.y1 = NA_real_,
  break.CI95.low = NA_real_,
  break.CI95.up = NA_real_,
  error = NA
)
# TODO - update fields
estimates.segreg <- data.frame(
  simID = true.parameters$simID,
  patientID = true.parameters$patientID,
  break.x1 = NA_real_,
  break.x1.sd = NA_real_,
  break.CI95.low = NA_real_,
  break.CI95.up = NA_real_,
  break.y1 = NA_real_,
  error = NA_character_
)


## Running the simulation ====

## FIRST LOOP on sample size
start.time <- Sys.time()
sim.nb = 0 # iterator count
for (ii in 1:nrow(scenario.data)){


  ## Scenario set up
  n.patients <- scenario.data[ii, "n.obs"]
  score.sd.slopes <- scenario.data[ii, "score.sd.slopes"]
  break01.sd <- scenario.data[ii, "break01.sd"]
  break12.sd <- scenario.data[ii, "break12.sd"]
  break12.mean <- scenario.data[ii, "break12.mean"]

  score.sd <- list("slope" = score.sd.slopes, "plateau" = plateau.sd)

  break.3 <- data.frame(
    pattern = c(1, 0, NA),
    bp.x = cumsum(c(0, 77, break12.mean)),  # psi, time duration of each phase (from reproduced study)
    bp.y = c(0, 9.3, 9.3),                  # height of breakpoints
    bp.x.sd = c(0, break01.sd, break12.sd), # noise levels
    bp.y.sd = c(0, 1, 0)
  )
  true.traj <- trueTraj(times[["value"]], break.3[1:3])

  ## REP LOOP on repetitions
  for (jj in 1:n.sim){

    ## increment
    sim.nb = sim.nb + 1

    # # debug - go to a specific iteration
    # if (sim.nb!=11) next

    ## 0. update 'progress bar' information
    if (sim.nb %% 100 == 0){ print(paste(sim.nb, "/", nsimAll)) }

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

    # store individual level truth
    true.parameters[
      true.parameters$simID == sim.nb, c("patientID", "break.x1", "break.y1")
    ] <- sim.data$sim.gen.model[,c("ID", "break.x1", "break.y1")]


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
        estimates.seglme[estimates.seglme$simID == sim.nb, 'error'] <<- e$message
      }
    )

    #  3.2. segmented (with apply)
    tryCatch(
      {
        mod.segreg <- sim.data$sim.dataset %>%
          group_by(ID) %>%
          mutate(model = list(segreg(score ~  0 + seg(time, npsi = 1, est = c(1,0)),
                                     data=pick(everything())))) %>%
          dplyr::select(ID, model) %>%
          ungroup() %>%
          distinct()

      },
      error = function(e){
        print("Error during estimation of 'segreg' model.")
        # browser()
        estimates.segreg[estimates.segreg$simID == sim.nb, 'error'] <<-
          paste(c(e$body, e$parent$message), collapse=" ==> ")
      }
    )


    ## 4. store results / errors
    # 4.1. segmented.lme results
    if(all(is.na(estimates.seglme[estimates.seglme$simID==sim.nb, 'error']))){

      res.seg.lme <- normalize(mod.seg.lme)
      # population level values
      estimates.seglme[
        estimates.seglme$simID==sim.nb,
        paste0("break.", c("avg", "x1.sd", "x1.random", "CI95.low", "CI95.up"))
      ] <- matrix(unlist(
        res.seg.lme[
          paste0("breaks.", c("avg", "sd", "random.sd","CI95.low", "CI95.up"))
        ], use.names = F
      ), nrow = n.patients, ncol = 5, byrow = T)
      # individual values
      estimates.seglme[
        estimates.seglme$simID==sim.nb,
        paste0("break.", c("x1", "y1"))
      ] <- unlist(res.seg.lme[c("breaks.ind", "intercept")], use.names = F)
    }

    # 4.2. segreg results
    if(all(is.na(estimates.segreg[estimates.segreg$simID==sim.nb, 'error']))){

      res.segreg <- normalize(mod.segreg)

      # only individual level values
      estimates.segreg[
        estimates.segreg$simID==sim.nb,
        paste0("break.", c("x1", "x1.sd", "CI95.low", "CI95.up", "y1"))
      ] <- res.segreg[,paste0("breaks.", c("x1", "x1.sd", "CI95.low", "CI95.up", "y1"))]
    }

  }
  # END REP LOOP

}
ending.seed <- .Random.seed
end.time <- Sys.time()
rm(ii, jj, sim.nb, break.3, score.sd, true.traj, sim.data)
rm(res.segreg, res.seg.lme, mod.segreg, mod.lme, mod.seg.lme)

# END

print("Simulation study lasted for:")
print(end.time-start.time)
# print("Average time/iteration:")
# print((end.time-start.time)/nsimAll)
# print("Estimated simulation study duration (in hours):")
# print((end.time-start.time)/nsimAll*n.scenario*559/60)
# Estimated ~8h of simulation


## WRITING RESULTS IN FILE
library(writexl)
write_xlsx(estimates.seglme, "data/esimates_seglme.xlsx")
write_xlsx(estimates.segreg, "data/esimates_segreg.xlsx")
write_xlsx(data.frame(ending.seed), "data/random_seed.xlsx")

