## Import useful scripts / packages ====
library(segmented)
library(dplyr)

source("R/trueTraj.R")
source("R/noiseTraj.R")
source("R/normalize.R")


## Simulation parameters ====

# FIXED PARAMETERS
n.sim <- 1L # 594L nb of repetition for each scenario
times <- list("value" = 20, "sd" = 2) # time of measurement in simulated experiments
na.prob <- 1 / 20 # measurement probability to be missing
n.trail <- -1L # nb of trailing observations (default is -1)
break.min.dist <- list(x = 45, y = 0) # min between-breakpoints distance
plateau.sd <- .25


# VARYING PARAMETERS
# value: range of parameters
# ref.ind: index of case-based value in vector of value
n.obs <- list(value = c(10, 15, 20, 25, 50, 100), ref.ind = 4) # nb of subjects in the datasets
score.sd.slopes <- list(value = c(0.5, 1, 2), ref.ind = 2) # level of noise during slopes
break.sd <- list(value = c(5, 15, 33), ref.ind = 2) # level of noise on break1 time coordinate
peak.duration.mean <- list(value = c(32, 78, 143), ref.ind = 2) # avg distance between break1 and end
outlier.prob <- list(value = c(0,.5,1), ref.ind = 2) # ind prob of having exactly one outlier 



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
  scenario.data <- scenario.data %>% 
    distinct() %>%
    mutate(across(-fac.var, as.numeric),
           scenarID = 1:n.scenario)
  
  return(scenario.data)
}

scenario.data <- one.at.a.time(
  n.obs = n.obs, 
  score.sd.slopes = score.sd.slopes, 
  break.sd = break.sd, 
  peak.duration.mean = peak.duration.mean,
  outlier.prob = outlier.prob
)
# mixed data frame integration ? 
# syntax with scenarID: integer if homogeneous scenar
# ID.group if mixed groups
scenario.data <- rbind(scenario.data, data.frame(
  fac.var = rep("mixed.group", 2),
  n.obs = rep(25, 2),
  score.sd.slopes = rep(1, 2),
  break.sd =c(15, 33),
  peak.duration.mean = c(32, 143),
  outlier.prob = rep(1/2, 1),
  scenarID = max(scenario.data$scenarID) + c(1.1, 1.2)
))

rm(n.obs, score.sd.slopes, outlier.prob, break.sd, peak.duration.mean)

## Handling mixed or homogeneous scenarios
# TODO - CHECK IF it works nicely integrated in the loop
# TODO - add mixed group label integration in the loop
generateData <- function(true.traj, n.patients, score.sd, times.sd, 
                         breakpoints.sd, break.min.dist, outlier.prob, na.prob,
                         n.trail, mixed = F){
  if (!mixed) {
    res <- noiseTraj(
      true.traj, n.patients,
      score.sd, times.sd,
      breakpoints.sd,
      
      break.min.dist,
      outlier.prob,
      na.prob,
      n.trail
    )
  } else if (mixed) {
    # two inhomogeneous groups in the data
    # browser()
    df1 <- noiseTraj(
      true.traj[[1]], n.patients[1], score.sd[[1]], times.sd[1], breakpoints.sd[[1]], 
      break.min.dist, outlier.prob[1], na.prob[1], n.trail[1]
    )
    df2 <- noiseTraj(
      true.traj[[2]], n.patients[2], score.sd[[2]], times.sd[2], breakpoints.sd[[2]], 
      break.min.dist, outlier.prob[2], na.prob[2], n.trail[2]
    )
    
    res <- list()
    res$sim.dataset <- rbind(df1$sim.dataset %>% mutate(group=1), 
                             df2$sim.dataset %>% mutate(group=2))
    res$sim.gen.model <- rbind(df1$sim.gen.model %>% mutate(group=1), 
                               df2$sim.gen.model %>% mutate(group=2))
  }
  
  return(res)
  
}


## Seed management ====
# nb of scenarios: product of length of parameter.range
n.scenario <- ceiling(max(scenario.data$scenarID)) # length(unique(floor(scenario.data$scenarID)))

set.seed(20041418)
nsimAll <- n.sim * n.scenario
allseeds <- sample.int(n = .Machine$integer.max, size = nsimAll, replace=FALSE) 


## Allocating memory ====

# true values datasets
true.parameters <- data.frame(
  scenarID = rep(scenario.data$scenarID, times = n.sim*as.numeric(scenario.data[,"n.obs"])),
  simID = rep(1:nsimAll, times = rep(as.numeric(scenario.data[,"n.obs"]), each = n.sim)),
  patientID = sequence(rep(as.numeric(scenario.data[,"n.obs"]), each = n.sim)),
  break.x1 = NA_real_,
  break.y1 = NA_real_
)

# estimates dataset - mixed model
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
# estimates dataset - plateau model
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
sim.nb <- 0 # iterator count
mixed.scenar.ii <- 0 # rows in scenario.data being part of mixed scenario (init)
for (ii in 1:nrow(scenario.data)){
  
  # skip to not iterate multiple times over mixed group setting
  if(ii %in% mixed.scenar.ii) next

  ## Scenario set up (taking mixed scenario into account)
  mixed.scenar <- scenario.data[ii, "n.obs"] != floor(scenario.data[ii, "n.obs"])
  if (!mixed.scenar){
    # homogeneous scenario setting
    n.patients <- scenario.data[ii, "n.obs"]
    score.sd.slopes <- scenario.data[ii, "score.sd.slopes"]
    break01.sd <- scenario.data[ii, "break.sd"]
    break12.sd <- scenario.data[ii, "break.sd"] # same value for both breaks sd
    break12.mean <- scenario.data[ii, "peak.duration.mean"]
    outlier.prob <- scenario.data[ii, "outlier.prob"]
    
    score.sd <- list("slope" = score.sd.slopes, "plateau" = plateau.sd)

    break.3 <- data.frame(
      pattern = c(1, 0, NA),
      bp.x = cumsum(c(0, 77, break12.mean)),  # psi, time duration of each phase (from reproduced study)
      bp.y = c(0, 9.3, 9.3),                  # height of breakpoints
      bp.x.sd = c(0, break01.sd, break12.sd), # noise levels
      bp.y.sd = c(0, 1, 0)
    )
    true.traj <- trueTraj(times[["value"]], break.3[1:3])
    
  } else if (mixed.scenar){
    # mixed scenario setting
    mixed.scenar.ii <- scenario.data$scenarID > floor(scenario.data[ii, "scenarID"]) & 
      scenario.data$scenarID < ceiling(scenario.data[ii, "scenarID"])
    mixed.data <- scenario.data[mixed.scenar.ii ,]
    
    n.patients <- scenario.data[mixed.scenar.ii, "n.obs"]
    score.sd.slopes <- scenario.data[mixed.scenar.ii, "score.sd.slopes"]
    break01.sd <- scenario.data[mixed.scenar.ii, "break.sd"]
    break12.sd <- scenario.data[mixed.scenar.ii, "break.sd"] # same value for both breaks sd
    break12.mean <- scenario.data[mixed.scenar.ii, "peak.duration.mean"]
    outlier.prob <- scenario.data[mixed.scenar.ii, "outlier.prob"]
    
    score.sd <- list(
      list("slope" = score.sd.slopes[1], "plateau" = plateau.sd),
      list("slope" = score.sd.slopes[2], "plateau" = plateau.sd)
    )
    
    g1.break <- data.frame(
      pattern = c(1, 0, NA),
      bp.x = cumsum(c(0, 77, break12.mean[1])), 
      bp.y = c(0, 9.3, 9.3),                  
      bp.x.sd = c(0, break01.sd[1], break12.sd[1]),
      bp.y.sd = c(0, 1, 0)
    )
    g1.true.traj <- trueTraj(times[["value"]], g1.break[1:3])
    g2.break <- data.frame(
      pattern = c(1, 0, NA),
      bp.x = cumsum(c(0, 77, break12.mean[2])), 
      bp.y = c(0, 9.3, 9.3),                  
      bp.x.sd = c(0, break01.sd[2], break12.sd[2]),
      bp.y.sd = c(0, 1, 0)
    )
    g2.true.traj <- trueTraj(times[["value"]], g2.break[1:3])
    true.traj <- list(g1.true.traj, g2.true.traj)
    break.3 <- list(g1.break, g2.break)
    
  }
  

  

  ## REP LOOP on repetitions
  for (jj in 1:n.sim){

    ## increment
    sim.nb = sim.nb + 1

    # # debug - go to a specific iteration
    # if (sim.nb!=11) next

    ## 0. update 'progress bar' information
    if (sim.nb %% 100 == 0 | sim.nb %in% c(2, nsimAll)){ print(paste(sim.nb, "/", nsimAll)) }

    ## 1. set seed
    set.seed(allseeds[sim.nb])

    ## 2. simulate data and store the truth
    browser()
    # different cases depending on mixture of population in scenario
    if(!is.data.frame(break.3)) {
      break.3 <- lapply(break.3, function(x) x[4:5])
    } else {
      
    }
    sim.data <- generateData(
      true.traj, n.patients,
      score.sd, times.sd = times[["sd"]],
      breakpoints.sd =,
      
      break.min.dist = break.min.dist,
      outlier.prob = outlier.prob,
      na.prob = na.prob,
      n.trail = n.trail,
      mixed = mixed.scenar
    )  
    
    if (!mixed) {
      sim.data <- noiseTraj(
        true.traj, n.patients,
        score.sd, times.sd = times[["sd"]],
        breakpoints.sd = break.3[4:5],
  
        break.min.dist = break.min.dist,
        outlier.prob = outlier.prob,
        na.prob = na.prob,
        n.trail = n.trail,
        mixed = is.list(n.obs)
      )    
      
      # store individual level truth
      true.parameters[
        true.parameters$simID == sim.nb, c("patientID", "break.x1", "break.y1")
      ] <- sim.data$sim.gen.model[,c("ID", "break.x1", "break.y1")]
    } else if (mixed) {
      
    }




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
        
        # browser()
        mod.seg.lme <- segmented.lme(
          mod.lme, ~time,
          random = list(ID = pdDiag(~1+G0))
        )
      },
      error = function(e){
        # print("Error during estimation of 'segmented.lme' model.")
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
        # print("Error during estimation of 'segreg' model.")
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
rm(n.patients, score.sd.slopes, break01.sd, break12.sd, break12.mean)
rm(res.segreg, res.seg.lme, mod.segreg, mod.lme, mod.seg.lme)

# END

print("Simulation study lasted for:")
print(end.time-start.time)
# Time difference of 9.825832 hours
# print("Average time/iteration:")
# print((end.time-start.time)/nsimAll)
# print("Estimated simulation study duration (in hours):")
# print((end.time-start.time)/nsimAll*n.scenario*559/60)
# Estimated ~8h of simulation


## WRITING RESULTS IN FILE
library(writexl)
write_xlsx(estimates.seglme, "data/esimates_seglme.xlsx")
write_xlsx(estimates.segreg, "data/esimates_segreg.xlsx")
write_xlsx(list(true.parameters=true.parameters, scenario.data = scenario.data), 
           "data/true_parameters.xlsx")
write_xlsx(list(all.seeds=data.frame(allseeds), ending.seed=data.frame(ending.seed)), 
                "data/random_seed.xlsx")