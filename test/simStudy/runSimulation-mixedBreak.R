# ## Import useful scripts / packages ====
# library(segmented)
# library(dplyr)
# library(rlang) # Error in names(res[[1]]) %||% "" : could not find function "%||%"
# 
# source("R/trueTraj.R")
# source("R/noiseTraj.R")
# source("R/normalize.R")
# source("R/confint.R")
# source("R/print.R")
# source("R/plot.R")
# source("R/mixed.break1.R")
# source("R/confint.R")
# 
# 
# ## Simulation parameters ====
# 
# # FIXED PARAMETERS
# n.sim <- 594L # 594L nb of repetition for each scenario
# print(paste("n.sim =", n.sim, "repetitions per DGM"))
# times <- list("value" = 20, "sd" = 2) # time of measurement in simulated experiments
# na.prob <- 1 / 20 # measurement probability to be missing
# n.trail <- -1L # nb of trailing observations (default is -1)
# break.min.dist <- list(x = 45, y = 0) # min between-breakpoints distance
# plateau.sd <- .25
# 
# 
# # VARYING PARAMETERS
# # value: range of parameters
# # ref.ind: index of case-based value in vector of value
# n.obs <- list(value = c(10, 15, 20, 25, 50, 100), ref.ind = 4) # nb of subjects in the datasets
# score.sd.slopes <- list(value = c(0.5, 1, 2), ref.ind = 2) # level of noise during slopes
# break.sd <- list(value = c(5, 15, 33), ref.ind = 2) # level of noise on break1 time coordinate
# peak.duration.mean <- list(value = c(32, 78, 143), ref.ind = 2) # avg distance between break1 and end
# outlier.prob <- list(value = c(0,.5,1), ref.ind = 2) # ind prob of having exactly one outlier 
# 
# 
# 
# # SIMULATION GRID (varying parameters only - could add fixed ones also)
# one.at.a.time <- function(...){
#   factors <- list(...)
#   n.scenario <- 1
#   for (ii in seq_along(factors)){
#     n.scenario = n.scenario + length(factors[[ii]]$value) - 1
#   }
#   
#   # initialize scenario data
#   scenario.data <- matrix(0, ncol = length(factors)+1, nrow = n.scenario)
#   scenario.data <- data.frame(scenario.data)
#   colnames(scenario.data) <- c("fac.var", names(factors))
#   
#   # complete scenario data
#   for (ii in seq_along(factors)){
#     fixed <- rep(NA, length(factors))
#     fixed[-ii] <- unlist(lapply(factors, function(x) x$value[x$ref.ind]))[-ii]
#     for (jj in seq_along(factors[[ii]]$value)){
#       fixed[ii] <- factors[[ii]]$value[jj]
#       fac.var <- if_else(jj==factors[[ii]]$ref.ind, "REFERENCE", names(factors)[[ii]])
#       scenario.data[
#         c(0, cumsum(unlist(lapply(factors, function(x) length(x$value)))))[ii] + jj,
#       ] <- c(fac.var, fixed)
#     }
#   }
#   
#   # distinct() to not repeat reference scenario multiple times
#   scenario.data <- scenario.data %>% 
#     distinct() %>%
#     mutate(across(-fac.var, as.numeric),
#            scenarID = 1:n.scenario)
#   
#   return(scenario.data)
# }
# 
# scenario.data <- one.at.a.time(
#   n.obs = n.obs, 
#   score.sd.slopes = score.sd.slopes, 
#   break.sd = break.sd, 
#   peak.duration.mean = peak.duration.mean,
#   outlier.prob = outlier.prob
# )
# rm(n.obs, score.sd.slopes, outlier.prob, break.sd, peak.duration.mean)
# 
# 
# 
# ## Seed management ====
# # nb of scenarios: product of length of parameter.range
# n.scenario <- ceiling(max(scenario.data$scenarID)) # length(unique(floor(scenario.data$scenarID)))
# 
# set.seed(20041418)
# nsimAll <- n.sim * n.scenario
# allseeds <- sample.int(n = .Machine$integer.max, size = nsimAll, replace=FALSE) 
# 
# 
# ## Allocating memory ====
# 
# # true values datasets
# true.parameters <- data.frame(
#   scenarID = rep(scenario.data$scenarID, times = n.sim*as.numeric(scenario.data[,"n.obs"])),
#   simID = rep(1:nsimAll, times = rep(as.numeric(scenario.data[,"n.obs"]), each = n.sim)),
#   patientID = sequence(rep(as.numeric(scenario.data[,"n.obs"]), each = n.sim)),
#   break.x1 = NA_real_,
#   break.y1 = NA_real_
# )
# 
# # estimates dataset - mixed model
# estimates.11 <- data.frame(
#   simID = true.parameters$simID,
#   patientID = true.parameters$patientID,
#   break.x1.avg = NA_real_,
#   break.x1.sd = NA_real_,
#   break.x1.ind = NA_real_,
#   break.x1.random.sd = NA_real_,
#   break.y1.ind = NA_real_,
#   break.CI95.low = NA_real_,
#   break.CI95.up = NA_real_,
#   slope2 = NA_real_,
#   slope2.sd = NA_real_,
#   error = NA
# )
# # estimates dataset - plateau model
# estimates.10 <- data.frame(
#   simID = true.parameters$simID,
#   patientID = true.parameters$patientID,
#   break.x1.avg = NA_real_,
#   break.x1.sd = NA_real_,
#   break.x1.ind = NA_real_,
#   break.x1.random.sd = NA_real_,
#   break.y1.ind = NA_real_,
#   break.CI95.low = NA_real_,
#   break.CI95.up = NA_real_,
#   error = NA_character_
# )
# 
# 
# ## Running the simulation ====
# 
# ## FIRST LOOP on sample size
# start.time <- Sys.time()
# sim.nb <- 0 # iterator count
# for (ii in 1:nrow(scenario.data)){
#   # browser()
#   
#   # homogeneous scenario setting
#   n.patients <- scenario.data[ii, "n.obs"]
#   score.sd.slopes <- scenario.data[ii, "score.sd.slopes"]
#   break01.sd <- scenario.data[ii, "break.sd"]
#   break12.sd <- scenario.data[ii, "break.sd"] # same value for both breaks sd
#   break12.mean <- scenario.data[ii, "peak.duration.mean"]
#   outlier.prob <- scenario.data[ii, "outlier.prob"]
#   
#   score.sd <- list("slope" = score.sd.slopes, "plateau" = plateau.sd)
#   # # cheating - noise free data
#   # score.sd <- list("slope" = 0.01, "plateau" = 0.01)
#   # outlier.prob <- 0
#   
#   break.3 <- data.frame(
#     pattern = c(1, 0, NA),
#     bp.x = cumsum(c(0, 77, break12.mean)),  # psi, time duration of each phase (from reproduced study)
#     bp.y = c(0, 9.3, 9.3),                  # height of breakpoints
#     bp.x.sd = c(0, break01.sd, break12.sd), # noise levels
#     bp.y.sd = c(0, 1, 0)
#   )
#   true.traj <- trueTraj(times[["value"]], break.3[1:3])
#   
#   
#   ## REP LOOP on repetitions
#   for (jj in 1:n.sim){
#     
#     ## increment
#     # browser()
#     sim.nb = sim.nb + 1
#     
#     # # debug - go to a specific iteration
#     # if (sim.nb!=11) next
#     
#     ## 0. update 'progress bar' information
#     if (sim.nb %% 10 == 0 | sim.nb %in% c(2, nsimAll)){ print(paste(sim.nb, "/", nsimAll)) }
#     
#     ## 1. set seed
#     set.seed(allseeds[sim.nb])
#     
#     ## 2. simulate data and store the truth
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
#     # store individual level truth
#     true.parameters[
#       true.parameters$simID == sim.nb, c("patientID", "break.x1", "break.y1")
#     ] <- sim.data$sim.gen.model[,c("ID", "break.x1", "break.y1")]
#     
#     
#     ## 3. evaluate models
#     #  3.1. mixedBreak11
#     tryCatch(
#       {
#         # browser()
#         mod.11 <- suppressMessages(mixed.break1(
#           score ~ 0 + time + (0+time|ID), 
#           data = sim.data$sim.dataset,
#           psi0 = NA, # default psi initialization
#           it.max = 10L,
#           pattern = "11"
#         ))
#       },
#       error = function(e){
#         # print("Error during estimation of 'mixedBreak11' model.")
#         # browser()
#         estimates.11[estimates.11$simID == sim.nb, 'error'] <<- e$message
#       }
#     )
#     
#     #  3.2. mixedBreak10
#     tryCatch(
#       {
#         # browser()
#         mod.10 <- suppressMessages(mixed.break1(
#           score ~ 0 + time + (0+time|ID), 
#           data = sim.data$sim.dataset,
#           psi0 = NA, # default psi initialization
#           it.max = 10L,
#           pattern = "10"
#         ))
#       },
#       error = function(e){
#         # print("Error during estimation of 'mixedBreak 10' model.")
#         # browser()
#         estimates.10[estimates.10$simID == sim.nb, 'error'] <<- e$message
#       }
#     )
#     
#     
#     ## 4. store results / errors
#     # 4.1. mixedBreak11 results
#     # browser()
#     if(all(is.na(estimates.11[estimates.11$simID==sim.nb, 'error']))){
#       
#       res.11 <- normalize(mod.11)
#       # population level values
#       col.11 <- c(
#         paste0("break.", c("x1.avg", "x1.sd", "x1.random.sd", "CI95.low", "CI95.up")), 
#         paste0("slope2", c("", ".sd"))
#       )
#       estimates.11[estimates.11$simID==sim.nb, col.11] <- matrix(
#         unlist(res.11[col.11], use.names = F), 
#         nrow = n.patients, ncol = length(col.11), byrow = T
#       )
#       # individual breakpoint coordinates
#       estimates.11[estimates.11$simID==sim.nb, paste0("break.", c("x1", "y1"), ".ind")] <- 
#         unlist(res.11[paste0("break.", c("x1", "y1"), ".ind")], use.names = F)
#     }
#     
#     # 4.2. mixedBreak10 results
#     if(all(is.na(estimates.10[estimates.10$simID==sim.nb, 'error']))){
#       
#       # browser()
#       res.10 <- normalize(mod.10)
#       
#       # # only individual level values
#       # col.10 <- paste0("break.", c("x1.avg", "x1.sd", "CI95.low", "CI95.up", "y1.ind"))
#       # estimates.10[estimates.10$simID==sim.nb, col.10] <- res.10[,col.10]
#       
#       col.10 <- c(
#         paste0("break.", c("x1.avg", "x1.sd", "x1.random.sd", "CI95.low", "CI95.up"))
#       )
#       estimates.10[estimates.10$simID==sim.nb, col.10] <- matrix(
#         unlist(res.10[col.10], use.names = F), 
#         nrow = n.patients, ncol = length(col.10), byrow = T
#       )
#       # individual breakpoint coordinates
#       estimates.10[estimates.10$simID==sim.nb, paste0("break.", c("x1", "y1"), ".ind")] <- 
#         unlist(res.10[paste0("break.", c("x1", "y1"), ".ind")], use.names = F)
#     }
#     
#   }
#   # END REP LOOP
#   
# }
# ending.seed <- .Random.seed
# end.time <- Sys.time()
# rm(ii, jj, sim.nb, break.3, score.sd, true.traj, sim.data)
# rm(n.patients, score.sd.slopes, break01.sd, break12.sd, break12.mean)
# rm(res.10, res.11, mod.10, mod.lme, mod.11)
# 
# # END
# 
# print("Simulation study lasted for:")
# print(end.time-start.time)
# # Time difference of 10.48273 hours
# # print("Average time/iteration:")
# # print((end.time-start.time)/nsimAll)
# # print("Estimated simulation study duration (in hours):")
# # print((end.time-start.time)/nsimAll*n.scenario*594/60)
# # Estimated ~8.5h of simulation
# 
# 
# # # ## WRITING RESULTS IN FILE
# # library(writexl)
# # write_xlsx(estimates.11, "data/estimates_11.xlsx")
# # write_xlsx(estimates.10, "data/estimates_10.xlsx")
# # write_xlsx(list(true.parameters=true.parameters, scenario.data = scenario.data),
# #            "data/true_parameters_mixedBreak.xlsx")
# # write_xlsx(list(all.seeds=data.frame(allseeds), ending.seed=data.frame(ending.seed)),
# #            "data/random_seed_mixedBreak.xlsx")
# # # last 2 should be the same as before, right ?