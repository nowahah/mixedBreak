### simData10.R ---
## * simData10 (documentation)
##' @title Simulate "10"-alike trajectories with known parameters
##' @description Simulates trajectories of SDI in patients after psilocybin intake.
##' Any "10"-alike pattern can be specified with breakpoints coordinates and noise level.
##' The SDI scale is an integer subjective scale ranging from 0-10.
##'
##' @param n.obs [integer >= 1] number of observations to simulate in the dataset.
##' @param score.sd [numeric>0] noise level (sd) on the measurement.
##' TODO - later could be a vector to specify noise level in each phase of the experiment
##' (typically lower during plateau phase)
##' @param pattern 
##' @param times [named list or integer vector] times and noise (if any) of simulated measurements,
##' for regular strategies.
##' Can be a vector with values of time measurements to specify other irregular strategies.
##'
##' @param breakpoints [data.frame] data.frame with \code{length(pattern)+2} rows
##' and 5 columns: \itemize{
##' \item "pattern" shape of the pattern. Should only contain 1s or 0s, 
##' never two 0 in a row. The last value is ignored and set to \code{NA}.
##' The pattern of the trajectory has as many segments as \code{length(pattern)},
##' and 1 breakpoint more than the length of pattern (beginning and ending are considered breakpoints). 
##' Specifying 1 means a slope coefficients is expected to be non-zero, and 0 means a plateau (slope=0). 
##' \item "bp.x" every breakpoints x coordinates on the time scale (sorted,
##' should be between 0 and the end of the trajectory). 
##' \item "bp.y" every breakpoints y coordinates on the SDI scale (between 0 and 10). 
##' If \code{pattern[i]==0}, then there should be \code{bp.y[i-1] == bp.y[i]}.
##' \item "bp.x.sd" noise (sd) parameter for breakpoints' x coordinates.
##' \item "bp.y.sd" noise (sd) parameter for breakpoints' y coordinates.
##' If \code{pattern[i]==0}, there should be \code{bp.y.sd[i]==0}
##' }
##' In case of conflict between pattern and breakpoints specifications, the pattern prevails
##' and the previous breakpoint specification will override the next one. 
##' A warning message will be displayed.
##'
##' @param break.min.dist [list of numeric] Minimum between-breakpoints distance on
##' x & y coordinates. A list with 2 items: \itemize{
##' \item "x" for x lower bound on distance
##' \item "y" *IGNORED* for y lower bound on distance
##' }
##' @param na.prob [0<numeric<1] proportion of missing at random values in the dataset.
##' @param outlier.prob [0<numeric<1] proportion of outlier in the dataset.
##' Outliers are simulated as a random drop of 1-3 points in SDI.
##' @param n.trail [integer >= 1] number of trailing observations at the end of the trajectory
##'
##'
##' @return A dataset with n.obs trajectories of SDI, simulated according to specified parameters
##' It contains ID, time of measurements, Data Generation Model, simulated measurements
##' and other characteristics such as outlier or missing value flag


simData10 <- function(n.obs, score.sd, times = list("value" = 20, "sd" = 0),
                      breakpoints = NULL, break.min.dist = list(x = 0, y = NA),
                      outlier.prob = 0, na.prob = 0, n.trail = 0L) {
  require(dplyr)
  
  ## check entries of the function
  # TODO
  pattern = breakpoints[["pattern"]]
  
  ## normalize parametrization of trajectories
  # TODO
  
  ## simulate breakpoints coordinates
  n.break <- nrow(breakpoints) ## nb of breakpoints
  break.x.value = matrix(rnorm(n.obs*n.break, breakpoints[["bp.x"]], breakpoints[["bp.x.sd"]]), 
                         nrow=n.break,dimnames = list(paste0("break.x", 1:n.break), 1:n.obs))
  
  # check breakpoint validity - sorted in time for every individual
  # instead of 0 here we could add a minimal distance between breakpoints
  while(any(diff(break.x.value)<break.min.dist$x)){
    # which ind have unsorted break points
    ind.retry <- which(diff(break.x.value)<break.min.dist$x, arr.ind = T)[, 2]
    n.retry <- length(ind.retry)
    # for those ind, re-simulate a vector of breakpoints
    break.x.value[, ind.retry] = matrix(rnorm(n.retry*n.break, breakpoints[["bp.x"]], 
                                              breakpoints[["bp.x.sd"]]), 
                                        nrow = n.break,
                                        dimnames = list(paste0("break.x", 1:n.break), 
                                                        ind.retry))
  }
  
  # simulate height / y-coordinate of breakpoints
  break.y.value = matrix(rnorm(n.obs*n.break, breakpoints[["bp.y"]], breakpoints[["bp.y.sd"]]), 
                         nrow=n.break,dimnames = list(paste0("break.y", 1:n.break), 1:n.obs))
  # force height of breakpoints after plateau
  if(any(pattern==0)){
    plateau.i <- which(pattern==0) # which segments are plateaus
    
    ## warning in case of over-parametrization
    # case 1: sd is specified after plateau (but height is forced)
    over.param <- breakpoints[["bp.y.sd"]][plateau.i+1]>0
    # case 2: pre-plateau and post-plateau means are different
    over.param <- over.param | (breakpoints[["bp.y"]][plateau.i+1]!=breakpoints[["bp.y"]][plateau.i])
    if(any(over.param)){
      warning(paste0("Over-parametrization: y-distribution (mean, sd) for breakpoint(s) {", 
                     paste((plateau.i+1)[which(over.param)], collapse = ", "),
                     "} is/are ignored because pattern '", 
                     paste(pattern[1:(n.break-1)], collapse = ""),
                     "' indicates a plateau for the considered segment(s)."))
    }
    
    # plateau overrides over-specified distribution
    break.y.value[plateau.i+1,] <- break.y.value[plateau.i,]
  }
  break.y.value <- pmax(0, pmin(10, break.y.value))
  break.y.value <- matrix(break.y.value, nrow = n.break)
  
  ## Compute slopes coefficients
  slopes <- diff(break.y.value) / diff(break.x.value)
  slopes <- slopes[plateau.i-1,] # remove plateau (slope=0)
  slopes <- as.data.frame(t(slopes))
  names(slopes) <- paste0("beta.", plateau.i-1)
  
  ## Compute ending time of trajectories
  ending.times <- break.x.value[n.break,]
  max.time <- max(ending.times)
  
  
  ## Derive times of measurements
  if (length(times[["value"]]) == 1) { # regular measurements
    time <- rep(seq(0, max.time + times[["value"]], by = times[["value"]]), n.obs)
    time.step <- times[["value"]]
  } else { # irregular measurements
    time <- rep(times[["value"]], n.obs)
    time.step <- min(diff(times[["value"]]))
  }
  n.time <- length(time) / n.obs
  # Adding noise to the time measurements
  if (times[["sd"]] > 0) {
    time[time != 0] <- round(pmax(time[time != 0] + rnorm(length(time[time != 0]), sd = times[["sd"]]), 0))
  }
  
  
  ## Simulated dataset
  sim.dataset <- data.frame(
    ID = as.factor(rep(1:n.obs, each = n.time)),
    time = time
  )
  ## Data Generation Model
  sim.gen.model <- data.frame(
    ID = as.factor(1:n.obs)
  )
  # add breakpoints coordinates
  sim.gen.model <- cbind(sim.gen.model, t(break.x.value), t(break.y.value))
  names(sim.gen.model)[2:11] <- paste0(c("start", paste0("break", 1:(n.break-2)), "end"), 
                                       rep(c(".x", ".y"), each=n.break))
  
  ## Simulate the "true" and "measured" trajectory (noise added)
  # TODO
  # sim.dataset <- sim.dataset %>%
  #   mutate(
  #     # compute "true perfect" trajectory (101 pattern)
  #     # TODO rewrite with a case_when
  #     truth.101 = if_else(time < break.1, b.onset * (time / 20),
  #                         if_else(time < break.2, plateau, plateau + b.return * (time - break.2) / 20)
  #     ),
  #     noised.101 = pmin(10, pmax(0, truth.101 + rnorm(n(), sd = score.sd))),
  #     # truncated end of trajectory
  #     score = if_else(time < end.time + (1 + n.trail) * time.step, pmin(10, pmax(0, round(noised.101))), NA),
  #     # trailing at the end: end.time + (1+n.trail)*time.step (to derive)
  #     truth.101 = pmin(10, pmax(0, truth.101))
  #   )
  
  
  ## Adding outliers - patient went to the bathroom
  # Can happen at most once per patient
  # sim.dataset <- sim.dataset |>
  #   group_by(ID) |>
  #   mutate(
  #     # outliers location and draft
  #     outliers = (row_number() == sample.int(n(), 1)) * rbinom(n(), 1, outlier.prob) == 1, 
  #     # outliers is a drop in truth of 1-3 points (still bounded 0<=score<=10)
  #     score = pmin(10, pmax(0, score - sample(1:3, 1, prob = c(1, 4, 4)) * outliers)),
  #     outliers = if_else(is.na(score), NA, outliers)
  #   )
  
  
  ## Adding missing values (Missing Completely At Random)
  # sim.dataset <- sim.dataset %>%
  #   mutate(
  #     na.value = if_else(is.na(score), NA, rbinom(n(), 1, na.prob) == 1),
  #     score = if_else(na.value, NA, score)
  #   )
  
  
  ## export
  class(sim.dataset) <- append("trajData", class(sim.dataset))
  return(list(sim.dataset = sim.dataset, sim.gen.model = sim.gen.model))
}

# breakpoints |> mutate(duration = c(diff(bp.x), NA), slope = c(diff(bp.y)/diff(bp.x), NA))
