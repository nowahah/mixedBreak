### noiseTraj.R ---
## * noiseTraj (documentation)
##' @title Simulate "10"-alike trajectories from known parameters and noise level
##' @description Simulates trajectories of SDI in patients after psilocybin intake from
##' a dataset that contains the true values of parameters. 
##' Noise level can be specified for every parameter specified in the reference dataset.
##'
##' @param n.obs [integer >= 1] number of observations to simulate in the dataset.
##' @param score.sd [numeric>0] noise level (sd) on the measurement/score.
##' TODO - later could be a vector to specify noise level in each phase of the experiment
##' (typically lower during plateau phase)
##' @param times.sd [numeric>0]
##'
##' @param breakpoints.sd [data.frame] data.frame with \code{length(pattern)+1} rows
##' and 2 columns: \itemize{
##' \item "bp.x.sd" noise (sd) parameter for breakpoints' x coordinates.
##' \item "bp.y.sd" noise (sd) parameter for breakpoints' y coordinates.
##' If \code{pattern[i]==0}, there should be \code{bp.y.sd[i]==0}.
##' }
##' In case of conflict between pattern and breakpoints specifications, the pattern prevails
##' and the previous breakpoint specification will override the next one.
##' A warning message will be displayed.
##'
##' @param break.x.dist [numeric>0] Minimum between-breakpoints distance on
##' x coordinates. Default is 0.
##' @param na.prob [0<numeric<1] proportion of missing at random values in the dataset.
##' @param outlier.prob [0<numeric<1] proportion of outlier in the dataset.
##' Outliers are simulated as a random drop of 1-3 points in SDI, with higher probability
##' on 
##' @param n.trail [integer] number of trailing observations at the end of the trajectory.
##' Default is 0. Has a meaning when positive. 
##' If negative, it will remove the last observations from the trajectory.
##' Ignored if pattern ends with a 0 / plateau.
##'
##'
##' @return A list of 2 datasets: \itemize{
##' \item "sim.dataset" contains ID, time, true, noised measurements and score,
##' as well as outlier / missing / trailing flags
##' \item "sim.gen.model" contains ID and patients-specific breakpoints coordinates
##' }


noiseTraj <- function(true.traj, n.obs, score.sd, times.sd,
                      breakpoints.sd = NULL, break.x.dist = 0,
                      outlier.prob = 0, na.prob = 0, n.trail = 0L){
  require(dplyr)
  require(tidyr)
  
  ## check entries of the function
  # TODO
  # check type of entries, length
  
  # if(length(score.sd)==1){
  #   score.sd <- setNames(rep(score.sd, 2), c("slope","plateau"))
  # }else if(length(score.sd)==2){
  #   if(is.null(names(score.sd))){
  #     names(score.sd) <- c("slope","plateau")
  #   }else if(any(sort(names(score.sd))!=c("slope","plateau"))){
  #     stop("Argument \'score.sd\' should have name \"slope\" and \"plateau\". \n")
  #   }
  # }else if(length(score.sd)>3){
  #   stop("Argument \'score.sd\' too long (length>2). \n")
  # }
  
  # extract information contained in 'trueTraj' object
  breakpoints <- cbind(true.traj$breakpoints, breakpoints.sd)
  times <- list(value = true.traj$times, sd = times.sd)
  
  pattern <- breakpoints[["pattern"]]
  if(pattern[length(pattern)-1]==0 & n.trail > 0L){ 
    warning(paste0("Trailling is ignored because pattern '",
                   str_sub(paste(pattern, collapse = ""), end=-3),
                   "' ends with a plateau."))
    n.trail <- -1L 
  }
  
  ## normalize parametrization of trajectories
  # TODO
  
  ## simulate between-breakpoints distances to compute coordinates
  n.break <- nrow(breakpoints) ## nb of breakpoints
  mu <- diff(breakpoints[["bp.x"]])
  std <- breakpoints[["bp.x.sd"]][2:n.break]
  
  ## If we would like to simulate correlated breakpoints  
  ## -> mvtnorm::rmvnorm
  break.x.dist <- matrix(rlnorm(n.obs * (n.break-1), log(mu^2/sqrt(mu^2+std^2)), sqrt(log(1+(std/mu)^2))),
                         nrow = n.break-1, dimnames = list(paste0("break.x", 1:(n.break-1)), 1:n.obs)
  )
  
  ## check distance validity - bounded by break.min.dist$x
  while (any(break.x.dist < break.min.dist$x)) {
    # which ind have unsorted break points
    ind.retry <- which(break.x.dist < break.min.dist$x, arr.ind = T)[, 2]
    n.retry <- length(ind.retry)
    # for those ind, re-simulate a vector of breakpoints
    break.x.dist[, ind.retry] <- matrix(
      rlnorm(
        n.retry * (n.break-1), log(mu^2/sqrt(mu^2+std^2)), sqrt(log(1+(std/mu)^2))),
      nrow = n.break-1, dimnames = list(paste0("break.x", 1:(n.break-1)), ind.retry)
    )
  }
  
  break.x.value <- rbind(
    rep(0, n.obs),
    apply(break.x.dist, 2, cumsum)
  )
  
  
  ## simulate height / y-coordinate of breakpoints 
  # w\ truncated normal distribution
  slopes.i  <- which(pattern == 1) # which segments are slopes
  n2sim <- length(slopes.i)
  break.y.tmp <- matrix(
    EnvStats::rnormTrunc(n.obs * n2sim, min=0, max=10,
                         mean = breakpoints[["bp.y"]][slopes.i+1], 
                         sd = breakpoints[["bp.y.sd"]][slopes.i+1]),
    nrow = n2sim, dimnames = list(paste0("break.y", slopes.i), 1:n.obs)
  )
  
  
  ## force height of breakpoints after plateau
  plateau.i <- which(pattern == 0) # which segments are plateaus
  if (any(pattern == 0, na.rm = T)) {
    ## warning in case of over-parametrization
    # case 1: sd is specified after plateau (but height is forced)
    over.param <- breakpoints[["bp.y.sd"]][plateau.i + 1] > 0
    # case 2: pre-plateau and post-plateau means are different
    over.param <- over.param | (breakpoints[["bp.y"]][plateau.i + 1] != breakpoints[["bp.y"]][plateau.i])
    if (any(over.param)) {
      warning(paste0(
        "Over-parametrization: y-distribution (mean, sd) for breakpoint(s) {",
        paste((plateau.i + 1)[which(over.param)], collapse = ", "),
        "} is/are ignored because pattern '",
        paste(pattern[1:(n.break - 1)], collapse = ""),
        "' indicates a plateau for the considered segment(s)."
      ))
    }
  }
  ## TODO - plateau height is forced for its last breakpoint
  ## careful to take into account all possible patterns
  # DOESN't WORK for pattern 1010 - repair
  break.y.value <- matrix(0, nrow = n.break, ncol = n.obs,
                          dimnames = list(paste0("break.y", 1:n.break-1), 1:n.obs))
  break.y.value[-c(1, 1+plateau.i),] <- break.y.tmp
  break.y.value[plateau.i + 1, ] <- break.y.value[plateau.i, ]
  
  ## Compute slopes coefficients
  slopes <- diff(break.y.value) / diff(break.x.value)
  slopes.i <- which(pattern == 1)
  slopes <- slopes[slopes.i, ] # remove plateau(s) (slope=0)
  slopes <- as.data.frame(t(slopes))
  names(slopes) <- paste0("beta.", slopes.i)
  
  ## Compute ending time of trajectories
  ending.times <- break.x.value[n.break, ]
  max.time <- max(ending.times)
  
  
  ## Derive times of measurements
  if (length(times[["value"]]) == 1) { # regular measurements
    time.step <- times[["value"]]
    
    # nb of measurements per ID (assuming starting at time=0)
    n.time <- ceiling(ending.times/time.step + n.trail + 1)
    time <- sequence(n.time, from = 0)*time.step
    
  } else { # irregular measurements
    time.step <- min(diff(times[["value"]]))
    time <- rep(c(times[["value"]], last(times[["value"]]) + seq(1, n.trail)*time.step), n.obs)
    
    # nb of measurements per ID (assuming starting at time=0)
    n.time <- length(times[["value"]])+n.trail
    
  }
  # Adding noise to the time measurements
  if (times[["sd"]] > 0) {
    time[time != 0] <- round(pmax(time[time != 0] + rnorm(length(time[time != 0]), sd = times[["sd"]]), 0))
  }
  
  
  ## Simulated dataset
  sim.dataset <- data.frame(
    ID = as.factor(rep(1:n.obs, times = n.time)),
    time = time
  )
  ## Data Generation Model
  sim.gen.model <- data.frame(
    ID = as.factor(1:n.obs)
  )
  # add breakpoints coordinates
  sim.gen.model <- cbind(sim.gen.model, t(break.x.value), t(break.y.value), slopes)
  names(sim.gen.model)[1+1:(2*n.break)] <- paste0(
    rep("break", n.break),
    paste0(rep(c(".x", ".y"), each = n.break), (1:n.break - 1))
  )
  
  # compute segments number
  sim.dataset <- sim.dataset %>%
    left_join(sim.gen.model %>%
                select(c(1, starts_with("break."))), by = "ID") %>%
    rowwise() %>% # row-wise cut on time against break x-values
    mutate(
      segment = cut(time,
                    breaks = c_across(starts_with("break.x")),
                    labels = 1:(n.break - 1), include.lowest = T
      ),
      pattern = pattern[segment]
    )
  
  # compute true and noised trajectories
  sim.dataset <- sim.dataset %>%
    pivot_longer(cols = starts_with("break.")) %>%
    group_by(ID, segment) %>%
    filter(name %in% paste0(rep(paste0("break.", c("x", "y")), each=2), 
                            as.numeric(segment)+c(-1,0)) |
             (is.na(segment) & (name %in% paste0(rep(paste0("break.", c("x", "y")), each=2), 
                                                 n.break - c(2, 1))))) %>%
    mutate(
      x1 = first(value), x2 = nth(value, 2), 
      y1 = nth(value, 3), y2 = nth(value, 4),
      true.traj = case_when(
        is.na(pattern) ~ y2,
        pattern == 0 ~ y1,
        pattern == 1 ~ approx(c(x1, x2), c(y1, y2), xout = time, ties = "ordered")$y
      )
    ) %>%
    ungroup() %>%
    select(-c(name, value, x1, x2, y1, y2)) %>%
    distinct()
  
  
  # Adding noised and measured trajectories (integer-round)
  sim.dataset <- sim.dataset %>%
    mutate(
      noised.traj = true.traj + rnorm(n(), sd = score.sd),
      # rounded score, including trailing observations
      score = if_else(time < ending.times[ID] + (1 + n.trail) * time.step, 
                      pmin(10, pmax(0, round(noised.traj))), NA)
      # trailing is checked bcs due to noise on time of measurements, the nb of trailing may change
      # TODO - should trailling be noised as well ?
    )
  
  
  ## Adding outliers - patient went to the bathroom
  # Can happen at most once per patient
  sim.dataset <- sim.dataset %>%
    group_by(ID) %>%
    mutate(
      # outliers location and draft
      outliers = (row_number() == sample.int(n(), 1)) * rbinom(n(), 1, outlier.prob) == 1,
      # outliers is a drop in truth of 1-3 points (still bounded 0<=score<=10)
      score = pmin(10, pmax(0, score - sample(1:3, 1, prob = c(1, 4, 4)) * outliers)),
      outliers = if_else(is.na(score), NA, outliers)
    ) %>%
    ungroup()
  
  
  ## Adding missing values (Missing Completely At Random)
  sim.dataset <- sim.dataset %>%
    mutate(
      na.flag = if_else(is.na(score), NA, rbinom(n(), 1, na.prob) == 1),
      score = if_else(na.flag, NA, score),
      type = case_when(
        outliers ~ "outlier",
        is.na(pattern) ~ "trailing",
        .default = "signal"
      )
    )
  
  
  ## export
  sim.data <- list(sim.dataset = sim.dataset, sim.gen.model = sim.gen.model)
  class(sim.data) <- append("trajData", class(sim.dataset))
  return(sim.data)
}

