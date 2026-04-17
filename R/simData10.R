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
##' @param times [named list or integer vector] times and noise (if any) of simulated measurements,
##' for regular strategies.
##' Can be a vector with values of time measurements to specify other irregular strategies.
##'
##' @param breakpoints [data.frame] data.frame with \code{length(pattern)+1} rows
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
##' \item "x" [numeric>0] for x lower bound on distance. Default is 0.
##' \item "y" *IGNORED* [numeric>0] for y lower bound on distance. Default is NA.
##' }
##' @param na.prob [0<numeric<1] proportion of missing at random values in the dataset.
##' @param outlier.prob [0<numeric<1] proportion of outlier in the dataset.
##' Outliers are simulated as a random drop of 1-3 points in SDI.
##' @param n.trail [integer] number of trailing observations at the end of the trajectory.
##' Default is 0. Has a meaning when positive. 
##' If negative, it will remove the last observations from the trajectory.
##' Ignored if pattern ends with a 0 / plateau
##'
##'
##' @return A dataset with n.obs trajectories of SDI, simulated according to specified parameters
##' It contains ID, time of measurements, Data Generation Model, simulated measurements
##' and other characteristics such as outlier or missing value flag

##' @export
simData10 <- function(n.obs, score.sd, times = list("value" = 20, "sd" = 0),
                      breakpoints = NULL, break.min.dist = list(x = 0, y = NA),
                      outlier.prob = 0, na.prob = 0, n.trail = 0L) {
  require(dplyr)
  require(tidyr)

  ## check entries of the function
  # TODO
  pattern <- breakpoints[["pattern"]]
  if(pattern[length(pattern)-1]==0 & n.trail > 0L){ 
    warning(paste0("Trailling is ignored because pattern '",
                   str_sub(paste(pattern, collapse = ""), end=-3),
                   "' ends with a plateau."))
    n.trail <- -1L 
  }

  ## normalize parametrization of trajectories
  # TODO

  ## simulate breakpoints coordinates
  n.break <- nrow(breakpoints) ## nb of breakpoints
  break.x.value <- matrix(rnorm(n.obs * n.break, breakpoints[["bp.x"]], breakpoints[["bp.x.sd"]]),
    nrow = n.break, dimnames = list(paste0("break.x", 1:n.break), 1:n.obs)
  )

  # check breakpoint validity - sorted in time for every individual
  # instead of 0 here we could add a minimal distance between breakpoints
  while (any(diff(break.x.value) < break.min.dist$x)) {
    # which ind have unsorted break points
    ind.retry <- which(diff(break.x.value) < break.min.dist$x, arr.ind = T)[, 2]
    n.retry <- length(ind.retry)
    # for those ind, re-simulate a vector of breakpoints
    break.x.value[, ind.retry] <- matrix(
      rnorm(
        n.retry * n.break, breakpoints[["bp.x"]],
        breakpoints[["bp.x.sd"]]
      ),
      nrow = n.break,
      dimnames = list(
        paste0("break.x", 1:n.break),
        ind.retry
      )
    )
  }

  # simulate height / y-coordinate of breakpoints
  break.y.value <- matrix(rnorm(n.obs * n.break, breakpoints[["bp.y"]], breakpoints[["bp.y.sd"]]),
    nrow = n.break, dimnames = list(paste0("break.y", 1:n.break), 1:n.obs)
  )
  # force height of breakpoints after plateau
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

    # plateau overrides over-specified distribution
    break.y.value[plateau.i + 1, ] <- break.y.value[plateau.i, ]
  }
  break.y.value <- pmax(0, pmin(10, break.y.value))
  break.y.value <- matrix(break.y.value, nrow = n.break)

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
  
  
  ## population level truth - breakpoints
  pop.lvl.truth <- list(
    times = times[["value"]],
    breakpoints = breakpoints[,c("pattern", "bp.x", "bp.y")]
  )
  class(pop.lvl.truth) <- append("trajTruth", class(pop.lvl.truth))
  
  
  ## export
  sim.data <- list(
    sim.dataset = sim.dataset, 
    sim.gen.model = sim.gen.model,
    pop.lvl.truth = pop.lvl.truth
  )
  class(sim.data) <- append("trajData", class(sim.dataset))
  return(sim.data)
}

# breakpoints |> mutate(duration = c(diff(bp.x), NA), slope = c(diff(bp.y)/diff(bp.x), NA))
