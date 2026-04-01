### simData.R ---
## * simData (documentation)
##' @title Simulate SDI trajectories with known truth and noise
##' @description Simulates trajectories of SDI in patients after psilocybin intake.
##' The simulated trajectories have 3 separate linear phases: onset, peak/plateau,
##' return to normal consciousness. There are many different parametric representations of such
##' trajectories. In case it is over specified, *the first parameters will prevail*.
##' The SDI scale is an integer subjective scale ranging from 0-10.
##'
##' @param n.obs [integer >= 1] number of observations to simulate in the dataset.
##' @param score.sd [numeric>0] noise level (sd) on the measurement.
##' TODO - later could be a vector to specify noise level in each phase of the experiment
##' (typically lower during the plateau phase)
##' @param times [named list or integer vector] times and noise (if any) of simulated measurements,
##' for regular strategies.
##' Can be a vector with values of time measurements to specify other irregular strategies.
##'
##' @param plateau [named list] named list with \itemize{
##' \item "value" for the height / altitude of the plateau (e.g. its 'y' coordinate).
##' \item "value.sd" for the noise level of height of plateau.
##' It is *highly recommended* to specify plateau's height given its constraint (\leq 10),
##' to avoid any unintended re-computation of \code{b.onset} or \code{break.1}. 
##' They would then not match the pre-specified distribution for those parameters.
##' \item "duration" for the duration of the plateau on the same time scale as \code{times}.
##' If \code{plateau$duration==0}, then it is assumed that there is no plateau, and the 
##' second breakpoint value is ignored if provided.
##' \item "duration.sd" for the noise level of duration of plateau (on the same time scale 
##' as \code{times}).
##' \item "duration.min" Minimum duration of simulated plateau, if any (always >= 0)
##' \item "duration.max" Maximum duration of simulated plateau, if any
##' }
##' @param b.onset [named list] named list with \itemize{
##' \item "value" for coefficient value,
##' \item "sd" for the noise level,
##' \item "min" for minimum value simulated.
##' }
##' @param b.return [named list] named list with \itemize{
##' \item "value" for coefficient value,
##' \item "sd" for the noise level,
##' \item "max" for maximum value simulated.
##' }
##' @param breakpoints [data.frame] data.frame with \itemize{
##' \item "bp" first column containing the two breakpoints (sorted,
##' should be between 0 and the end of the trajectory)
##' \item "bp.sd" second column with noise parameter for each breakpoint.
##' }
##' @param ending.times [named list] named list with \itemize{
##' \item "value" for average ending time of effects on the patients
##' \item "sd" for the noise level of parameter
##' \item "max" for maximum value simulated.
##' }
##'
##' @param na.prob [0<numeric<1] proportion of missing at random values in the dataset.
##' @param outlier.prob [0<numeric<1] proportion of outlier in the dataset.
##' Outliers are simulated as a random drop of 1-3 points in SDI.
##' @param n.trail [integer >= 1] number of trailing observations at the end of the trajectory
##'
##'
##' @return A dataset with n.obs trajectories of SDI, simulated according to specified parameters
##' It contains ID, time of measurements, Data Generation Model, simulated measurements
##' and other characteristics such as outlier or missing value flag


simData <- function(n.obs, score.sd, times = list("value" = 20, "sd" = 0),
                    plateau = list("value" = 9, "value.sd" = 3,
                                   "duration" = 82.4, "duration.sd" = 33.1), 
                    b.onset = NULL, b.return = NULL, 
                    breakpoints = NULL, ending.times = NULL,
                    outlier.prob = 0, na.prob = 0, n.trail = 0L) {
  require(dplyr)

  ## check entries of the function
  # TODO

  ## normalize parametrization of trajectories
  # TODO
  stopifnot("value" %in% names(plateau), "Plateau duration is mandatory to parametrize trajectories")
  stopifnot(!is.null(b.onset) | NA, "b.onset or break.1 is mandatory to parametrize the curve")


  ## simulate breakpoints value
  n.break <- nrow(breakpoints) ## nb of breakpoints
  bp.list <- lapply(1:n.break, function(i) {
    rnorm(n.obs, mean = breakpoints$bp[i], sd = breakpoints$bp.sd[i])
  })
  break.value <- as.data.frame(setNames(bp.list, paste0("break.", 1:n.break))) # every bp value once


  ## Simulate slopes coefficients
  coef.onset <- pmax(b.onset[["value"]] + rnorm(n.obs, sd = b.onset[["sd"]]), b.onset[["min"]])
  coef.return <- pmin(b.return[["value"]] + rnorm(n.obs, sd = b.return[["sd"]]), b.return[["max"]])


  ## check breakpoint validity
  # If plateau is computed above 10 - the upper bound of the scale here -
  # we recompute the first breakpoint so the plateau's value is exactly 10.
  peaks <- coef.onset / 20 * break.value[[1]]
  break.value[[1]][peaks > 10] <- 10 / (coef.onset[peaks > 10] / 20)


  ## Compute ending time of trajectories
  # !!! Only works for patterns "101" and "11"
  ending.times <- break.value[[ncol(break.value)]] - coef.onset * break.value[[1]] / coef.return
  end.time <- max(ending.times)


  ## Derive times of measurements
  if (length(times[["value"]]) == 1) { # regular measurements
    time <- rep(seq(0, end.time + times[["value"]], by = times[["value"]]), n.obs)
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
    time = time,
    b.onset = rep(coef.onset, each = n.time),
    b.return = rep(coef.return, each = n.time),
    return.0 = rep(ending.times, each = n.time)
  )
  # add and rep breakpoints values and ending times
  sim.dataset <- cbind(sim.dataset, apply(break.value, 2, rep, each = n.time))


  ## Simulate the "true" and "measured" trajectory (noise added)
  # !!! Only for 101 model
  sim.dataset <- sim.dataset %>%
    mutate(
      plateau = pmin(10, b.onset * (break.1 / 20)),
      # compute "true perfect" trajectory (101 pattern)
      # TODO rewrite with a case_when
      truth.101 = if_else(time < break.1, b.onset * (time / 20),
        if_else(time < break.2, plateau, plateau + b.return * (time - break.2) / 20)
      ),
      noised.101 = pmin(10, pmax(0, truth.101 + rnorm(n(), sd = score.sd))),
      # truncated end of trajectory
      score = if_else(time < return.0 + (1 + n.trail) * time.step, pmin(10, pmax(0, round(noised.101))), NA),
      # trailing at the end: return.0 + (1+n.trail)*time.step (to derive)
      truth.101 = pmin(10, pmax(0, truth.101))
    )


  ## Adding outliers - patient went to the bathroom
  # Can happen at most once per patient
  sim.dataset <- sim.dataset |>
    group_by(ID) |>
    mutate(
      outliers = (row_number() == sample.int(n(), 1)) * rbinom(n(), 1, outlier.prob) == 1, # outliers location and draft
      # outliers is a drop in truth of 1-3 points (still bounded 0<=score<=10)
      score = pmin(10, pmax(0, score - sample(1:3, 1, prob = c(1, 4, 4)) * outliers)),
      outliers = if_else(is.na(score), NA, outliers)
    )


  ## Adding missing values (Missing Completely At Random)
  sim.dataset <- sim.dataset %>%
    mutate(
      na.value = if_else(is.na(score), NA, rbinom(n(), 1, na.prob) == 1),
      score = if_else(na.value, NA, score)
    )


  ## export
  class(sim.dataset) <- append("trajData", class(sim.dataset))
  return(sim.dataset)
}
