### trueTraj.R ---
## * trueTraj.R (documentation)
##' @title Simulate "10"-alike trajectories with known true parameters
##' @description Specify the simulation model for SDI trajectories in patients after psilocybin intake.
##' Any "10"-alike pattern can be specified with breakpoints coordinates and noise level.
##' The SDI scale is an integer subjective scale ranging from 0-10.
##'
##' @param times [numeric] timestep of simulated measurements, for regular strategies.
##' Can also be a vector with values of time measurements to specify other irregular strategies.
##' @param breakpoints [data.frame] data.frame with \code{length(pattern)+1} rows
##' and 3 columns: \itemize{
##' \item "pattern" shape of the pattern. Should only contain 1s or 0s,
##' never two 0 in a row. The last value is ignored and set to \code{NA}.
##' The pattern of the trajectory has as many segments as \code{length(pattern)},
##' and 1 breakpoint more than the length of pattern (beginning and ending are considered breakpoints).
##' Specifying 1 means a slope coefficients is expected to be non-zero, and 0 means a plateau (slope=0).
##' \item "bp.x" every breakpoints x coordinates on the time scale (sorted,
##' should be between 0 and the end of the trajectory).
##' \item "bp.y" every breakpoints y coordinates on the SDI scale (between 0 and 10).
##' }
##' A warning message will displayed in case of conflict between pattern and average breakpoint value.
##' Pattern prevails on specified distributions.
##'
##'
##' @return A dataset with n.obs trajectories of SDI, simulated according to specified parameters
##' It contains ID, time of measurements, Data Generation Model, simulated measurements
##' and other characteristics such as outlier or missing value flag


trueTraj <- function(times = 20L, breakpoints = NULL){

    traj.truth <- list(times = times, breakpoints = breakpoints)
  
  class(traj.truth) <- append("trajTruth", class(traj.truth))
  return(traj.truth)
}

