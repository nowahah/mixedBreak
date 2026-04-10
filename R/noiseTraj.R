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


noiseTraj <- function(trueTraj, n.obs, score.sd, times.sd,
                      breakpoints.sd = NULL, break.x.dist = 0,
                      outlier.prob = 0, na.prob = 0, n.trail = 0L){
  
}