subset.trajData <- function(x, observed = TRUE, default = FALSE, subset = NULL){
  if(default){
    subset.data.frame(x, subset = subset)
  }
  
  ## Default method to subset the only information that would have been measured 
  ## if we really did the experiment
  if (observed){
    x$sim.dataset |> select(c(ID, time, score))
  }
  
  # TODO - needs to handle more cases
}