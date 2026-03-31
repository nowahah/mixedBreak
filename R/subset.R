subset.trajData <- function(x, observed = TRUE){
  if(!observed){
    subset.data.frame(x)
  }
  
  ## Default method to subset the only informations that would have been measured 
  ## if we really did the experiment
  x |> select(-c(3:8))
}