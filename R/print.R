library(dplyr)

print.trajData <- function(x, default = FALSE, n.line = 0L, dgm = T){
  if(default){
    print.data.frame(x)
  }else{
    # Prints simulated data
    if(n.line>0L){
      print(paste0("Simulated data head (", n.line,") :"))
      print.data.frame(x |> select(-c(3:8)) |> head(n.line))
    }
    
    # Prints Data Generation Model
    if(dgm){
      print("Data Generation Model:")
      print.data.frame(x |> select(c(1, 3:8)) |> distinct())
    }
  }
}

