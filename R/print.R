library(dplyr)

print.trajData <- function(sim.data, default = FALSE, dgm = T, n.sim = 10L,
                           cluster = 1:nrow(sim.data$sim.dataset)){
  if(default){
    print.data.frame(sim.data)
  }else{
    # Prints simulated data
    if(n.sim>0L){
      print(paste0("Simulated data head (", n.sim," lines) :"))
      print.data.frame(sim.data$sim.dataset |> head(n.sim))
    }
    
    # Prints Data Generation Model
    if(dgm){
      print("Data Generation Model:")
      
      # breakpoints coordinates 
      break.data <- sim.data$sim.gen.model %>%
        filter(ID %in% cluster) %>%
        select(-starts_with("beta.")) %>%
        pivot_longer(cols = starts_with("break.x"), names_to = "break.x", values_to = "x.coord") %>%
        pivot_longer(cols = starts_with("break.y"), names_to = "break.y", values_to = "y.coord") %>%
        filter(str_sub(break.x, -1L) == str_sub(break.y, -1L))
      
      break.data %>%
        group_by(ID) %>%
        mutate(bp.coord = paste0("(", round(x.coord, 1), ", ", round(x.coord, 1), ")")) %>%
        summarise(breakpoints = paste(bp.coord, collapse = ", "))
      
    }
  }
}

