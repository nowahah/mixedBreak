library(dplyr)

print.trajData <- function(sim.data, default = FALSE, dgm = T, n.lines = 10L,
                           cluster = levels(sim.data$sim.dataset$ID)){
  if(default){
    print.data.frame(sim.data)
  }else{
    # Prints simulated data
    if(n.lines>0L){
      print(cat(paste0("\nSimulated data head(", n.lines,"):")))
      print.data.frame(
        sim.data$sim.dataset %>% 
          filter(ID %in% cluster) %>% 
          head(n.lines)
        )
    }
    
    # Prints Data Generation Model
    if(dgm){
      pattern <- sim.data$sim.dataset
      print(cat("\nData Generation Model - pattern "))
      
      # breakpoints coordinates 
      break.data <- sim.data$sim.gen.model %>%
        filter(ID %in% cluster) %>%
        select(-starts_with("beta.")) %>%
        pivot_longer(cols = starts_with("break.x"), names_to = "break.x", values_to = "x.coord") %>%
        pivot_longer(cols = starts_with("break.y"), names_to = "break.y", values_to = "y.coord") %>%
        filter(str_sub(break.x, -1L) == str_sub(break.y, -1L))
      
      # printed result
      break.data %>%
        filter(ID %in% cluster) %>%
        group_by(ID) %>%
        mutate(
          bp.coord = paste0("(", round(x.coord, 1), ", ", round(y.coord, 1), ")"),
          bp.nb = paste0("break.", str_extract(break.x, "\\d"))
        ) %>%
        pivot_wider(id_cols=ID, names_from = "bp.nb", values_from = "bp.coord") %>%
        print()
      
    }
  }
}

