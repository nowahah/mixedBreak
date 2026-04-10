library(dplyr)

print.trajData <- function(sim.data, default = FALSE, dgm = T, n.lines = 10L,
                           cluster = levels(sim.data$sim.dataset$ID)){
  if(default){
    print.data.frame(sim.data)
  }else{
    traj.data <- sim.data$sim.dataset %>%
      filter(ID %in% cluster)
    
    # Prints simulated data
    if(n.lines>0L){
      cat(paste0("\nSimulated data head(", n.lines,"):\n"))
      print.data.frame(traj.data %>% head(n.lines))
    }
    
    # Prints Data Generation Model
    if(dgm){
      pattern <- traj.data %>% 
        filter(ID==cluster[1]) %>% 
        group_by(segment) %>% 
        summarise(pattern = pattern[1])
      pattern <- pattern$pattern[!is.na(pattern)]
      pattern <- paste(pattern, collapse="")
      cat(paste("\nData Generation Model - pattern", pattern, "\n"))
      
      # breakpoints coordinates 
      break.data <- sim.data$sim.gen.model %>%
        select(-starts_with("beta.")) %>%
        pivot_longer(cols = starts_with("break.x"), names_to = "break.x", values_to = "x.coord") %>%
        pivot_longer(cols = starts_with("break.y"), names_to = "break.y", values_to = "y.coord") %>%
        filter(str_sub(break.x, -1L) == str_sub(break.y, -1L))
      
      # printed result
      break.data %>%
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

