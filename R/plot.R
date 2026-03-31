## graphical review of simulated datasets

library(ggplot2)
library(dplyr)

plot.trajData <- function(sim.dataset){
  # TODO add factor variable to code and display normal / outlier / missing
  ggplot(sim.dataset %>% mutate(label = paste("ID:", ID, "- Plateau = ", round(plateau, 1))), 
         aes(x=time, y=score, colour=outliers)) +
    geom_point() + 
    facet_wrap(~label) +
    
    # adding breakpoints and lines
    geom_point(aes(x=break.1, y=0), colour = "black", shape = 18, size=3) +
    geom_vline(aes(xintercept = break.1), linetype = "dashed", colour = "grey50") +
    geom_point(aes(x=break.2, y=0), colour = "black", shape = 18, size=3) +
    geom_vline(aes(xintercept = break.2), linetype = "dashed", colour = "grey50") +
    geom_point(aes(x=return.0, y=0), colour = "black", shape = 3, size=2) +
    
    
    # custom y ticks
    scale_y_continuous(breaks = seq(0,10,by=2), limits = c(0, 10)) + 
    
    # labels
    labs(title = "Simulated Subjective Drug Intensity over time for all patients",
         x = "Time since drug intake (minutes)", y = "SDI")
  
}