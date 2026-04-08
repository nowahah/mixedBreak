## graphical review of simulated datasets

require(ggplot2)
require(dplyr)
require(stringr)

plot.trajData <- function(sim.data, breakpoints = T, lines = T, default = F,
                          cluster = 1:nrow(sim.data$sim.dataset),
                          alpha = .75, true.color = "green4") {
  if(default){
    print.default(sim.data)
  }
  
  traj.data <- sim.data$sim.dataset %>%
    filter(ID %in% cluster)
  
  # TODO add factor variable to code and display normal / outlier / missing
  p <- ggplot(traj.data, aes(x = time, y = score)) + # , color=outliers
    geom_point() +
    facet_wrap(~ID) +

    # custom y ticks
    scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10)) +

    # labels
    labs(
      title = "Simulated Subjective Drug Intensity over time for all patients",
      x = "Time since drug intake (minutes)", y = "SDI"
    )
    
  # breakpoints coordinates
  break.data <- sim.data$sim.gen.model %>%
    filter(ID %in% cluster) %>%
    select(-starts_with("beta.")) %>%
    pivot_longer(cols = starts_with("break.x"), names_to = "break.x", values_to = "x.coord") %>%
    pivot_longer(cols = starts_with("break.y"), names_to = "break.y", values_to = "y.coord") %>%
    filter(str_sub(break.x, -1L) == str_sub(break.y, -1L))

  # Add Breakpoints on the graph
  if(breakpoints){
    
    p <- p +
      geom_point(
        data = break.data, mapping = aes(x = x.coord, y = y.coord),
        colour = true.color, shape = 18, size = 3, alpha = alpha
      )
  }  
  
  # Add the true trajectory
  if(lines){
    # TODO - add breakpoints in the data to have proper angles at breaks
    p <- p +
      geom_line(
        data = traj.data %>%
          add_row(ID=break.data$ID, time=break.data$x.coord, true.traj=break.data$y.coord) %>%
          arrange(ID, time),
        mapping = aes(x = time, y = true.traj), colour = true.color, alpha = alpha
      )
  }
  
  return(p)
}
