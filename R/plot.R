## graphical review of simulated datasets

require(ggplot2)
require(dplyr)
require(stringr)

plot.trajData <- function(sim.data, breakpoints = T, lines = T, default = F,
                          cluster = 1:nrow(sim.data$sim.dataset),
                          alpha = .65, true.color = "green4") {
  if(default){
    print.default(sim.data)
  }
  
  traj.data <- sim.data$sim.dataset %>%
    filter(ID %in% cluster)
  
  # deriving pattern for display
  pattern <- traj.data %>% 
    filter(ID==cluster[1]) %>% 
    group_by(segment) %>% 
    summarise(pattern = pattern[1])
  pattern <- pattern$pattern
  pattern <- paste(pattern[!is.na(pattern)], collapse="")
  
  # points customization
  cols = c("signal" = "blue", "outlier" = "red", "trailing" = "black")
  pchs = c("signal" = 19, "outlier" = 17, "trailing" = 1)
  
  # main plot
  p <- ggplot(traj.data, aes(x = time, y = score, color = type, shape = type)) + 
    geom_point() +
    facet_wrap(~ID) +

    # custom y ticks and points style
    scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10)) +
    scale_color_manual(values = cols) + # custom point colors according to measurement type
    scale_shape_manual(values = pchs) + # same but for shape

    # labels
    labs(
      title = paste("Simulated Subjective Drug Intensity over time per patient - pattern",
                    pattern),
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
    # adding breakpoints coord to have proper angles
    p <- p +
      geom_line(
        data = traj.data %>%
          add_row(ID=break.data$ID, 
                  time=break.data$x.coord, 
                  true.traj=break.data$y.coord,
                  type = "breakpoints") %>%
          arrange(ID, time),
        mapping = aes(x = time, y = true.traj), 
        inherit.aes = F, colour = true.color, alpha = alpha
      )
  }
  
  return(p)
}
