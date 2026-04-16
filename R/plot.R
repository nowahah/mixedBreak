
plot <- function(object){
  UseMethod("plot")
}

## Method for object of class 'trajTruth'
plot.trajTruth <- function(true.traj, breakpoints = T, lines = T,
                           alpha = .65, true.color = "orange2"){
  require(ggplot2)
  
  max.time <- max(true.traj$breakpoints$bp.x)
  pattern <- true.traj$breakpoints$pattern
  pattern <- paste(pattern[!is.na(pattern)], collapse = "")
  ggplot(true.traj$breakpoints, mapping = aes(x=bp.x, y=bp.y)) +
    geom_line(colour = true.color, alpha = alpha) + 
    geom_point(colour = true.color, alpha = alpha, shape = 18, size = 3) +
    
    # custom y ticks and points style
    scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
    scale_x_continuous(breaks = seq(0, max.time, by = 60), limits = c(0, max.time)) +
    
    labs(
      title = paste("Pop.-level 'True' Subjective Drug Intensity over time - pattern", pattern),
      x = "Time since drug intake (minutes)", y = "'True' SDI"
    )
}


## Method for object of class 'trajData'
# graphical review of simulated dataset
plot.trajData <- function(sim.data, breakpoints = T, lines = T, default = F,
                          cluster = 1:nrow(sim.data$sim.dataset),
                          alpha = .65, true.color = "green4") {
  require(ggplot2)
  require(dplyr)
  require(stringr)
  
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
