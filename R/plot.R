## Method for object of class 'trajTruth'
##' @export
plot.trajTruth <- function(true.traj, breakpoints = T, lines = T,
                           alpha = .65, true.color = "orange2"){
  require(ggplot2)
  require(dplyr)
  
  max.time <- max(true.traj$breakpoints$bp.x)
  pattern <- true.traj$breakpoints$pattern
  pattern <- paste(pattern[!is.na(pattern)], collapse = "")
  
  # ## TODO - theoretical measurements points
  # if (length(true.traj$times)==1){
  #   time <- seq(0, max.time+true.traj$times, by = true.traj$times)
  # } else if (length(true.traj$times>1){
  #   time <- times$value
  # }
  # time.point <- data.frame(
  #   time = c(time, true.traj$breakpoints$bp.x),
  #   seg.lim = c(rep(F, length(time)), rep(T, length(true.traj$breakpoints$bp.x)))
  # ) %>%
  # arrange(time, -seg.lim)
  
  
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
##' @export
plot.trajData <- function(sim.data, breakpoints = T, lines = T, default = F,
                          cluster = levels(sim.data$sim.dataset$ID),
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
    dplyr::select(-starts_with("beta.")) %>%
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


## Method for object of class 'break.lm'
##' @export
plot.break.lm <- function(z, breaks = T, fit = T, default = F,
                          fit.color = "#DF546B", lwd = 2, cex = 1, breaks.ci = F) {
  if (default){
    print.default(z)
  } else {
    x <- z$x
    y <- z$y
    peak.y <- z$peak.y
    psi <- z$psi
    
    plot(x, y, pch = 20, ylim = c(0, max(10, peak.y)),
         main = "101 Model's fit",
         xlab = "Time since drug intake (minutes)",
         ylab = "Subjective Drug Intensity")
    
    if (breaks) points(psi, rep(peak.y, 2), pch = 17, col = fit.color, cex = cex)
    if (fit) {
      # TODO - change when intercept is nonnull (or pattern != 101)
      break.data <- data.frame(
        x = c(0, psi, max(x)),
        y = c(0, rep(peak.y, 2), peak.y + z$coefficients[2]*(max(x) - psi[2]))
      )
      lines(break.data$x, break.data$y, col = fit.color, lwd = lwd)
    }
    if (breaks.ci) {
      warning("Parameter 'breaks.ci' is ignored at the moment")
    }
  }
}
