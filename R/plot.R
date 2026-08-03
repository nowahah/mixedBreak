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
    psi <- z$psi.i
    
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


## Method for object of class 'mixedBreak1'
##' @export
plot.mixedBreak1 <- function(
    z, breaks = TRUE, fit = TRUE, fit.color = "orange2", lwd = 2, cex = 1, 
    breaks.ci = FALSE, alpha = .65, y.lim = c(NA, NA)
) {
  #TODO - check the script especially call to `[[`
  # require(ggplot2)
  require(dplyr)
  
  pattern <- z$pattern
  vars <- z$var.name
  model.plot <- z$model[,unique(c(vars$response, vars$segmented, vars$group))]
  names(model.plot)[names(model.plot)==vars$group] <- "ID"
  names(model.plot)[names(model.plot)==vars$response] <- "yy"
  model.plot$fitted <- z$fitted
  model.plot$psi <- z$psi.i[model.plot$ID]
  model.plot$psi.y <- (z$psi.i * z$random[[1]])[model.plot$ID]
  
  p <- ggplot2::ggplot(model.plot, ggplot2::aes(x=time, y=yy)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ID) + 
    # ggplot2::scale_y_continuous(breaks = seq(0,14,by=2), limits = y.lim) +
    ggplot2::xlab("Time since drug intake (minutes)") +
    ggplot2::ylab("SDI") 
  
  fit.data <- data.frame(
    ID = factor(levels(model.plot$ID), levels = levels(model.plot$ID)),
    psi = z$psi.i,
    psi.y = z$psi.i * z$random[[1]],
    max.time = (model.plot %>% 
                  group_by(ID) %>% 
                  summarise(max.time = max(time)))$max.time
  )
  fit.data <- fit.data %>%
    mutate(
      beta.2 = if_else(
        rep(pattern=="11", nlevels(model.plot$ID)),
        z$random$U + z$random[[1]], 0
      ),
      yend = psi.y + (max.time-psi)*beta.2)

  if (breaks) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(x = psi, y = psi.y), data = fit.data,
                 colour = fit.color, shape = 18, size = 3, alpha = alpha) +
      ggplot2::annotate(ggplot2::GeomPoint, x = 0, y = 0, 
               colour = fit.color, shape = 18, size = 3, alpha = alpha) +
      ggplot2::geom_point(ggplot2::aes(x = max.time, y = yend), data = fit.data,
                 colour = fit.color, shape = 18, size = 3, alpha = alpha)
    
  }
  if (fit) {
    p <- p +
      ggplot2::geom_segment(ggplot2::aes(x=0, y=0, xend=psi, yend=psi.y), 
                            data = fit.data, colour = fit.color, lwd = 1, alpha = alpha)  +
      ggplot2::geom_segment(ggplot2::aes(x=psi, y=psi.y, xend=max.time, yend=yend), 
                            data = fit.data, colour = fit.color, lwd = 1, alpha = alpha)
  }
  if (breaks.ci) {
    warning("Parameter 'breaks.ci' is ignored at the moment")
  }
  
  return(p)
}


## Method for object of class 'mixedBreak2'
##' @export
##' # Careful here, when pattern == "101" and breakpoints are switched,
##' The plot has to be adapted
plot.mixedBreak2 <- function(
    z, breaks = TRUE, fit = TRUE, fit.color = "orange2", lwd = 2, cex = 1, 
    breaks.ci = FALSE, alpha = .75, breaks.vline = FALSE, y.lim = c(NA, NA)
    #, cluster = NA, 
) {
  # require(ggplot2)
  require(dplyr) # operator `%>%`
  
  if(any(diff(t(z$psi.i)) < 0)) {
    # browser()
    # ind.switched <- which(diff(t(z$psi.i)) < 0)
    ind.switched <- rownames(z$psi.i)[diff(t(z$psi.i)) < 0]
    warning(paste(
      "Some individuals have switched estimated breakpoints:\n  ",
      z$var.name$group, "= {", paste(ind.switched, collapse = ", "), "}\n  ",
      "Graph is incorrect for these groups."
    ))
  }
  
  # browser()
  pattern <- z$pattern
  n.psi <- dim(z$psi.history)[2]
  vars <- z$var.name
  model.plot <- z$model[,unique(c(vars$response, vars$segmented, vars$group))]
  names(model.plot)[names(model.plot)==vars$group] <- "ID"
  names(model.plot)[names(model.plot)==vars$response] <- "yy"
  names(model.plot)[names(model.plot)==vars$segmented] <- "time"
  model.plot$fitted <- z$fitted
  model.plot$psi <- z$psi.i[model.plot$ID, ]
  model.plot$psi.y1 <- (z$psi.i[,1] * z$random[[1]])[model.plot$ID]
  # if(!all(is.na(cluster))){
  #   model.plot <- model.plot[model.plot[[vars$group]] %in% cluster, ]
  # }
  
  p <- ggplot2::ggplot(model.plot, ggplot2::aes(x=time, y=yy)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ID) +
    ggplot2::xlab("Time since drug intake (minutes)") +
    ggplot2::ylab("SDI") + 
    # ggplot2::scale_y_continuous(breaks = seq(-25, 100, by=25), limits = y.lim) +
    ggplot2::geom_hline(yintercept = 0, lty = "dashed")
  
  fit.data <- data.frame(
    ID = factor(levels(model.plot$ID), levels = levels(model.plot$ID)),
    psi = unname(z$psi.i),
    psi.y1 = z$psi.i[,1] * z$random[[1]], # there are psi.y1 and psi.y2
    max.time = (model.plot %>% group_by(ID) %>% 
                  summarise(max.time = max(time)))$max.time,
    beta.2 = if_else(
      rep(pattern=="111", nlevels(model.plot$ID)),
      z$random[[1]] + z$random$U1, 0
    )
  )
  # browser()
  fit.data <- fit.data %>%
    dplyr::mutate(
      psi.y2 = psi.y1 + beta.2*(psi.2 - psi.1),
      beta.3 = beta.2 + z$random$U2,
      yend = psi.y2 + beta.3*(max.time-psi.2)
    )
  # browser()

  if (breaks) {
    p <- p +
      ggplot2::annotate(ggplot2::GeomPoint, x = 0, y = 0, 
               colour = fit.color, shape = 18, size = 3, alpha = alpha) +
      ggplot2::geom_point(ggplot2::aes(x = psi.1, y = psi.y1), data = fit.data,
                 colour = fit.color, shape = 18, size = 3, alpha = alpha) + # square: 2
      ggplot2::geom_point(ggplot2::aes(x = psi.2, y = psi.y2), data = fit.data,
                 colour = fit.color, shape = 18, size = 3, alpha = alpha) + # triangle?
      ggplot2::geom_point(ggplot2::aes(x = max.time, y = yend), data = fit.data,
                 colour = fit.color, shape = 18, size = 3, alpha = alpha)
    
  }
  if (fit) {
    p <- p +
      ggplot2::geom_segment(ggplot2::aes(x=0, y=0, xend=psi.1, yend=psi.y1), 
                            data = fit.data, colour = fit.color, lwd = 1, alpha = alpha) +
      ggplot2::geom_segment(ggplot2::aes(x=psi.1, y=psi.y1, xend=psi.2, yend=psi.y2), 
                            data = fit.data, colour = fit.color, lwd = 1, alpha = alpha) +
      ggplot2::geom_segment(ggplot2::aes(x=psi.2, y=psi.y2, xend=max.time, yend=yend), 
                            data = fit.data, colour = fit.color, lwd = 1, alpha = alpha)
  }
  if (breaks.ci) {
    warning("Parameter 'breaks.ci' is ignored at the moment")
  }
  
  p
}