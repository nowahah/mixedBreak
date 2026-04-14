# library(remotes)
# install_github("bozenne/lmbreak")

library(lmbreak)   # SDIpsilo data
library(segmented) # Segmented Relationships in Regression Models with Breakpoints
library(dplyr)     # data manipulation
library(ggplot2)   # data visualization

data(SDIpsilo, package = "lmbreak")
# SDIpsilo <- SDIpsilo[SDIpsilo$type %in% c("noise","trailing") == FALSE,]
str(SDIpsilo)

# general plot parameters
max.time = max(SDIpsilo$time)
cols = c("signal" = "blue", "noise" = "red", "trailing" = "black", "added" = "black")
pchs = c("signal" = 19, "noise" = 17, "trailing" = 1, "added" = 8)

##################################
## focus example on one patient ##
##################################
patient = 2
SDIpsilo.ind <- SDIpsilo %>% filter(id==patient & !is.na(score))
out.lm <- lm(score ~ 0 + time, data = SDIpsilo.ind) 

# starting model includes just 1 covariate: time
n.breaks = 2
o<-segmented(out.lm, seg.Z = ~time, npsi = n.breaks) # 2 breakpoints for time
summary(o)
(bp.ci = confint.segmented(o))

# plot the output of model
p <- ggplot(SDIpsilo.ind, aes(x=time, y=score, colour = type, shape = type)) +
  geom_point() +
  
  # custom x, y ticks and points shapes and colors
  scale_x_continuous(breaks = seq(0,max.time, by=40), limits = c(0,max.time)) + 
  scale_y_continuous(breaks = 0:10, limits = c(0, 10)) + # custom y ticks
  scale_color_manual(values = cols) + # custom point colors according to measurment type
  scale_shape_manual(values = pchs) + # same but for shape
  
  # labels
  labs(title = paste("Subjective Drug Intensity over time for patient:", patient),
       x = "Time (minutes)", y = "SDI")


# Adding breakpoints and 95% confidence intervals
p <- p +
  geom_point(aes(x=bp.ci[1,1], y=0), colour = "black", shape = 18, size=3) + # bp1
  geom_vline(xintercept = bp.ci[1,1], linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x=bp.ci[1,2], y=0, xend=bp.ci[1,3], yend=0), colour = "black") +
  geom_point(aes(x=bp.ci[2,1], y=0), colour = "black", shape = 18, size=3) + # bp2
  geom_vline(xintercept = bp.ci[2,1], linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x=bp.ci[2,2], y=0, xend=bp.ci[2,3], yend=0), colour = "black") 
  # geom_point(aes(x=bp.ci[3,1], y=0), colour = "black", shape = 18, size=3) + # bp3
  # geom_vline(xintercept = bp.ci[3,1], linetype = "dashed", colour = "grey50") +
  # geom_segment(aes(x=bp.ci[3,2], y=0, xend=bp.ci[3,3], yend=0), colour = "black") 

p

# Adding model's fit and confidence intervals
# TODO for CI
p +
  geom_segment(aes(x=bp.ci[1,1], y=predict(o, data.frame(time = bp.ci[1,1])),
                   xend=bp.ci[2,1], yend=predict(o, data.frame(time = bp.ci[2,1]))))
  # geom_point(aes(x=bp.ci[1,1], y=predict(o, data.frame(time = bp.ci[1,1]))), 
  #            colour = "orange", size=3) + # bp1 ?
  # geom_point(aes(x=time, y=o$fitted.values), 
  #            shape = 5, size = 2, color = "#008B45") +
  # geom_line(aes(x=time, y=o$fitted.values), color = "#008B45")

# segmented plot
plot.segmented(o, conf.level=.95, res=T, ylim = c(0,12), 
               pch = pchs[SDIpsilo.ind$type],
               res.col = cols[SDIpsilo.ind$type],
               ylab = "SDI", xlab = "Time (minutes)", 
               main = paste("Subjective Drug Intensity over time for patient:", patient))


# # residuals over time
# plot(SDIpsilo.ind$time[!is.na(SDIpsilo.ind$score)], o$residuals,
#      xlab = "Time (minutes)", ylab = "Residuals", 
#      main = "Residuals from model against time")
# # residuals against fitted values
# plot(o$fitted.values, o$residuals,
#      xlab = "Fitted values (SDI)", ylab = "Residuals", 
#      main = "Residuals against fitted values")


# facet graph of SDIpsilo data
ggplot(data = SDIpsilo, mapping = aes(x=time, y=score, colour=type, shape=type)) +
  geom_point() +
  facet_wrap(~id, nrow = 4) +
  
  # custom x, y ticks and points shapes and colors
  scale_x_continuous(breaks = seq(0,max.time, by=40), limits = c(0,max.time)) + 
  scale_y_continuous(breaks = 0:10, limits = c(0, 10)) + # custom y ticks
  scale_color_manual(values = cols) + # custom point colors according to measurment type
  scale_shape_manual(values = pchs) + # same but for shape
  
  # labels
  labs(title = "Subjective Drug Intensity over time for all patients",
       x = "Time (minutes)", y = "SDI")
