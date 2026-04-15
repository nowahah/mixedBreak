library(lmbreak)
library(ggplot2)

data(SDIpsilo, package = "lmbreak")

SDIpsilo$type <- as.factor(SDIpsilo$type)

## rough overall summary
summary(SDIpsilo)

## facet graph of SDIpsilo trajectories
ggplot(data = SDIpsilo, mapping = aes(x=time, y=score, colour=type, shape=type)) +
  geom_point() +
  facet_wrap(~id, nrow = 4) +
  
  # custom x, y ticks and points shapes and colors
  scale_x_continuous(breaks = seq(0,max.time, by=60), limits = c(0,max.time)) + 
  scale_y_continuous(breaks = seq(0,10, by=2), limits = c(0, 10)) + # custom y ticks
  scale_color_manual(values = cols) + # custom point colors according to measurement type
  scale_shape_manual(values = pchs) + # same but for shape
  
  # labels
  labs(title = "Subjective Drug Intensity over time for all patients",
       x = "Time (minutes)", y = "SDI")

## table and histogram of SDI scores
table(SDIpsilo$score) # 3.5 ??
ggplot(SDIpsilo, aes(x=score)) + 
  geom_histogram(binwidth=.5) +
  
  # adding counts on top
  stat_count(binwidth = 1, 
             geom = 'text', 
             color = 'white', 
             aes(label = after_stat(count)),
             position = position_stack(vjust = 0.5)) +
  
  # custom x ticks 
  scale_x_continuous(breaks = seq(0,10, by=1), limits = c(-.5, 10.5)) + # custom x ticks

  # labels
  labs(title = "Subjective Drug Intensity histogram for all patients",
       x = "SDI score", y = "Count")
# I wouldn't call it normal... 
# ~ U-shaped discrete distribution, bounded between 0-10


## histogram of SDI scores per patient
ggplot(SDIpsilo, aes(x=score)) + 
  geom_histogram(binwidth=.5) +
  facet_wrap(~id) +
  
  # custom x ticks 
  scale_x_continuous(breaks = seq(0,10, by=1), limits = c(-.5, 10.5)) + # custom x ticks
  
  # labels
  labs(title = "Subjective Drug Intensity histogram for every patients",
       x = "SDI score", y = "Count")
# most patient we could assume they either peak of trail at the mode of their 
# discrete distribution. The distinction could be mad simply by the value of the mode
