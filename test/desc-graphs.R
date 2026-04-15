
## load the datasets with DataLoader file
source("scripts/DataLoader.R")

library(ggplot2)
library(dplyr)

## Basel whole dataset
ggplot(data = basel, mapping = aes(x=time.num, y=score, colour=Study)) +
  geom_point() + facet_wrap(~ID) +
  
  # custom y ticks 
  scale_y_continuous(breaks = seq(0,100,by=25), limits = c(0, 100)) + # custom y ticks
  
  # labels
  labs(title = "Patient's feeling of the experiment over time",
       x = "Time since drug intake (minutes)", y = "Feeling on Visual Analog Scale (0-100)")

## Basel psilo
ggplot(data = basel %>% filter(Study=="SPS"), mapping = aes(x=time.num, y=score, colour=Study)) +
  geom_point() + facet_wrap(~ID) +
  
  # custom x y ticks 
  scale_x_continuous(breaks = seq(0,420,by=60), limits = c(0, 420)) + # custom x ticks for SPS
  scale_y_continuous(breaks = seq(0,100,by=25), limits = c(0, 100)) + # custom y ticks
  
  # labels
  labs(title = "Patient's feeling of the experiment over time",
       x = "Time since drug intake (minutes)", y = "Feeling on Visual Analog Scale (0-100)")


## NRU Psilo
max.time = max(SDIpsilo$time, na.rm=T)
ggplot(data = SDIpsilo %>% filter(!is.na(score) & time>=0),
       mapping = aes(x=time, y=score)) +
  geom_point() +
  facet_wrap(~ID) +
  
  # custom x, y ticks
  scale_x_continuous(breaks = seq(0,max.time, by=60), limits = c(0,max.time)) + 
  scale_y_continuous(breaks = seq(0,10,by=2), limits = c(0, 10)) + # custom y ticks
  
  # labels
  labs(title = "Subjective Drug Intensity over time for all patients",
       x = "Time since drug intake (minutes)", y = "SDI")
