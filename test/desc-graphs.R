
## load the datasets with DataLoader file
source("test/DataLoader.R")

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
max.time <- basel %>% filter(Study=="SPS" & !is.na(score)) %>% summarise(max(time.num))
ggplot(data = basel %>% filter(Study=="SPS"), mapping = aes(x=time.num, y=score, colour=Study)) +
  geom_point() + facet_wrap(~ID) +
  
  # custom x y ticks 
  scale_x_continuous(breaks = seq(0,max.time[[1]],by=60), limits = c(0, max.time[[1]])) + # custom x ticks for SPS
  scale_y_continuous(breaks = seq(0,100,by=25), limits = c(0, 100)) + # custom y ticks
  
  # labels
  labs(title = "Patient's feeling of the experiment over time",
       x = "Time since drug intake (minutes)", y = "Feeling on Visual Analog Scale (0-100)")


## NRU psilo ???
head(SDIpsilo.ext)

## NRU Psilo
data(SDIpsilo, package = "lmbreak")
SDIpsilo <- SDIpsilo[SDIpsilo$type %in% c("noise","trailing") == FALSE,] # article's analysis
max.time = max(SDIpsilo$time, na.rm=T)
cols <- c("signal" = "blue", "noise" = "red", "trailing" = "black", "added" = "black")
pchs <- c("signal" = 19, "noise" = 17, "trailing" = 1, "added" = 8)
cluster <- levels(SDIpsilo$id) # c(1, 9, 13)
ggplot(data = SDIpsilo %>% filter(!is.na(score) & id %in% cluster),
       mapping = aes(x=time, y=score, colour = type, shape = type)) + # 
  geom_point() +
  facet_wrap(~id) +
  
  # custom x, y ticks
  scale_x_continuous(breaks = seq(0,max.time, by=120), limits = c(0,max.time)) + 
  scale_y_continuous(breaks = seq(0,10,by=2), limits = c(0, 10)) + # custom y ticks
  
  # custom points
  scale_color_manual(values = cols) + # custom point colors according to measurment type
  scale_shape_manual(values = pchs) + # same but for shape
  
  # labels
  labs(title = "Subjective Drug Intensity over time for all patients") +
  xlab("Time since drug intake (minutes)") + 
  ylab("SDI")

# ggsave("inst/figs/SDI-13-nofit.png", width = 148, height = 105, units = "mm")
