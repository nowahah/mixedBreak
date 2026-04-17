library(lmbreak)
library(ggplot2)
library(dplyr)

data(SDIpsilo, package = "lmbreak")

SDIpsilo <- SDIpsilo %>% 
  group_by(id) %>%
  mutate(type = as.factor(type),
         score.diff1 = c(diff(score), NA),
         score.diff2 = c(diff(score, 2), NA, NA))


## rough overall summary
summary(SDIpsilo)
# plot parameters
max.time = max(SDIpsilo$time)
cols = c("signal" = "blue", "noise" = "red", "trailing" = "black", "added" = "black")
pchs = c("signal" = 19, "noise" = 17, "trailing" = 1, "added" = 8)

## facet graph of SDIpsilo trajectories
# classic
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

## OBSERVATIONS:
# onset phase: patients' SDI is increasing, nearly strictly
# offset phase: SDI is decreasing nearly strictly. 
#   More "repeated scores" (2+ consecutive identical values). 
#   Trailing may occur during this phase (at least 3 id values).

# diff at lag 1
ggplot(data = SDIpsilo, mapping = aes(x=time, y=score.diff1)) +
  geom_point() +
  facet_wrap(~id, nrow = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color='blue') +
  
  # custom x, y ticks and points shapes and colors
  scale_x_continuous(breaks = seq(0,max.time, by=60), limits = c(0,max.time)) + 
  scale_y_continuous(breaks = seq(-4,4, by=2), limits = c(-4, 4)) + # custom y ticks

  # labels
  labs(title = "Differentiated Subjective Drug Intensity over time for all patients",
       x = "Time (minutes)", y = "SDI")

## OBSERVATIONS:
# onset phase: patients' SDI is positive, nearly strictly (no consecutive equal scores)
# offset phase: SDI is negative, but it is often less clear cut. 
#   More "repeated scores" (2+ consecutive identical values). 
#   Trailing may occur during this phase (2-3 consecutive zeros)
# plateaus: diff is null, up to a single non-zeros value

## ALGORITHM:
# To roughly estimate breakpoints, we could split the trajectory by isolating:
#   - onset phase: longest positive/increasing segment
#   - plateau: longest null/constant segment, up to single non-zero value
#   - offset phase: longest negative/decreasing segment
#   - trailing: same as plateau but from the end of the trajectory


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
  scale_x_continuous(breaks = seq(0,10, by=1)) + # custom x ticks

  # labels
  labs(title = "Subjective Drug Intensity histogram for all patients",
       x = "SDI score", y = "Count")
# ~ U-shaped discrete distribution, bounded between 0-10
# => strong ceiling effect, and strong-moderate floor effect
# Normality assumption is inappropriate, but breaks allows for relaxing hypothesis:
#   they only have to be respected inside a given segment


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
# discrete distribution. The distinction could be made simply by the value of the mode
