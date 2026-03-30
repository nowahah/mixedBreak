### simData-test.R ----

source("/indirect/student/noahdegoulange/scripts/simData.R")

library(ggplot2)

n.obs <- 4L
breakpoints <- data.frame(
  bp = c(90, 180),
  bp.sd = c(10, 10)
)
b.onset <- list("value" = 2*1.4, "sd" = 2*0.5, "min" = 0.5) # onset coefficient per 20 min
b.return <- list("value" = 2*-0.7, "sd" = 2*0.3, "max" = -0.2) # return to normal coefficient per 20 min
score.sd <- 1 # noise level on the measurements
# different time specifications
time.reg <- list("value" = 20, "sd" = 0) # default regular "perfect" measurement strategy
time.noise <- list("value" = 20, "sd" = 2) # regular "perfect" measurement strategy with noise
time.spec <- list("value" = c(0,10,20,30,40,60,120), "sd" = 1)   # specific measurement strategy

outlier.prob <- 0.05
na.prob <- 0.10
n.trail <- 0L


sim.dataset <- simData(n.obs, breakpoints, b.onset, b.return, score.sd,
                       times = time.noise,
                       outlier.prob = outlier.prob, na.prob = na.prob,
                       n.trail = n.trail)
summary(sim.measurements <- sim.dataset %>% select(-c(3:8)))
(sim.patients <- sim.dataset %>% select(c(1, 3:8)) %>% distinct())


## graphical review
# TODO add factor variable to code and display normal / outlier / missing
ggplot(sim.dataset %>% mutate(label = paste("ID:", ID, "- Peak = ", round(peak, 1))), 
       aes(x=time, y=score, colour=outliers)) +
  geom_point() + 
  facet_wrap(~label) +
  
  # adding breakpoints and lines
  geom_point(aes(x=break.1, y=0), colour = "black", shape = 18, size=3) +
  geom_vline(aes(xintercept = break.1), linetype = "dashed", colour = "grey50") +
  geom_point(aes(x=break.2, y=0), colour = "black", shape = 18, size=3) +
  geom_vline(aes(xintercept = break.2), linetype = "dashed", colour = "grey50") +
  geom_point(aes(x=return.0, y=0), colour = "black", shape = 3, size=2) +
  
  
  # custom y ticks
  scale_y_continuous(breaks = seq(0,10,by=2), limits = c(0, 10)) + 
  
  # labels
  labs(title = "Simulated Subjective Drug Intensity over time for all patients",
       x = "Time since drug intake (minutes)", y = "SDI")


# It works but onset relation ship is non monotonic VS real ones tends to be during a given phase
# Higher variance during plateau (especially if peak <≈ 10) -> not realistic (need lower variance during plateau)
