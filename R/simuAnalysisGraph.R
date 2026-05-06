library(dplyr)
library(ggplot2)


## Reading results in files
library(readxl)
estimates.seglme <- read_excel("data/esimates_seglme.xlsx")
estimates.segreg <- read_excel("data/esimates_segreg.xlsx", guess_max = 10000L)
true.parameters <- read_excel("data/true_parameters.xlsx", sheet = "true.parameters")
scenario.data <- read_excel("data/true_parameters.xlsx", sheet = "scenario.data")
all.seeds <- read_excel("data/random_seed.xlsx", sheet = "all.seeds")$allseeds
# ending.seed <- read_excel("data/random_seed.xlsx", sheet = "ending.seed")$ending.seed
# restore with .GlobalEnv$.Random.seed <- ending.seed

ref.ind <- which(scenario.data$fac.var=="REFERENCE")
scenario.data$fac.var[ref.ind] <- scenario.data$fac.var[ref.ind-1]
scenarLabel <- character(nrow(scenario.data))
for(i in 1:nrow(scenario.data)){
  scenarLabel[i] <- paste(scenario.data$fac.var[i],"=", scenario.data[i, scenario.data$fac.var[i]])
}
scenarLabel[ref.ind] <- paste(scenarLabel[ref.ind], "(REF)")# should I ?
scenario.data$scenarLabel <- scenarLabel

## 1. Graphical explanatory analysis ====
# 1.1. Distribution of theta and sd.theta - for each method and data generating model
# LOOK FOR: outliers, overall distribution of data

# creating appropriate dataset with labels for the plot
dataset <- estimates.seglme %>% 
  select(simID, break.avg) %>% 
  distinct() %>% 
  left_join(true.parameters[, c("simID", "scenarID")], by="simID") %>% 
  left_join(scenario.data[, c("scenarID", "scenarLabel")], by="scenarID") %>%
  distinct()
summary(dataset)
ggplot(dataset , aes(x=break.avg)) +
  geom_boxplot() +
  facet_wrap(vars(scenarLabel)) +
  geom_vline(xintercept = 77, colour = "orange2", linetype = "dashed") +
  labs(title="Breakpoint.x distribution for every data-generating mechanism") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
# order is a bit odd

# 1.2. Bivariate plot of theta
# LOOK FOR: correlation between methods, systematic differences


# 1.3. Zip plot: CI ranked by |zi|=|thetai-theta|/sd.theta
# LOOK FOR: issue with coverage

