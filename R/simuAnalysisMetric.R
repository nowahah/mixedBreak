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


## 2. Performance measure and Monte Carlo SE ====
# 2.1. segmented.lme


# 2.2. segreg
