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

## 0. Missing values - reported at the study level for both methods ====
# 0.1. segmented.lme
estimates.seglme %>% 
  select(simID, error) %>% 
  distinct() %>% 
  summarise(error.rate = sum(!is.na(error))/n(),
            error.nb = sum(!is.na(error)))
# rate of 1357/8316 (16.3%), but under allocated threshold of 20%, so it's OK

# what kind of errors ?
error.lme <- estimates.seglme %>% 
  select(simID, error) %>% 
  distinct() %>%
  select(error)
error.lme <- as.factor(error.lme$error)
levels(error.lme) # aggregate identical error (up to a certain irrelevant value)
levels(error.lme) <- c("nlminb problem, singular convergence (7)", 
                       rep(substr(levels(error.lme)[2], 0, 34), nlevels(error.lme)-1))
table(error.lme)

# 0.2. segreg
estimates.segreg %>% 
  select(simID, error) %>% 
  distinct() %>% 
  summarise(error.rate = sum(!is.na(error))/n(),
            error.nb = sum(!is.na(error)))
# good 219/8316 (2.6%)

# what kind of error
error.segreg <- estimates.segreg %>% 
  select(simID, error) %>% 
  distinct() %>%
  select(error)
error.segreg <- as.factor(error.segreg$error)
levels(error.segreg) # truncated ID for which error occurred
error.segreg2 <- factor(sub(".*==> ", "", error.segreg))
levels(error.segreg2) <- c(
  "$ operator is invalid for atomic vectors",
  "Not enough information for the fit: too many psi/small sample/replicated data",
  "psi values too close each other. Please change (decreases number of) starting values",
  "starting psi out of the range"
)
table(error.segreg2) 
# At least two errors that could have been prevented with a robust code:
# 1. $ operator is invalid for atomic vectors
# 2. starting psi outside of the range

# which are common to both methods ?
error.distrib <- estimates.seglme %>%
  select(simID, error) %>% 
  distinct() %>% 
  inner_join(estimates.segreg %>% select(simID, error) %>% distinct(), by = join_by(simID)) %>%
  mutate(error.common = case_when(
    !is.na(error.x) & !is.na(error.y) ~ "both",
    !is.na(error.x) ~ "segmented.lme",
    !is.na(error.y) ~ "segreg",
    .default = "None"
  )) %>% 
  filter(error.common!="None")
table(error.distrib$error.common)
table(error.segreg2, error.lme, useNA = "ifany")
addmargins(table(error.segreg2, error.lme, useNA = "ifany"))
# none of computationally singular was shared with segreg
# most interestingly, the vast majority of the mixed model errors did not affect
# the patient-wise model

# 0.3.0. Explore datasets from which an error has occurred
error.seed <- all.seeds[error.distrib$simID]