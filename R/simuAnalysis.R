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

## 0. Missing values - reported at the study level for both methods ====
# 0.1. segmented.lme
estimates.seglme %>% 
  select(simID, error) %>% 
  distinct() %>% 
  summarise(error.rate = sum(!is.na(error))/n(),
            error.nb = sum(!is.na(error)))
# prop of missing estimates is higher than expected 1351/7826 (17.3% > 15%)
# => need to extend the simulation. Restore the ending seed with
# .GlobalEnv$.Random.seed <- ending.seed

# what kind of errors ?
error.lme <- estimates.seglme %>% 
  select(simID, error) %>% 
  distinct() %>%
  select(error)
error.lme <- as.factor(error.lme$error)
levels(error.lme) # aggregate identical error (up to a certain irrelevant value)
levels(error.lme) <- c("nlminb problem, singular convergence (7)", rep(substr(levels(error.lme)[2], 0, 34), 4))
table(error.lme)

# 0.2. segreg
estimates.segreg %>% 
  select(simID, error) %>% 
  distinct() %>% 
  summarise(error.rate = sum(!is.na(error))/n(),
            error.nb = sum(!is.na(error)))
# good 202/7826 (2.6%)

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
  "Not enought information for the fit",
  "psi values too close each other",
  "starting psi out of the range"
)
table(error.segreg2) # uh ? I see at least two errors that could have been prevented

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

# 0.3.0. Explore datasets from which an error has occurred
error.seed <- all.seeds[error.distrib$simID]


## 1. Graphical explanatory analysis ====
# 1.1. Distribution of theta and sd.theta - for each method and data generating model
# LOOK FOR: outliers, overall distribution of data

# 1.2. Distribution of theta VS sd.theta - for each method and data generating model
# LOOK FOR: bivariate outliers

# 1.3. Bivariate plot of theta
# LOOK FOR: correlation between methods, systematic differences

# 1.4. Bland-Altmann plot - with comparator being the true value - for each methods
# LOOK FOR: agreement between methods and truth

# 1.5. Zip plot: CI ranked by |zi|=|thetai-theta|/sd.theta
# LOOK FOR: issue with coverage


## 2. Performance measure and Monte Carlo SE ====
# 2.1. segmented.lme


# 2.2. segreg
