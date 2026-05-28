library(dplyr)
library(ggplot2)
library(rsimsum)
library(tidyr)

## Reading results in files
library(readxl)
estimates.seglme <- read_excel("data/esimates_seglme.xlsx")
estimates.segreg <- read_excel("data/esimates_segreg.xlsx", guess_max = 10000L)
true.parameters <- read_excel("data/true_parameters.xlsx", sheet = "true.parameters")
scenario.data <- read_excel("data/true_parameters.xlsx", sheet = "scenario.data")
all.seeds <- read_excel("data/random_seed.xlsx", sheet = "all.seeds")$allseeds
# ending.seed <- read_excel("data/random_seed.xlsx", sheet = "ending.seed")$ending.seed
# restore with .GlobalEnv$.Random.seed <- ending.seed
TRUE.PSI <- 77

ref.ind <- which(scenario.data$fac.var=="REFERENCE")
scenario.data$fac.var[ref.ind] <- scenario.data$fac.var[ref.ind-1]
scenarLabel <- character(nrow(scenario.data))
for(i in 1:nrow(scenario.data)){
  scenarLabel[i] <- paste(scenario.data$fac.var[i],"=", 
                          scenario.data[i, scenario.data$fac.var[i]])
}
scenarLabel[ref.ind] <- paste(scenarLabel[ref.ind], "(REF)")
reorder.index <- c(2,1,12,11,13,14,9,10, 3,5,6,7,8,4)
scenario.data$scenarLabel <- as.factor(scenarLabel)
scenario.data$scenarLabel <- factor(scenario.data$scenarLabel, 
                                    levels(scenario.data$scenarLabel)[reorder.index])
levels(scenario.data$scenarLabel) <- paste0("(", 1:nrow(scenario.data), ") ", 
                                            levels(factor(scenario.data$scenarLabel)))


## creation of unique dataset with all results for both methods
seglme.res <- estimates.seglme %>%
  group_by(simID) %>%
  summarise(break.avg = first(break.avg),
            break.x1.sd = first(break.x1.sd),
            break.x1.random = first(break.x1.random),
            break.CI95.low = first(break.CI95.low),
            break.CI95.up = first(break.CI95.up))
segreg.res <- estimates.segreg %>% 
  group_by(simID) %>%
  summarise(break.avg = mean(break.x1),
            break.x1.sd = sd(break.x1),
            break.CI95.low = mean(break.CI95.low),
            break.CI95.up = mean(break.CI95.up))
dataset <- seglme.res %>%
  bind_rows(segreg.res, .id = "method") %>%
  mutate(method = if_else(method==1, "segmented.lme", "segreg")) %>%
  left_join(true.parameters %>% 
              select(simID, scenarID) %>% 
              distinct(), 
            by = "simID") %>%
  left_join(scenario.data %>%
              select(scenarID, scenarLabel), 
            by= "scenarID")

metrics <- simsum(
  dataset, 
  estvarname = "break.avg", 
  se = "break.x1.sd",
  true = TRUE.PSI, # true value of estimated parameter 
  methodvar = "method", # variable containing the methods to compare
  by = "scenarLabel", # factor of data-generating mechanism
  ci.limits = c("break.CI95.low", "break.CI95.up"),
  x = TRUE, # necessary for zip plots
  control = list(mcse = TRUE)
)

# every metrics in summary
summ <- summary(metrics, stats = c("thetamean", "se2mean", "bias", "empse", 
                                   "cover", "becover"))
# adding publication-ready columns with est (mcse)
formatted.res <- summ$summ %>% 
  select(-c("lower", "upper")) %>%
  filter(stat %in% c("bias", "empse", "cover")) %>%
  pivot_wider(names_from = method, values_from = c(est:mcse)) %>%
  arrange(stat) %>%
  mutate(round.lme = -floor(log10(abs(mcse_segmented.lme))),
         segmented.lme = paste0(round(est_segmented.lme, round.lme), " (", 
                                round(mcse_segmented.lme, round.lme), ")"),
         round.seg = -floor(log10(abs(mcse_segreg))),
         segreg = paste0(round(est_segreg, round.seg), " (", 
                         round(mcse_segreg, round.seg), ")")) %>%
  select(stat, scenarLabel, segmented.lme, segreg)

formatted.res %>%
  rename(seglme = segmented.lme) %>%
  pivot_wider(names_from = stat, values_from = c("seglme", "segreg"))%>%
  gt(
    rowname_col = "scenarLabel",
    row_group_as_column = TRUE
  ) %>%
  tab_header(title = "Performance measurements per method and DGM") %>%
  tab_stubhead(label = md("**DGM**")) %>%
  cols_align(align = "left", columns = scenarLabel) %>%
  tab_spanner(label=md("**Bias**"), columns=ends_with("_bias")) %>%
  tab_spanner(label=md("**Coverage**"), columns=ends_with("_cover")) %>%
  tab_spanner(label=md("**EmpSE**"), columns=ends_with("_empse"))


# bias lollipop plot
autoplot(summ, stats = "bias", type = "lolly") # GOOD
autoplot(summ, stats = "cover", type = "lolly") 
autoplot(summ, stats = "empse", type = "lolly") 