# library(remotes)
# install_github("bozenne/lmbreak")

library(lmbreak) # Brice's package for patient-wise linear regression with breakpoints
library(dplyr) # data manipulation
library(ggplot2) # data visualization

data(SDIpsilo, package = "lmbreak")
SDIpsilo <- SDIpsilo %>% filter(type %in% c("noise", "trailing") == FALSE) # rm noise and trailing
str(SDIpsilo)

# focus example on a single patient
patient <- 2
SDIpsilo.ind <- SDIpsilo %>% filter(id == patient)

# different patterns for the regression
e.XP101 <- lmbreak(score ~ 0 + bp(time, "101"), data = SDIpsilo.ind)
e.XP111 <- lmbreak(score ~ 0 + bp(time, "111"), data = SDIpsilo.ind)
e.XP11 <- lmbreak(score ~ 0 + bp(time, "11"), data = SDIpsilo.ind)
e.XP1010 <- lmbreak(score ~ 0 + bp(time, "1010"), data = SDIpsilo.ind)

# plot results from the models
plot(e.XP101, ylim = c(0, 12))
plot(e.XP111, ylim = c(0, 12))
plot(e.XP11, ylim = c(0, 12))
plot(e.XP1010, ylim = c(0, 12))

# table summary of model
model.tables(e.XP101)
model.tables(e.XP111)
model.tables(e.XP11)
model.tables(e.XP1010)
rm(e.XP101, e.XP11, e.XP111, e.XP1010)

# multiple patients: 101 VS 1010 visual comparison
e.XPall101 <- mlmbreak(score ~ 0 + bp(time, "101"),
  cluster = "id",
  data = SDIpsilo, trace = FALSE
)
e.XPall1010 <- mlmbreak(score ~ 0 + bp(time, c("1010", "101", "11")),
  cluster = "id",
  data = SDIpsilo, trace = FALSE
)
plot(e.XPall101, linewidth = 0, breakpoint = F)
plot(e.XPall1010) # doing 1010 then 101 if no convergence seems appropriate
plot(e.XPall1010, ylim = c(0, 10), cluster = 2:4) # focus on some clusters

summary(e.XPall101)
summary(e.XPall1010)

# Breakpoints:
# id pattern   cv continuity        R2                     breakpoint      maxVs
# 1    1010 TRUE       TRUE 0.9892315 84.50704, 175.78947, 274.03509    < 1e-07
# 2    1010 TRUE       TRUE 0.9916394  55.55743, 92.50575, 232.49985 0.00025244
# 3     101 TRUE       TRUE 0.9915031            65.14286, 166.48148    < 1e-07
# 4    1010 TRUE       TRUE 0.9897900   100.4808, 230.0000, 311.6667    < 1e-07
# ...


# exploration of break points distribution
coeffs <- coef(e.XPall1010) %>% mutate(id = as.numeric(id)) %>%
  group_by(id) %>% 
  mutate(name = paste0("bp.", row_number()))
breaks <- coeffs %>% pivot_wider(id_cols = id, values_from = breakpoint)


