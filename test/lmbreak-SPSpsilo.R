library(lmbreak) # Brice's package for patient-wise linear regression with breakpoints
library(dplyr) # data manipulation
library(tidyr)
library(ggplot2) # data visualization

source("R/DataLoader.R")
rm(SDIpsilo.ext)

SPSpsilo <- basel %>% filter(Study=="SPS")
summary(SPSpsilo)

## multiple patients: 101 VS 1010 visual comparison
# UNSTABLE RESULTS
set.seed(1217170426) # unstable despite seed setting
e.XPall101 <- mlmbreak(
  score ~ 0 + bp(time.num, "101"), 
  cluster = "ID",
  data = SPSpsilo, trace = FALSE
)
e.XPall1010 <- mlmbreak(
  score ~ 0 + bp(time.num, c("1010", "101", "11")), 
  cluster = "ID",
  data = SPSpsilo, trace = FALSE
)
plot(e.XPall101)
# summary(e.XPall101)
summary(coef(e.XPall101, type="R2")$R2)
plot(e.XPall1010)
# summary(e.XPall1010)
summary(coef(e.XPall1010, type="R2")$R2)

# better fit by considering trailing as a potential plateau
# main plateau identification feels visually better with this rescue strategy


## exploration of break points distribution
coeffs <- coef(e.XPall101) %>% 
  group_by(ID) %>% 
  mutate(name = paste0("bp.", row_number()))
breaks <- coeffs %>% pivot_wider(id_cols = ID, values_from = breakpoint)
mean(breaks$bp.1); sd(breaks$bp.1); range(breaks$bp.1)
# average (sd) [range] time to plateau:  
# - pattern 101 : 81.3 (43.9) min [15.7 ; 161.8]
# - pattern 1010: 89.5 (57.6) min [15.1 ; 237.6]
# RESULTS ARE UNSTABLE (up to tenths of a minute)

mean(breaks$bp.2 - breaks$bp.1); sd(breaks$bp.2 - breaks$bp.1)
range(breaks$bp.2 - breaks$bp.1)
# average (sd) plateau duration: 
# - pattern 101 : 48.0 (36.8) min [6.2  ; 131.9]
# - pattern 1010: 75.9 (41.4) min [18.0 ; 170.4] BIG DIFFERENCES


## slopes distributions
coeffs <- coef(e.XPall101, type=c("slope")) %>%
  group_by(ID) %>% 
  mutate(name = paste0("beta.", row_number())) %>%
  filter(slope!=0)
slopes <- coeffs %>% pivot_wider(id_cols = ID, values_from = slope)
round(10*c(mean(slopes$beta.1), sd(slopes$beta.1), range(slopes$beta.1)), 2)
# !!! SPS scale range is [0 - 100]
# mean (sd) [range] onset coefficient for 10 min: 
# - pattern 101 : 16.88 (15.00) SPS/10min [3.37 ; 54.00]
# - pattern 1010: 16.84 (15.04) SPS/10min [3.02 ; 54.00] VERY SIMILAR

round(10*c(mean(slopes$beta.3), sd(slopes$beta.3), range(slopes$beta.3)), 2)
# mean (sd) [range] offest coeff for 10 min:
# - pattern 101 : -3.23 (1.01) SPS/10min [-4.88  ; -0.59]
# - pattern 1010: -9.33 (7.14) SPS/10min [-29.74 ; -0.59] BIG DIFFERENCES
