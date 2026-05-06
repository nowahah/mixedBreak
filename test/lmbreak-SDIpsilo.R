# library(remotes)
# install_github("bozenne/lmbreak")

library(lmbreak) # Brice's package for patient-wise linear regression with breakpoints
library(dplyr) # data manipulation
library(tidyr)
library(ggplot2) # data visualization

data(SDIpsilo, package = "lmbreak")
SDIpsilo <- SDIpsilo %>% filter(type %in% c("noise", "trailing") == FALSE) # rm noise and trailing
str(SDIpsilo)

## multiple pattern: 101 VS 1010 visual comparison ====
e.XPall101 <- mlmbreak(score ~ 0 + bp(time, "101"),
  cluster = "id",
  data = SDIpsilo, trace = FALSE
)
e.XPall111 <- mlmbreak(score ~ 0 + bp(time, "111"),
                       cluster = "id",
                       data = SDIpsilo, trace = FALSE
)
e.XPall1010 <- mlmbreak(score ~ 0 + bp(time, c("1010", "101", "11")),
  cluster = "id",
  data = SDIpsilo, trace = FALSE
)
cluster = c(12,13,11,4)
plot(e.XPall101, cluster=cluster)
plot(e.XPall111, cluster=cluster)
plot(e.XPall1010) # doing 1010 then 101 if no convergence seems appropriate
# ggsave("../../figures/SDI/gg-101-exs.png", width = 21, height = 14.8, units = "cm")


summary(e.XPall101)
# summary(e.XPall1010)


## exploration of break points distribution ====
coeffs <- coef(e.XPall101) %>% 
  mutate(id = as.numeric(id)) %>%
  group_by(id) %>% 
  mutate(name = paste0("bp.", row_number()))
breaks <- coeffs %>% pivot_wider(id_cols = id, values_from = breakpoint)
mean(breaks$bp.1); sd(breaks$bp.1)
# average (sd) time to plateau:  77.2 (38.2) min - 
# higher than paper:  71.4 (NA) min


mean(breaks$bp.2 - breaks$bp.1); sd(breaks$bp.2 - breaks$bp.1)
# average (sd) plateau duration: 83.0 (33.3) min - 
# ~ same as paper: 82.4 (33.1) min


## slopes distributions ====
coeffs <- coef(e.XPall101, type=c("slope")) %>%
  mutate(id = as.numeric(id)) %>%
  group_by(id) %>% 
  mutate(name = paste0("beta.", row_number())) %>%
  filter(slope!=0)
slopes <- coeffs %>% pivot_wider(id_cols = id, values_from = slope)
round(10*c(mean(slopes$beta.1), sd(slopes$beta.1), range(slopes$beta.1)), 2)
# mean (sd) [range] onset coefficient for 10 min: 1.44 (0.58) [0.63 ; 2.50]
#                  rounded values from the paper: 1.4  (0.5)  [0.6 ; 2.3] != SD, MAX

round(10*c(mean(slopes$beta.3), sd(slopes$beta.3), range(slopes$beta.3)), 2)
# mean (sd) [range] offest coeff for 10 min: -0.67 (0.26) [-1.25 ; -0.26]
#             rounded values from the paper: -0.7  (0.3)  [-1.3 ; -0.3] OK


# ==================================================

## reproducing last figure with breakpoint average instead of prediction average
p <- plot(e.XPall101, scales = "none", ylim = c(0, 10), color=c("cadetblue3", "springgreen4"),
          linewidth = .5, breakpoint=T, size=c(1, 0))

# gather breakpoint avg
mean(breaks$bp.1)
mean(breaks$bp.2)
coeffs <- coef(e.XPall101, type="intercept")
breaks <- as.data.frame(breaks)
breaks$intercept.peak <- coeffs[4*(1:15-1)+2, 2]
breaks$intercept.end <- coeffs[4*(1:15), 2]
summary(breaks)
break.end <- SDIpsilo %>%
  group_by(id) %>%
  summarise(break.end = max(time[!is.na(score)], na.rm = T)) %>%
  summarise(break.end.mean = mean(break.end))

avg.break <- data.frame(
  x = c(0, mean(breaks$bp.1), mean(breaks$bp.2), as.numeric(break.end)),
  y = c(0, rep(mean(breaks$intercept.peak), 2), mean(breaks$intercept.end))
)

p$plot + 
  geom_point(data=avg.break, mapping=aes(x=x, y=y), shape = 18, size = 3) + 
  geom_line(data=avg.break, mapping=aes(x=x, y=y), size = 1)
ggsave("../../figures/SDI/gg-101-break-avg.png", width = 21, height = 14.8, units = "cm")
