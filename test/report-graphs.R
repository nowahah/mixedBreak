library(lmbreak)
library(dplyr)
library(ggplot2)

data(SDIpsilo, package = "lmbreak")
SDIpsilo <- SDIpsilo %>% filter(type %in% c("noise", "trailing") == FALSE) # article's analysis



ind13.fit <- lmbreak(score ~ 0 + bp(time, "101"), data = SDIpsilo %>% filter(id==13),
                     trace = FALSE)
summary(ind13.fit)
p <- plot(ind13.fit, title = "")
p$plot +
  scale_y_continuous(breaks = seq(0,10,by=2), limits = c(0, 10)) + # custom y ticks
  ylab("SDI") +
  xlab("Time since drug intake (minutes)")
ggsave("inst/figs/SDI-13-fit.png", width = 148, height = 105, units = "mm")
