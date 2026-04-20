## Reproduce the analysis made from Stenbæck et al. 2020 article
# From lmbreak github page: https://github.com/bozenne/lmbreak/tree/main

library(lmbreak)


data(SDIpsilo, package = "lmbreak")
SDIpsilo <- SDIpsilo[SDIpsilo$type %in% c("noise","trailing") == FALSE,] # labelled by hand
str(SDIpsilo)


## Multiple clusters (patients)
e.XPall <- mlmbreak(score ~ 0 + bp(time, "101"), cluster = "id", data = SDIpsilo,
                    trace = FALSE)
summary(e.XPall)
tbl <- model.tables(e.XPall)
# plot(e.XPall)

breakpoints <- tbl[rep((1:15-1)*4,each=2)+2:4, 2]
breakpoints <- matrix(breakpoints, ncol = 3, byrow = T) # breakpoints
breakpoints <- as.data.frame(breakpoints)
names(breakpoints) = paste0("break.x", 1:3)
intercept <- coef(e.XPall, type="intercept")[rep((1:15-1)*4)+2, 2]
b.onset <- tbl[(1:15-1)*4+1, 5] # onset slopes
b.return <- tbl[(1:15-1)*4+3, 5] # return to normal slopes
peak.duration <- tbl[(1:15-1)*4+3, 2] # duration of plateau

library(rgl)
plot3d(x=b.onset, y=b.return, z=peak.duration,
       type="s", radius = 3,
       rgl.printRglwidget = TRUE)
rglwidget() # click on 'Show in new window' in the plot pane (next to the broom)

## 0 - breakpoints and intercept
summary(breakpoints$break.x1) ; sd(breakpoints$break.x1)
# mean (sd) [range] = 77.2 (38.2) [32.5 ; 157.9]
summary(breakpoints$break.x2 - breakpoints$break.x1)
sd(breakpoints$break.2 - breakpoints$break.1)
# mean (sd) [range] = 77.5 (33.3) [29.0 ; 142.8]
summary(breakpoints$break.x3 - breakpoints$break.x2)
sd(breakpoints$break.x3 - breakpoints$break.x2)
# mean (sd) [range] = 131.6 (39.3) [90.0 ; 204.1]
# NOTE: break.3 is the time of last observation during experiment,
# not the time of return to SDI=0 like reported in the article.

# intercept at plateau
summary(intercept); sd(intercept)
# mean (sd) [range] = 9.3 (1.0) [7.5 ; 10.0]

## 1 - onset VS return
plot(b.onset, b.return)
cor.test(b.onset, b.return) # non significant correlation
## 2 - onset VS plateau duration
plot(b.onset, peak.duration)
cor.test(b.onset, peak.duration) # significant correlation
## 3 - return VS plateau duration
plot(peak.duration, b.return)
cor.test(b.return, peak.duration) # significant correlation

## correlogram reproduced
library(corrplot)
corrplot(cor(data.frame(b.onset, peak.duration, b.return)), 
         type = "upper", diag = F, method = "ellipse", 
         col = rev(COL2('RdBu', 200)))
