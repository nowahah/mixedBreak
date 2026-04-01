## Reproduce the analysis made from Stenbæck et al. 2020 article
# From lmbreak github page: https://github.com/bozenne/lmbreak/tree/main

library(lmbreak)


data(SDIpsilo, package = "lmbreak")
SDIpsilo <- SDIpsilo[SDIpsilo$type %in% c("noise","trailing") == FALSE,]
str(SDIpsilo)


## Single pattern, Single dataset
SDIpsilo13 <- SDIpsilo[SDIpsilo$id==13,]

e.XP111 <- lmbreak(score ~ 0 + bp(time, "111"), data = SDIpsilo13)
e.XP101 <- lmbreak(score ~ 0 + bp(time, "101"), data = SDIpsilo13)
e.XP11 <- lmbreak(score ~ 0 + bp(time, "11"), data = SDIpsilo13)

plot(e.XP111, ylim = c(0,12)) ## left panel
plot(e.XP101, ylim = c(0,12)) ## middle panel
plot(e.XP11, ylim = c(0,12)) ## right panel

model.tables(e.XP101)

coef(e.XP101, type = "auc", interval = c(0,300))

fit.XP101 <- predict(e.XP101, newdata = data.frame(time = seq(0,440,by=1)))
cbind(head(fit.XP101), "",tail(fit.XP101))

e.XP01 <- lmbreak(score ~ 0 + bp(time, "01"), data = SDIpsilo13)

## Multiple patterns
e.XPrescue <- lmbreak(score ~ 0 + bp(time, c("01","11")), data = SDIpsilo13)
coef(e.XPrescue,c("pattern","breakpoint"))
rm(e.XPrescue, e.XP01, e.XP101, e.XP11, e.XP111, fit.XP101, SDIpsilo13)


## Multiple datasets
e.XPall <- mlmbreak(score ~ 0 + bp(time, "101"), cluster = "id", data = SDIpsilo,
                    trace = FALSE)
# That is where I will be working from
summary(e.XPall)
tbl <- model.tables(e.XPall)
# plot(e.XPall, cluster=11)

n.obs <- 15L
breakpoints <- tbl[rep((1:15-1)*4,each=2)+2:3, 2]
breakpoints <- matrix(breakpoints, ncol = 2, byrow = T) # breakpoints
b.onset <- tbl[(1:15-1)*4+1, 5] # onset slopes
b.return <- tbl[(1:15-1)*4+3, 5] # return to normal slopes
peak.duration <- tbl[(1:15-1)*4+3, 2] # duration of plateau

library(rgl)
plot3d(x=b.onset, y=b.return, z=peak.duration,
       type="s", radius = 3,
       rgl.printRglwidget = TRUE)
rglwidget()

# 1 - onset VS return
plot(b.onset, b.return)
cor.test(b.onset, b.return) # non significant correlation
# 2 - onset VS plateau duration
plot(b.onset, peak.duration)
cor.test(b.onset, peak.duration) # significant correlation
# 3 - return VS plateau duration
plot(peak.duration, b.return)
cor.test(b.return, peak.duration) # significant correlation

library(corrplot)
corrplot(cor(data.frame(b.onset, peak.duration, b.return)), 
         type = "upper", diag = F, method = "ellipse", 
         col = rev(COL2('RdBu', 200)))
