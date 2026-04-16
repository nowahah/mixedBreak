library(lmbreak)
library(segmented)
library(corrplot)

data("SDIpsilo")
SDIpsilo$type <- as.factor(SDIpsilo$type)

## Comparison in estimates from packages lmbreak and segmented on:
#   - breakpoints coordinates
#   - slopes estimates
n.breaks = 2
breaks.lmb <- matrix(0, nrow = nlevels(SDIpsilo$id), ncol = 2*n.breaks,
                     dimnames = list(rep(1:nlevels(SDIpsilo$id)), 
                                     paste0(rep(c("x", "y"), n.breaks), 
                                            rep(1:n.breaks, each = 2))))
breaks.seg <- matrix(0, nrow = nlevels(SDIpsilo$id), ncol = 2*n.breaks,
                     dimnames = list(rep(1:nlevels(SDIpsilo$id)), 
                                     paste0(rep(c("x", "y"), n.breaks), 
                                            rep(1:n.breaks, each = 2))))
slopes.lmb <- matrix(0, nrow = nlevels(SDIpsilo$id), ncol = 2,
                     dimnames = list(1:nlevels(SDIpsilo$id), c("onset", "return")))
slopes.seg <- matrix(0, nrow = nlevels(SDIpsilo$id), ncol = 2,
                     dimnames = list(1:nlevels(SDIpsilo$id), c("onset", "return")))


# patient-wise estimations
for(id in levels(SDIpsilo$id)){
  sdi.ind <- SDIpsilo[SDIpsilo$id==id,] # patient subset
  
  ## lmbreak
  mod.lmb <- lmbreak(score ~ 0 + bp(time, "101"), data = sdi.ind)
  breaks.lmb[id, 1+2*0:(n.breaks-1)] <- mod.lmb$breakpoint$value # time coord
  breaks.lmb[id, 2*1:n.breaks] <- mod.lmb$phase$intercept[1+1:n.breaks] # intercept
  slopes.lmb[id,] <- mod.lmb$phase$slope[c(1,3)] # slopes
  
  ## segmented
  mod.seg <- segreg(score ~  0 + seg(time, npsi=n.breaks, est=c(1,0,1)), data=sdi.ind)
  breaks.seg[id, 1+2*0:(n.breaks-1)] <- mod.seg$psi[,2] # time coord
  slopes.seg[id,] <- mod.seg$coefficients[1:n.breaks] # slopes
  breaks.seg[id, 2*1:n.breaks] <- slopes.seg[id,1]*breaks.seg[id,1] # intercept
}

## Comparison
# onset slopes
sum(sqrt((slopes.lmb[,1]-slopes.seg[,1])^2)) # 0.03878301
# corrplot(cor(slopes.lmb), method = "ellipse", type = "upper", diag = F) 
cor(slopes.lmb) # 0.3297237

# offset slopes
sum(sqrt((slopes.lmb[,2]-slopes.seg[,2])^2)) # 0.01180891
cor(slopes.seg) # 0.2474174 - more moderate

# break 1 coordinates
plot(breaks.lmb[,1] - breaks.seg[,1], pch = 18, xlab = "Patient id",
     main = "Difference in breakpoint 1 time estimates", ylab = "relative difference") 
abline(h=0, lty="dashed")
# important differences on ind 2, 7, 8, 15?

plot(breaks.lmb[,2] - breaks.seg[,2], pch = 18, xlab = "Patient id",
     main = "Difference in breakpoint 1 height estimates", ylab = "relative difference") 
abline(h=0, lty="dashed")
# idem except on 15

cor(breaks.lmb)


# break 2 coordinates
plot(breaks.lmb[,3] - breaks.seg[,3], pch = 18, xlab = "Patient id",
     main = "Difference in breakpoint 2 time estimates", ylab = "relative difference") 
abline(h=0, lty="dashed")
# important difference on ind 2, 7, 8

plot(breaks.lmb[,4] - breaks.seg[,4], pch = 18, xlab = "Patient id",
     main = "Difference in breakpoint 2 height estimates", ylab = "relative difference") 
abline(h=0, lty="dashed")
# ind 2, 7, 8, 11, 15


## graphical comparison on ind 2, 7, 8
cluster = c(2, 7, 8)
sdi.cluster <- SDIpsilo[SDIpsilo$id %in% cluster,]
mod.lmb <- mlmbreak(score ~ 0 + bp(time, "101"), cluster = "id", data = sdi.cluster)
p <- plot(mod.lmb)$plot
p + 
  geom_point(
    aes(x=x, y=y), 
    data=data.frame(
      # id label is not exact if cv != TRUE, TRUE - ok for 2, 7, 8
      id=as.factor(paste(rep(cluster,2), "(cv = TRUE,TRUE)")),
      x=c(breaks.seg[cluster,c(1,3)]),
      y=c(breaks.seg[cluster,c(2,4)])
    ), shape=18, size=3, color="blue"
  )

