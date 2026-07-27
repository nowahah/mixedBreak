# ghit::install_github("bozenne/lmbreak")
library(lmbreak)
data("SDIpsilo", package = "lmbreak")
SDIpsilo$id <- as.factor(as.numeric(SDIpsilo$id))
data10 <- SDIpsilo[SDIpsilo$time < 190,]

# Model estimation ====
library(segmented)
# library(lme4)

psi0 <- mean(range(data10$time))
data10$U0 <- pmax(0, data10$time - psi0)
pre.fit <- lme(
  score ~ 0 + time + U0 , 
  data = data10,
  random = ~time|id,
  na.action = na.omit
)
segmented.fit <- segmented.lme(
  pre.fit,  ~ time,
  random = list(id = pdDiag(~ time + U + G0 - 1)),
  psi.link = "logit"
)
plot(segmented.fit, vline = TRUE)

range(segmented.fit$psi.i) ## 54.92810 95.67344
plot(segmented.fit, vline = TRUE, xlim = c(0,75)) ## All breakpoints around between 60-75

# Here the displayed individuals breakpoint values are close to identical
# whereas in the actual estimation psi.i, they are much more spread out.

# The problem appears to come from lines 65-67 in plot.segmented.lme function :
# when extracting individuals psi from the fitted object with obj$misc$matrix.psi[, 1]
# # HERE - psi used are != obj$psi.i
# # psi <- if (ncol(obj$misc$matrix.psi) <= 1) 
# #   obj$misc$matrix.psi[, 1]
# # else rowSums(obj$misc$matrix.psi[, 1:(level + 1)])
# psi <- obj$psi.i # correction

# But then what is obj$misc$matrix.psi ?