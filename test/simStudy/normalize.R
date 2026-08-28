normalize <- function(object){
  UseMethod("normalize")
}

## method for object of class 'segmented'
# TODO - to adapt to call ?apply function
normalize.segmented <- function(x){
  res <- list(
    breaks.ind = x$psi[,2],
    breaks.sd = x$psi[,3],
    breaks.CI95.low = confint(x)[,2],
    breaks.CI95.up = confint(x)[,3],
    intercept = intercept(x)$time[2]
  )
  
  return(res)
}

## method for object of class 'tbl_df' (containing segmented models)
normalize.tbl_df <- function(x){
  tmp <- apply(x, 1, function(z) normalize(z$model))
  res <- matrix(unlist(tmp), ncol = 5, nrow = nrow(x), byrow = T)
  res <- data.frame(res)
  colnames(res) <- paste0("breaks.", c("x1", "x1.sd", "CI95.low", "CI95.up", "y1"))
  
  return(res)
}

## method for object of class 'segmented.lme'
normalize.segmented.lme <- function(x){
  tmp <- summary(x$lme.fit)
  res <- list(
    breaks.avg = x$lme.fit$coefficients$fixed['G0'],
    breaks.ind = as.numeric(x$psi.i),
    breaks.sd = sqrt(tmp$varFix["G0", "G0"]),
    breaks.random.sd = as.numeric(VarCorr(x$lme.fit)['G0', 'StdDev']),
    breaks.CI95.low = confint(x)[1,'G0'],
    breaks.CI95.up = confint(x)[2,'G0'],
    intercept = x$lme.fit$coefficients$random$id[,"(Intercept)"] + 
      x$lme.fit$coefficients$fixed["time"] * as.numeric(x$psi.i),
    slope2 = sum(x$lme.fit$coefficients$fixed[c("time", "U")])
    # break.intercept = intercept.rd + time*break.x.rd
  )
  
  return(res)
}

## method for object of class 'mlmbreak'
normalize.mlmbreak <- function(x){
  res <- list(
    breaks.ind = unlist(lapply(x$breakpoint, function(x) x["value"]), use.names = F),
    breaks.sd = NA,
    intercept = unlist(lapply(x$phase, function(x) x[2, "intercept"]))
  )
  
  return(res)
}


## method for object of class 'mixedBreak'
normalize.mixedBreak <- function(x){
  # browser()
  summ <- summary(x$lme.fit)
  n.coef <- length(x$lme.fit@beta)
  a <- x$psi.range
  expit <- function(bp) { return((a[1]+a[2]*exp(bp)) / (1+exp(bp))) }
  d.expit <- function(bp) { return( exp(bp)*diff(a) / (1+exp(bp))^2 ) }
  break.x1.sd <- summ$coefficients[n.coef, "Std. Error"] * d.expit(summ$coefficients[n.coef, "Estimate"])
  M <- 1e4
  eta.rd.sample <- rnorm(M, x$lme.fit@beta[length(x$lme.fit@beta)], attr(lme4::VarCorr(x$lme.fit)$ID, "stddev")['G1'])
  break.x1.random.sd <- sd(expit(eta.rd.sample))
  res <- list(
    break.x1.avg = expit(x$lme.fit@beta[length(x$lme.fit@beta)]),
    break.x1.ind = x$psi.i,
    break.x1.sd = break.x1.sd,
    break.x1.random.sd = break.x1.random.sd,
    break.CI95.low = confint(x)[1, n.coef],
    break.CI95.up = confint(x)[2, n.coef],
    break.y1.ind = x$psi.i * x$random[[1]], 
    slope2 = ifelse(x$pattern == "11", sum(x$fixed[1:2, 1]), 0),
    slope2.sd = ifelse(x$pattern == "11", sum(x$fixed[1:2, 2])+ 2*summ$varcor$ID[1,2], NA)
    # break.intercept = intercept.rd + time*break.x.rd
  )
  
  return(res)
}