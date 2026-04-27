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

## method for object of class 'segmented.lme'
normalize.segmented.lme <- function(x){
  res <- list(
    breaks.avg = x$lme.fit$coefficients$fixed['G0'],
    breaks.ind = as.numeric(x$psi.i),
    breaks.random.sd = as.numeric(VarCorr(x$lme.fit)['G0', 'StdDev']),
    breaks.CI95.low = confint(x)[1,'G0'],
    breaks.CI95.up = confint(x)[2,'G0'],
    intercept = x$lme.fit$coefficients$random$id[,"(Intercept)"] + 
      x$lme.fit$coefficients$fixed["time"] * x$psi.i
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