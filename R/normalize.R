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
      x$lme.fit$coefficients$fixed["time"] * as.numeric(x$psi.i)
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