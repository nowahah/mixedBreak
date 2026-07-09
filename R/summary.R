## Method for object of class 'break.lm'
##' @export
summary.break.lm <- function(z, default = FALSE){
  if(default) {
    summary.default(z)
  } else {
    cat("\nCall:\n", paste(deparse(z$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    
    # Residuals:
    # View(summary.lm)
    # ans <- z[c("...")]
    # ...
    
    # Breakpoints:
    # Estimate, Std.Error, t value, Pr(>|t|) . * ** ***
    
    # Coefficients:
    
    # Residual std err
    # R^2
    # F-statistic
  }
}


## Method for object of class 'mixedBreak1'
##' @export
summary.mixedBreak1 <- function(z, pattern = "11", default = FALSE){
  if(default) {
    return(summary.default(z))
  }
  pattern <- z$pattern
  stopifnot(pattern %in% c("11", "10"))
  
  x <- z$lme.fit
  LL <- as.numeric(logLik(x))
  nameZ <- z$var.name$segmented
  resd <- resid(x, type = "pearson") # Muggeo uses lme.fit.noG assuming known break 
  if (length(resd) > 5) {
    resd <- quantile(resd, na.rm = TRUE)
    names(resd) <- c("Min", "Q1", "Med", "Q3", "Max")
  }
  
  cat("Segmented mixed-effects model fit by REML \n")
  print(data.frame(AIC = AIC(x), BIC = BIC(x), logLik = LL, 
                   row.names = " "))
  
  cat("\n")
  cat("Random effects:\n")
  tmp <- summary(x)$varcor[[z$var.name$group]]
  dim(attr(tmp, "stddev")) <- c(1, length(attr(tmp, "stddev")))
  if(pattern=="11"){
    colnames(attr(tmp, "stddev")) <- c("time.beta.1", "time.delta.2", "break.1")
  } else {
    colnames(attr(tmp, "stddev")) <- c("time.beta.1", "break.1")
  }
  rownames(attr(tmp, "stddev")) <- "StdDev:"
  print(attr(tmp, "stddev"))
  # FINAL - adapt names like below
  
  cat("\nFixed effects:\n")
  # differentiate the segment with good variable names only ?
  # or like M add ----left/rightS / breakpoint
  mod.coef <- summary(x)$coefficients
  row.names(mod.coef)[row.names(mod.coef)==nameZ] <- paste0(nameZ, ".beta.1")
  # all variables starting with U, G
  var.U <- substr(row.names(mod.coef), 1, 1) == "U"
  if(pattern=="11"){
    row.names(mod.coef)[var.U] <- paste0(nameZ, ".delta.", 2:(1+length(which(var.U))))
  } else {
    row.names(mod.coef)[var.U] <- nameZ
  }
  var.G <- substr(row.names(mod.coef), 1, 1) == "G"
  row.names(mod.coef)[var.G] <- paste0("(logit.)break.", 1:length(which(var.G)))
  mod.coef[var.G, 3] <- NA
  print(mod.coef) # TODO - psi is still on logit scale
  cat(" psi.link = logit\n")
  
  # TODO - check validity since Muggeo excluded breakpoint estimate from the fit here
  # cat("\nStandardized Within-Group Residuals:\n")
  # print(resd)
  cat("\nNumber of Observations:", x@devcomp$dims[["N"]])
  cat("\nNumber of Groups:", nlevels(x@frame[,names(x@cnms)]), "\n")
  
  invisible(z)
}