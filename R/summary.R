## Method for object of class 'mixedBreak'
##' @export
summary.mixedBreak <- function(z, default = FALSE){
  if(default) {
    return(summary.default(z))
  }
  require(rlang)
  pattern <- z$pattern

  # browser()
  x <- z$lme.fit
  LL <- as.numeric(logLik(x))
  nameZ <- z$var.name$segmented
  resd <- resid(x, type = "pearson") # Muggeo uses lme.fit.noG assuming known break 
  if (length(resd) > 5) {
    resd <- quantile(resd, na.rm = TRUE)
    names(resd) <- c("Min", "Q1", "Med", "Q3", "Max")
  }
  
  cat("Segmented mixed-effects model fit by REML \n")
  print(data.frame(LogLik = LL, REML = lme4::REMLcrit(x), row.names = " "))
  
  cat("\nRandom effects:\n")
  tmp <- lme4::VarCorr(x) # TODO - arrange names more nicely
  print(tmp)
  # FINAL - adapt names of variables for printing like below
  
  cat("\nFixed effects:\n")
  fixed <- z$fixed
  fixed[stringr::str_detect(rownames(fixed), "break"), "t value"] <- NA
  print(fixed)
  cat(" psi.link = logit\n")
  
  # browser()
  n.obs <- x@devcomp$dims[["N"]]
  n.group <- nlevels(x@frame[,names(x@cnms)])
  n.psi <- nchar(z$pattern)-1L
  t.val <- summary(z$lme.fit.check)$coefficients[, "t value"]
  t.val <- t.val[names(t.val) %in% paste0("VD", c("", 1:n.psi))]
  p.val <- pt(t.val, df = n.obs - (n.group + length(x@beta)) + 1) # significance of eta
  # browser()
  cat("\nSteady breakpoint linear predictor at convergence:\n") 
  cat(paste(paste0(" - eta.break.", 1:n.psi), round(p.val, 3), sep = ": p = ", 
            collapse = " \n"), "\n")
  cat("Low p-value indicates evidence of convergence issues.\n")
  #PSILINK
  
  cat("\nIndividual breakpoint(s) summary:\n")
  break.summ <- summary(z$psi.i[,1])
  if(n.psi > 1){
    for(i in 2:(n.psi))
    break.summ <- rbind(break.summ, summary(z$psi.i[,i]))
  }
  if(n.psi > 1) rownames(break.summ) <- paste0("break.", 1:n.psi)
  print(break.summ)
  
  # TODO - check validity since Muggeo excluded breakpoint estimate from the fit here
  # cat("\nStandardized Within-Group Residuals:\n")
  # print(resd)
  
  cat("\nPattern:", pattern)
  cat("\nApproximation:", z$approx)
  
  cat("\nNumber of Observations:", n.obs)
  cat("\nNumber of Groups:", n.group, "\n")
  
  invisible(z)
}
