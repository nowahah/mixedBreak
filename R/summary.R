## Method for object of class 'mixedBreak'
##' @export
summary.mixedBreak <- function(z, default = FALSE){
  if(default) {
    return(summary.default(z))
  }
  require(rlang)
  pattern <- z$pattern
  stopifnot(pattern %in% c("11", "10", "111", "101"))
  
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
  print(data.frame(AIC = AIC(x), BIC = BIC(x), logLik = LL, 
                   row.names = " "))
  
  cat("Pattern:", pattern, "\n\n")
  
  cat("Random effects:\n")
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
  if(all(p.val > 0.05)) {
    eta.flag <- "YES" # good
  } else {
    eta.flag <- "NO" # bad, we want estimated eta to be steady around former value tilde
  }
  cat("\nSteady breakpoint linear predictor at convergence ( eta≈tilde(eta) ):", 
      eta.flag, "- p.values:\n")
  cat(paste(paste0(" - eta.break.", 1:n.psi), round(p.val, 3), sep = ": ", 
            collapse = " \n"), "\n")
  #PSILINK
  
  # TODO - check validity since Muggeo excluded breakpoint estimate from the fit here
  # cat("\nStandardized Within-Group Residuals:\n")
  # print(resd)
  cat("\nNumber of Observations:", n.obs)
  cat("\nNumber of Groups:", n.group, "\n")
  
  invisible(z)
}
