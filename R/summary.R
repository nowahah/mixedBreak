## Method for object of class 'mixedBreak1'
##' @export
summary.mixedBreak1 <- function(z, default = FALSE){
  if(default) {
    return(summary.default(z))
  }
  require(rlang)
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
  tmp <- lme4::VarCorr(x) #[[z$var.name$group]]
  # # dim(attr(tmp, "stddev")) <- c(1, length(attr(tmp, "stddev")))
  # if(pattern=="11"){
  #   colnames(tmp) <- rownames(tmp) <- c("time.beta.1", "time.delta.2", "break.1")
  #   # colnames(attr(tmp, "stddev")) <- c("time.beta.1", "time.delta.2", "break.1")
  # } else {
  #   colnames(tmp) <- rownames(tmp) <- c("time.beta.1", "break.1")
  #   # colnames(attr(tmp, "stddev")) <- c("time.beta.1", "break.1")
  # }
  # # rownames(attr(tmp, "stddev")) <- "StdDev:"
  print(tmp)
  # print(attr(tmp, "stddev"))
  # FINAL - adapt names of variables for printing like below
  
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
  
  n.obs <- x@devcomp$dims[["N"]]
  n.group <- nlevels(x@frame[,names(x@cnms)])
  t.val <- summary(z$lme.fit.check)$coefficients["VD1", "t value"]
  p.val <- pt(t.val, df = n.obs - (n.group + length(x@beta)) + 1) # signifiance of eta
  if(p.val < 0.05) {
    eta.flag <- "YES"
  } else {
    eta.flag <- "NO"
  }
  cat("\ndelta*(eta-tilde(eta)) is significantly different from 0 at convergence:", eta.flag,
      "( p.value =", round(p.val, 3), ")\n")
  
  # TODO - check validity since Muggeo excluded breakpoint estimate from the fit here
  # cat("\nStandardized Within-Group Residuals:\n")
  # print(resd)
  cat("\nNumber of Observations:", n.obs)
  cat("\nNumber of Groups:", n.group, "\n")
  
  invisible(z)
}


## Method for object of class 'mixedBreak2'
##' @export
summary.mixedBreak2 <- function(z, default = FALSE){
  if(default) {
    return(summary.default(z))
  }
  require(rlang)
  pattern <- z$pattern
  stopifnot(pattern %in% c("111", "101"))
  
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
  
  cat("\n")
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
  n.psi <- ncol(z$psi.i)
  t.val <- summary(z$lme.fit.check)$coefficients[paste0("VD", 1:n.psi), "t value"]
  p.val <- pt(t.val, df = n.obs - (n.group + length(x@beta)) + 1) # significance of eta
  if(any(p.val > 0.05)) {
    eta.flag <- "NO" # good
  } else {
    eta.flag <- "YES" # bad, we want estimated eta to be steady around former value tilde
  }
  cat("\ndelta*(eta-tilde(eta)) != 0 significantly at convergence:", eta.flag, 
      "- p.values:\n")
  cat(paste(paste0(" eta.", 1:n.psi), round(p.val, 3), sep = ": ", collapse = " \n"), "\n")
  
  # TODO - check validity since Muggeo excluded breakpoint estimate from the fit here
  # cat("\nStandardized Within-Group Residuals:\n")
  # print(resd)
  cat("\nNumber of Observations:", n.obs)
  cat("\nNumber of Groups:", n.group, "\n")
  
  invisible(z)
}
