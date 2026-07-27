### confint.mixedBreak.R ---
## * confint.mixedBreak (documentation)
##' @title Confidence Intervals in segmented Mixed-Models
##' @description Computes confidence intervals for all regression parameters, 
##' including the breakpoint, in a fitted ‘segmented mixed’ model.

## Method for object of class 'mixedBreak'
##' @export
confint.mixedBreak <- function(z, level = 0.95){
  fixed <- lme4::fixef(z$lme.fit)
  fixed.se <- z$fixed[,"Std. Error"]
  
  alpha.lvl <- 1 - level
  ci <- matrix(
    fixed, nrow = 2, ncol = length(fixed), byrow = TRUE,
    dimnames = list(paste0(100*c(alpha.lvl/2, 1-alpha.lvl/2), "%"),
                    names(fixed.se))
  )
  ci <- ci + t(t(c(-1, 1)*qnorm(1 - alpha.lvl/2))) %*% fixed.se

  # rescale breakpoint parameters
  # TODO - only for a shared scale for both breakpoints
  n.psi <- nchar(z$pattern) - 1L
  a <- z$psi.range
  
  expit <- function(x) (a[1] + a[2]*  exp(x))/(1+exp(x))
  ci[,stringr::str_detect(colnames(ci), "break.")] <- 
    expit(ci[,stringr::str_detect(colnames(ci), "break.")])
  colnames(ci)[stringr::str_detect(colnames(ci), "break.")] <- 
    paste0("break.", 1:n.psi)
  
  ci
}