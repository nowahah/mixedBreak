### confint.mixedBreak.R ---
## * confint.mixedBreak (documentation)
##' @title Confidence Intervals in segmented Mixed-Models
##' @description Computes confidence intervals for all regression parameters, 
##' including the breakpoint, in a fitted ‘segmented mixed’ model.

## Method for object of class 'mixedBreak'
##' @export
confint.mixedBreak <- function(object, level = 0.95){
  fixed <- lme4::fixef(object$lme.fit)
  fixed.se <- object$fixed[,"Std. Error"]
  
  alpha.lvl <- 1 - level
  ci <- matrix(
    fixed, nrow = 2, ncol = length(fixed), byrow = TRUE,
    dimnames = list(paste0(100*c(alpha.lvl/2, 1-alpha.lvl/2), "%"),
                    names(fixed.se))
  )
  ci <- ci + t(t(c(-1, 1)*qnorm(1 - alpha.lvl/2))) %*% fixed.se
  
  ci
}