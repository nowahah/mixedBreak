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
