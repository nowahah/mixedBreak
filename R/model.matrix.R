
## * model.matrix.lmbreak.lme
##' @export
model.matrix.lmbreak.lme <- function(formula, data){
  expl.vars <- attr(terms.formula(formula), "term.label")
  data <- data[[expl.vars]] # only keep explanatory variable in the data
  data <- as.data.frame(unclass(XX), stringsAsFactors = TRUE)
  XX <- data[[is.numeric(data)]]
  
  browser()
}