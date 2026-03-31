# if there are different classes to handle,
# UseMethod("funcName") is the way to proceed for generic method
normalize <- function(object){
  UseMethod("normalize")
}

normalize.matrix <- function(object){
  
  print("a")
  
}

normalize.list <- function(object){
  
  print("b")
  
}