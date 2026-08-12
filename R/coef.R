coef.mixedBreak <- function(z, type = "all"){
  vars <- z$var.name
  ID <- levels(z$model[[vars$group]])
  breakpoints <- z$psi.i
  n.psi <- ncol(breakpoints)
  # browser()
  slopes <- z$random[,colnames(z$random) %in% paste0("U", c("", 1:(n.psi+1)))]
  slopes <- data.frame(slopes)
  names(slopes) <- paste(vars$segmented, which(strsplit(z$pattern, "")[[1]] == "1"),
                         sep = ".")
  # TODO - beta/delta naming in function of pattern
  
  res <- data.frame(
    ID, breakpoints, slopes
  )
  names(res) <- c(
    vars$group, paste0("break.", 1:n.psi), 
    paste(vars$segmented, which(strsplit(z$pattern, "")[[1]] == "1"), sep = ".")
  )
  rownames(res) <- 1:nrow(res)
  res
}