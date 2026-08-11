coef.mixedBreak <- function(z, type = "all"){
  vars <- z$var.name
  ID <- unique(z$model[[vars$group]])
  breakpoints <- z$psi.i
  slopes <- z$random[,colnames(z$random) %in% paste0("U", 1:(ncol(z$psi.i)+1))]
  names(slopes) <- paste(vars$segmented, which(strsplit(z$pattern, "")[[1]] == "1"),
                         sep = ".")
  # TODO - beta/delta naming in function of pattern
  
  res <- data.frame(
    ID, breakpoints, slopes
  )
  n.psi <- ncol(z$psi.i)
  names(res) <- c(
    vars$group, paste0("break.", 1:n.psi), 
    paste(vars$segmented, which(strsplit(z$pattern, "")[[1]] == "1"), sep = ".")
  )
  rownames(res) <- 1:nrow(res)
  res
}