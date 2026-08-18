#' Extract mixed.break Coefficients
#' 
#' Extract summary statistics and regression coefficients from the mixed-effects
#' segmented model.
#' 
#' @usage 
#' \usage{
#' \method{coef}{mixedBreak}(break.fit, type = "cluster")
#' }
#' 
#' @param z an object of class \code{mixedBreak}
#' @param type [vector of character] CURRENTLY UNUSED. a character vector indicating 
#' which type of coefficients to extract and return. Several options include:
#' \itemize{
#' \item auc: area under the curve\n
#' \item breakpoint: position of the breakpoint(s)\n
#' \item breakpoint.range: Minimum and Maximum theoretical value of the
#'  breakpoint(s)\n
#' \item cluster: all subjects-specific regression coefficients (slopes + 
#' breakpoints) \n
#' \item corr: \n
#' \item duration: duration of each phases (same unit as segmented variables)\n
#' \item slope: slopes associated with each phases \n
#' \item R2: Nakagawa's coefficient of determination for mixed-model fits 
#' (fraction of variance explained by fixed-effects only)\n
#' \item R2_cond: Conditional coefficient of determination for mixed-model fits
#' (fraction of variance explained by fixed+random-effects)\n
#' \item pattern: pattern used
#' }
#'
#' @return a data.frame with requested coefficients. Eventually a list of 2 
#' for subjects-level & population-level statistics.
#' @export
#'
#' @examples
#' data(SDIpsilo)
#' break.fit <- mixed.break(
#'   score ~ 0 + time + (0 + time|id),
#'   pattern = "10",
#'   data = SDIpsilo[SDIpsilo$time < 190, ]
#' )
#' coef(break.fit)
coef.mixedBreak <- function(z, type = "cluster"){
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