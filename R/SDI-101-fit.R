
# on purpose very simplified version of a fit taking as only formula
# score ~ 0 + time
# individual fit of 101 breakpoint model based on Muggeo's approach

SDI101.fit <- function(psi0 = rep(NA_real_, 2), data = NA, tol.psi = 1e-3, 
                       it.max = 100L, model = T, xx = F, yy = F){
  
  cl <- match.call() # call not so informative due to function's signature
  x <- data$time[!is.na(data$score)]
  y <- data$score[!is.na(data$score)]
  n.points <- length(x)
  n.psi <- length(psi0)
  
  if (any(is.na(psi0)) ){
    # default uniform breakpoint initialization
    psi0[is.na(psi0)] <- seq(min(x), max(x), length.out = 2 + n.psi)[2:(1+n.psi)]
  }
  psi0 <- sort(psi0)
  
  # main loop for the algorithm
  it <- 0
  psi.diff <- tol.psi + 1
  psi.history <- matrix(NA_real_, nrow = it.max+1, ncol = n.psi)
  psi.history[1,] <- psi0
  # L-Inf norm is used for convergence of breakpoints
  while(psi.diff > tol.psi & it < it.max){
    it = it + 1
    
    # 1. Fix Psi and calculate:
    U1 <- x - pmax(0, x - psi.history[it, 1])
    U2 <- pmax(0, x - psi.history[it, 2])
    V1 <- 1*(x > psi.history[it, 1])
    V2 <- - (x > psi.history[it, 2])
    
    # 2. Fit the model with the additional covariates U, V:
    mod.lin <- lm(y ~ 0 + U1 + U2 + V1 + V2, x = T)
    mod.coef <- as.list(mod.lin$coefficients)
    
    # 3. Improve breakpoints estimates:
    psi.history[it+1, 1] <- psi.history[it, 1] + mod.coef$V1/mod.coef$U1
    psi.history[it+1, 2] <- psi.history[it, 2] + mod.coef$V2/mod.coef$U2
    psi.diff <- max(abs(psi.history[it+1,] - psi.history[it,])) # L-Inf norm
  }
  
  # WARNING if non convergence (it.max reached)
  if(it==it.max) {
    warning(paste0("Algorithm did not converge: max number of iterations 'it.max=",
                   it.max,"' reached"))
  }
  
  XX <- mod.lin$x
  # varmat <- solve(t(XX) %*% XX) * sum((mod.lin$residuals)^2)/(n.points-(1+n.psi))
  z <- list(
    coefficients = c(seg1.time = mod.coef$U1, seg3.time = mod.coef$U2),
    coefficients.sd = sqrt(diag(solve(t(XX) %*% XX))[1:2]),
    psi = psi.history[it+1,],
    peak.y = mod.coef$U1*psi.history[it+1, 1],
    # psi.sd = NA_real_,
    psi.diff = psi.diff,
    psi.history = psi.history[1:(it+1),],
    n.iter = it,
    call = cl
  )
  
  # TODO - return model frame (from last lm) and data (optional)
  if(model) z$XX <- XX
  if(xx) z$x <- x
  if(yy) z$y <- y
  
  
  class(z) <- append("break.lm", class(z))
  return(z)
}

library(lmbreak)
library(ggplot2)
library(dplyr)
data(SDIpsilo, package = "lmbreak")


# algorithm is sensible to starting points, let's give the real values for initialization first
library(lmbreak)
e.XPall <- mlmbreak(score ~ 0 + bp(time, "101"), cluster = "id", data = SDIpsilo,
                    trace = FALSE, digits = 2)
plot(e.XPall)
breakpoints <- model.tables(e.XPall)[rep((1:15-1)*4,each=3)+2:4, 2]
breakpoints <- matrix(breakpoints, ncol = 3, byrow = T) # breakpoints
breakpoints


# start with right initialization
fit101.res <- list()
for(ID in levels(SDIpsilo$id)){
  SDI.ind <- SDIpsilo %>% filter(id==ID)
  ggplot(SDI.ind, aes(x=time, y=score)) + 
    geom_point()
  tryCatch({
    break101.fit <- SDI101.fit(data = SDI.ind, 
                               psi0 = breakpoints[as.numeric(ID), 1:2],
                               model = FALSE)
    fit101.res[[as.numeric(ID)]] <- break101.fit
  }, error = function(e){
    # browser()
    fit101.res[[as.numeric(ID)]] <<- NA

    print(paste("ID ", ID, "- ERROR:"))
    print(e$message)
  }
  # warning = function(w){
  #   fit101.res[[as.numeric(ID)]] <<- NA
  #   
  #   print(paste("ID ", ID, "- WARNING:"))
  #   # browser()
  #   print(w$message)
  # }, 
  # finally = {
  #   names(fit101.res)[length(fit101.res)] <- as.numeric(ID)
  # }
  )
}

# warnings for individuals 11 and 15 (again)
head(fit101.res[[11]]$psi.history)
head(fit101.res[[15]]$psi.history)
# both due to 2-periodic oscillations after init, in range of time domain and in order
# Figure out CV // Continuity in Brice's package
summary(e.XPall)


# to make it better:
# - bootstrap restart ? study what and how - bootstrap on the initial breakpoint? 
#   -> explained in the mixed model paper
# - look for a finer initialization than uniform (splines nadir, other...)
# - in case of oscillation: one solution might be not practically relevant
#   (negative or unordered breakpoints) while other could be ? heuristic

# Add plot, summary, print methods (3)
# would help simplify the analysis of results too