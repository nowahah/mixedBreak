
# on purpose very simplified version of a fit taking as only formula
# score ~ 0 + time
# individual fit of 101 breakpoint model based on Muggeo's approach

SDI101.fit <- function(psi0 = rep(NA_real_, 2), data = NA, tol.psi = 1e-3, it.max = 100L){
  
  cl <- match.call()
  x <- data$time[!is.na(data$score)]
  y <- data$score[!is.na(data$score)]
  n.points <- length(x)
  n.psi <- length(psi0)
  
  if (any(is.na(psi0)) ){
    psi0[is.na(psi0)] <- seq(min(x), max(x), length.out = 2 + n.psi)[2:(1+n.psi)]
    psi0 <- sort(psi0)
  }
  
  # main loop for the algorithm
  it <- 0
  psi.diff <- tol.psi + 1
  psi.old <- psi0
  psi.new <- psi0
  # L-Inf norm is used for convergence of breakpoints
  while(psi.diff > tol.psi & it < it.max){
    it = it + 1
    
    # 1. Fix Psi and calculate:
    U1 <- x - pmax(0, x - psi.new[1])
    U2 <- pmax(0, x - psi.new[2])
    V1 <- 1*(x > psi.new[1])
    V2 <- - (x > psi.new[2])
    
    # 2. Fit the model with the additional covariates
    mod.lin <- lm(y ~ 0 + U1 + U2 + V1 + V2, x = T)
    mod.coef <- as.list(mod.lin$coefficients)
    
    # 3. Improve breakpoints estimates
    # browser()
    psi.old <- psi.new
    psi.new[1] <- psi.old[1] + mod.coef$V1/mod.coef$U1
    psi.new[2] <- psi.old[2] + mod.coef$V2/mod.coef$U2
    psi.diff <- max(abs(psi.old - psi.new))
    
    # debug
    # print(psi.new)
    # print(psi.diff)
    
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
    psi = psi.new,
    peak.y = mod.coef$U1*psi.new[1],
    # psi.sd = NA_real_,
    psi.diff = psi.diff,
    n.iter = it,
    call = cl
  )
  class(z) <- c("break.lm")
  
  return(z)
}

library(lmbreak)
library(ggplot2)
data(SDIpsilo, package = "lmbreak")
for(ID in levels(SDIpsilo$id)){
  SDI.ind <- SDIpsilo %>% filter(id==ID)
  ggplot(SDI.ind, aes(x=time, y=score)) + 
    geom_point()
  tryCatch({
    break101.fit <- SDI101.fit(data = SDI.ind)
  }, error = function(e){
    print(ID)
  }, warning = function(w){
    print(ID)
  })
  
  # print(break101.fit)
}

# does not CV for 2, 4, 15, uniform initializing might be a bit shitty...