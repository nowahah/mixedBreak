### lme.break1.R ---
## * lme.break1 (documentation)
##' @title Mixed-Effects Models with 1 breakpoint estimate
##' @description This function fits a mixed-effects model with a single breakpoint estimate, 
##' in the formulation described in Muggeo (2014), 
##' allowing for random effects on the breakpoint estimate.
##'
##' @param formula a two-sided linear formula object describing both the fixed-effects 
##' and random-effects part of the model, with the response on the left of a ~ operator 
##' and the terms, separated by + operators, on the right. Random-effects terms are 
##' distinguished by vertical bars (|) separating expressions for design matrices from 
##' grouping factors. By default, non-scalar random effects (where the design matrix has 
##' more than one column, e.g. (1+x|f)) are fitted with unstructured (general positive 
##' semidefinite) covariance matrices.
##' @param data a data frame containing variables named in `formula`.
##' @param psi0 [numeric] an optional starting value for the breakpoint estimates (first guess). 
##' It must be a numeric value in the range of the `time` (segmented) variable. Default is the 
##' middle of the observed values for `time` (`mean(range(time))`).
##' @param it.max [integer>0] Maximum number of iterations allowed.
##' @param model,x,y [logicals] If TRUE the corresponding components of the fit 
##' (the model frame, the model matrix, the response) are returned.
##'
##' @return nothing


# On purpose a very specific function to fit segmented relationship with a mixed
# model and estimation a single breakpoint.
library(lme4)
library(rlang) # for some reaons is necessary for methods summary(lmer)

## * lme.break1
##' @export
lme.break1 <- function(
    formula, data, psi0 = NA, tol.logLik = 1e-3, it.max = 10L,
    psi.history = TRUE, model = TRUE, x = FALSE, y = FALSE)
{
  # Input handling ====
  ret.x <- x
  ret.y <- y
  ret.psi <- psi.history
  cl <- match.call() 
  
  yy <- data[,as.character(formula[[2]])] # LHS - model.response
  # TODO - need model.matrix(formula, data) method for making XX with other variate (factors!)
  # XX <- model.matrix(formula, data) # only works for num var, not factors
  # attr(terms.formula(formula), "term.label") # RHS of formula variables
  XX <- data.frame(response=yy, time=data[["time"]], ID=data[["id"]])
  # TODO - need to extract that from formula!!
  n.ind <- nlevels(XX$ID) # change with formula too ?
  n.psi <- 1 # TODO - change for multiple breakpoints version

  # 0. Initialization ====
  
  ## Fix a starting value for the change points
  # should be numeric in the range of segmented variable (time)
  if (!is.numeric(psi0) | min(XX[,"time"]) > psi0 | psi0 > max(XX[, "time"])) {
    psi0 <- mean(range(XX[,"time"])) # uniform initialization for psi0
    warning(sprintf(
      "psi0 not specified: using psi0 = mean(range(time)) = %s for all individuals", 
      format(psi0, digit=1))
    )
  } else if ( !(length(psi0) %in% c(1,n.ind))) {
    warning(sprintf(
      "Only a single breakpoint is allowed for `lme.break1`. \nUsing psi0[1] = %s as starting value",
      psi0[1])
    )
  }
  if (length(psi0)==n.psi) { psi0 <- matrix(rep(psi0, n.ind), ncol = n.psi, byrow = T) }
  
  ## compute the U variate
  XX$U <- pmax(0, XX$time - psi0[XX$ID,1])
  
  ## fit the initial lmm
  mod.init <- suppressWarnings(
    lme4::lmer(response ~ 0 + time + U + (0+time+U|ID), data = XX, REML = TRUE)
  )
  mod.work <- mod.init
  # delta <- mod.init@beta[2:length(mod.init@beta)] # fixed-effects coefficients
  # coef(mod.init) # SUM of fixed and random-effects coefficients
  # t(t(coef(mod.init)$ID) - mod.init@beta) # random-effects coefficients (time, U)
  delta.i <- coef(mod.init)$ID$U
  
  # time-scale transformation for breakpoints
  a1 <- min(XX$time)
  a2 <- max(XX$time)
  logit <- function(psi) return( log((psi-a1)/(a2-psi)) )
  expit <- function(eta) return( (a1 + a2 * exp(eta))/(1+exp(eta)) )
  eta.i <- logit(psi0)[,1] # 0 if default init for psi0 + select first col bcs1 breakpoint
  
  # FOR LOOP
  it <- 0
  logLik.diff <- tol.logLik + 1
  psi.history <- array(NA_real_, dim = c(it.max+1, n.psi, n.ind))
  psi.history[1,,] <- psi0
  while((logLik.diff > tol.logLik | logLik.diff < 0) & it < it.max){
    it <- it + 1
    
    # 1. Compute covariates U, G and O ====
    XX$U <- pmax(0, XX$time - psi.history[it,1,XX$ID])
    XX$V <- -1 * (XX$time > psi.history[it,1,XX$ID])
    XX$D <- ((exp(eta.i) * (a2 - a1)) / (1 + exp(eta.i))^2)[XX$ID]
    XX$G <- delta.i[XX$ID] * XX$V * XX$D
    XX$O <- -eta.i[XX$ID] * XX$G
    
    # 2. Fit the working LMM ====
    mod.work.prev <- mod.work
    mod.work <- suppressWarnings(suppressMessages(
      lme4::lmer(response ~ 0 + time + U + G + (0+time+U+G|ID), data = XX, REML = TRUE)
    )) 
    # message or warnings due to singular fits (induces by plateau setting, cor(time,U~=-1)
    logLik.diff <- REMLcrit(mod.work) - REMLcrit(mod.work.prev)
    # browser()
    
    # 3. Extract the breakpoint linear predictor (fixed+random), update breakpoint estimates ====
    #    extract the difference-in-slopes (fixed+random)
    eta.i <- coef(mod.work)$ID$G
    psi.history[it+1,1,] <- expit(eta.i)
    delta.i <- coef(mod.work)$ID$U
    
  }
  # END FOR
  
  # Return results ====
  # browser()
  summ <- summary(mod.work)
  rownames(summ$coefficients) <- c("time.beta.1", "time.delta.2", "psi.eta")
  z <- list(
    fixed = summ$coefficients,
    random = coef(mod.work)$ID[, 1:(1+n.psi)],
    psi = psi.history[it+1,1,],
    fitted = mod.work@resp$mu,
    logLik = logLik(mod.work),
    REML = REMLcrit(mod.work),
    coefficients = summ$coefficients,
    psi.history = psi.history[1:(it+1),1,],
    n.iter = it,
    call = cl
  )
  if(ret.y) z$y <- XX$yy
  if(model) z$model <- model.frame(mod.work)
  
  class(z) <- append("break.lme1", class(z))
  return(z)
}