### mixed.break1.R ---
## * mixed.break1 (documentation)
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
##' @param x,y [logicals] If TRUE the corresponding components of the fit 
##' (the model matrix, the response) are returned.
##'
##' @return TODO


# On purpose a very specific function to fit segmented relationship with a mixed
# model and estimation a single breakpoint.
library(lme4)
library(rlang) # for some reaons is necessary for methods summary(lmer)

## * mixed.break1
##' @export
mixed.break1 <- function(
    formula, data, psi0 = NA, tol.logLik = 1e-4, it.max = 10L, pattern = "11",
    psi.history = TRUE, x = FALSE, y = FALSE)
{
  if(!(pattern %in% c("11", "10")))
    stop(paste("'pattern' should either be: \n",
               "-'11' for two non-null slopes estimated) ;\n",
               "-'10' for one non-null slope followed by a constant segment"))
  
  # Input handling ====
  ret.x <- x
  ret.y <- y
  ret.psi <- psi.history
  cl <- match.call()
  vars <- list()
  vars$response <- as.character(formula[[2]])
  form.label <- attr(terms(formula), "term.labels")
  vars$fixed <- form.label[!stringr::str_detect(form.label, stringr::fixed("|"))]
  if(!attr(terms(formula), "intercept")){ 
    vars$fixed <- c(vars$fixed, "0") 
  }
  vars$segmented <- vars$fixed[1]
  vars$random <- form.label[stringr::str_detect(form.label, stringr::fixed("|"))]
  vars$random <- strsplit(vars$random, split = " | ", fixed=TRUE)[[1]]
  vars$group <- vars$random[2]
  vars$random <- strsplit(vars$random[1], " + ", fixed = TRUE)[[1]] 
  # one single grouping variable allowed
  
  # stopifnot(exprs = {
  #   is.logical(ret.x),
  #   is.logical(ret.y),
  #   rlang::is_formula(formula),
  # })
  
  # model frame
  XX <- data.frame(
    y = data[[vars$response]],
    time = data[[vars$segmented]],
    data[,setdiff(union(vars$fixed, vars$random), c("0", "-1", vars$segmented))],
    ID = as.factor(data[[vars$group]])
  )
  names(XX)[1] <- vars$response
  names(XX)[names(XX)=="ID"] <- vars$group
  n.ind <- nlevels(XX[[vars$group]])
  n.psi <- 1 # FINAL - change for multiple breakpoints version
  
  # 0. Initialization ====
  
  ## Fix a starting value for the change points
  # should be numeric in the range of segmented variable (time)
  if (!is.numeric(psi0) | 
      any(min(XX[[vars$segmented]]) > psi0) | 
      any(psi0 > max(XX[["time"]]))) {
    psi0 <- mean(range(XX[["time"]])) # uniform initialization for psi0
    message(sprintf(
      "psi0 not specified: using psi0 = mean(range(time)) = %s for all individuals", 
      format(psi0, digit=1))
    )
  } else if ( !(length(psi0) %in% c(1,n.ind))) {
    warning(sprintf(
      "Only a single breakpoint is allowed for `lme.break1`.
      Using psi0[1] = %s instead.", psi0[1])
    )
    psi0 <- psi0[1]
  }
  if (length(psi0)==n.psi) { psi0 <- matrix(rep(psi0, n.ind), ncol = n.psi, byrow = T) }
  
  ## compute the U variate
  XX$U1 <- pmax(0, XX[[vars$segmented]] - psi0[XX[[vars$group]],1])
  vars$fixed <- c(vars$fixed, "U1")
  vars$random <- c(vars$random, "U1")
  if(pattern=="10"){  # PLATEAU
    XX$U1 <- XX[[vars$segmented]] - XX$U1
    vars$fixed <- setdiff(vars$fixed, vars$segmented)
    vars$random <- setdiff(vars$fixed, vars$segmented)
  }  
  # FINAL - rename U1-like variables to show segment number clearly: delta.time.2 / delta.2 ?
  
  ## fit the initial lmm
  formula.init <- as.formula(paste(
    vars$response, "~", 
    paste(vars$fixed, collapse=" + "), "+ (",
    paste(paste(vars$random, collapse=" + "), vars$group, sep = " | "), ")")
  )
  mod.init <- suppressWarnings(
    lme4::lmer(formula.init, data = XX, REML = TRUE)
  )
  mod.work <- mod.init
  vars$fixed <- c(vars$fixed, "G1")
  vars$random <- c(vars$random, "G1")
  # FINAL - rename G-like variables to show segment / breakpoint index clearly: psi.1 ?
  
  # delta <- mod.init@beta[2:length(mod.init@beta)] # fixed-effects coefficients
  # coef(mod.init) # SUM of fixed and random-effects coefficients
  # t(t(coef(mod.init)[[vars$group]]) - mod.init@beta) # random-effects coefficients (time, U)
  delta.i <- coef(mod.init)[[vars$group]]$U1 # FINAL - U1 name
  
  # time-scale transformation for breakpoints
  a1 <- min(XX[[vars$segmented]])
  a2 <- max(XX[[vars$segmented]])
  logit <- function(psi) return( log((psi-a1)/(a2-psi)) )
  expit <- function(eta) return( (a1 + a2 * exp(eta))/(1+exp(eta)) )
  eta.i <- logit(psi0)[,1] # 0 if default init for psi0 + FINAL select col1 bcs 1 breakpoint
  
  # FOR LOOP
  it <- 0
  logLik.diff <- tol.logLik + 1
  psi.history <- array(NA_real_, dim = c(it.max+1, n.psi, n.ind))
  psi.history[1,,] <- psi0
  formula.it <- as.formula(paste(
    vars$response, "~", 
    paste(vars$fixed, collapse=" + "), "+ (",
    paste(paste(vars$random, collapse=" + "), vars$group, sep = " | "), ")")
  )
  while((logLik.diff > tol.logLik | logLik.diff < 0) & it < it.max){
    it <- it + 1
    
    # 1. Compute covariates U, G and O ====
    XX$U1 <- pmax(0, XX[[vars$segmented]] - psi.history[it,1,XX[[vars$group]]])
    XX$V1 <- -1 * (XX[[vars$segmented]] > psi.history[it,1,XX[[vars$group]]])
    XX$D1 <- ((exp(eta.i) * (a2 - a1)) / (1 + exp(eta.i))^2)[XX[[vars$group]]]
    XX$G1 <- delta.i[XX[[vars$group]]] * XX$V1 * XX$D1
    if(pattern=="10") { # PLATEAU
      XX$U1 <- XX[[vars$segmented]] - XX$U1
      XX$G1 <- (-1)*XX$G1
    }
    XX$Off1 <- -eta.i[XX[[vars$group]]] * XX$G1 # offset term
    XX[[vars$response]] <- data[[vars$response]] - XX$Off1
    # FINAL - RENAME U, V, G, O - like variables

    # 2. Fit the working LMM ====
    mod.work.prev <- mod.work
    mod.work <- suppressWarnings(suppressMessages(
      lme4::lmer(formula.it, data = XX, REML = TRUE)
    )) 
    # message or warnings due to singular fits (induces by plateau setting, cor(time,U~=-1)
    logLik.diff <- (REMLcrit(mod.work) - REMLcrit(mod.work.prev))/(REMLcrit(mod.work.prev)+.1)
    
    # 3. Extract the breakpoint linear predictor (fixed+random), update breakpoint estimates ====
    #    extract the difference-in-slopes (fixed+random)
    eta.i <- coef(mod.work)[[vars$group]]$G1
    psi.history[it+1,1,] <- expit(eta.i)
    delta.i <- coef(mod.work)[[vars$group]]$U1
    
  }
  # END FOR
  
  # Return results ====
  summ <- summary(mod.work)
  if(pattern=="10") { # PLATEAU
    rownames(summ$coefficients) <- c(vars$segmented, "(logit.)break.1") # PSILINK 
  } else {
    rownames(summ$coefficients) <- c(
      paste0(vars$segmented, c(".beta.1", paste0(".delta.", 2:nchar(pattern)))), 
      "(logit.)break.1" # PSILINK 
    ) 
  }
  # removing the offset term
  mod.work@frame[[vars$response]] <- data[[vars$response]][!is.na(data[[vars$response]])]
  z <- list(
    lme.fit = mod.work,
    fixed = summ$coefficients,
    random = coef(mod.work)[[vars$group]][, 1:(1+n.psi)],
    psi.i = psi.history[it+1,1,],
    fitted = mod.work@resp$mu,
    off = XX$Off[!is.na(data[[vars$response]])],
    logLik = logLik(mod.work),
    REML = REMLcrit(mod.work),
    psi.history = psi.history[1:(it+1),1,],
    n.iter = it,
    call = cl,
    var.name = vars,
    pattern = pattern
  )
  if(ret.y) z$y <- XX[[vars$response]]
  
  model <- model.frame(mod.work)
  if(pattern=="10"){
    model[[vars$segmented]] <- data[[vars$segmented]][!is.na(data[[vars$response]])]
  }
  z$model <- model
  
  class(z) <- append("mixedBreak1", class(z)) 
  return(z)
}
