### mixed.break.R ---
## * mixed.break (documentation)
##' @title Mixed-Effects Models with 1 or 2 breakpoint estimates
##' @description This function fits a mixed-effects model with one or two
##' breakpoint estimates, in the formulation described in Muggeo (2014), 
##' allowing for random effects on the breakpoint estimates.
##'
##' @param formula a two-sided linear formula object describing both the 
##' fixed-effects and random-effects part of the model, with the response on 
##' the left of `~` operator and the terms, separated by + operators, on the 
##' right. Random-effects terms are distinguished by vertical bars (|) 
##' separating expressions for design matrices from grouping factors. 
##' @param data a data frame containing variables named in `formula`.
##' @param psi0 [numeric] an optional starting value for the breakpoint 
##' estimates (first guess). 
##' It must be a numeric value in the range of the segmented variable. 
##' Default is the middle of the observed segmented variable values 
##' (`mean(range(time))`).
##' @param pattern [character] One of the two following options: \itemize{
##' \item '11' for a two-phases relationship with both slopes estimated 
##' (assumed non-null) ;\n
##' \item '10' for a two-phases relationship with a non-null slope followed by
##'  a constant segment (second-phase slope is forced to 0) ;\n
##' \item '111' for a three-phases relationship with all 3 slopes estimated 
##' (assumed non-null) ;\n
##' \item '101' for a three-phases relationship with one non-null slope 
##' followed by a constant segment (second phase slope is forced to 0) and 
##' another non-null slope.
##' }
##' @param it.max [integer>0] Maximum number of iterations allowed. Default to
##' 10. A warning is emitted is the maximum number of iterations is reached.
##' @param x,y [logical] If TRUE the corresponding components of the fit 
##' (the model matrix, the response) are returned.
##'
##' @return TODO


# On purpose a very specific function to fit segmented relationship with a mixed
# model and estimation a single breakpoint. Pattern can be either of four 
# options in c("11", "10", "101", "111").

## * mixed.break
##' @export
mixed.break <- function(
    formula, data, pattern, psi0 = NA, tol.logLik = 1e-4, it.max = 10L, 
    psi.history = TRUE, x = FALSE, y = FALSE, dev.step = 0)
{
  if(!(pattern %in% c("11", "10", "101", "111"))){
    stop(paste(
      "Argument 'pattern' should be one of c(\"11\", \"10\", \"101\", \"111\") \n"
      # "-'111' for three non-null slopes estimated (unconstrained 3 segments) ;\n",
      # "-'101' for one non-null slope followed by a constant segment,",
      # "then another non-null segment"
    ))
  }
  require(dplyr)
  browser()
  
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
  # one single grouping variable assumed
  
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
  n.psi <- nchar(pattern) - 1L

  # 0. Initialization ====
  
  ## Fix a starting value for the change points
  # should be numeric in the range of segmented variable (time)
  psi0.ok <- TRUE
  if (!is.numeric(psi0) | 
      any(min(XX[[vars$segmented]]) > psi0) | 
      any(max(XX[[vars$segmented]]) < psi0)) {
    psi0.ok <- FALSE
  } else if ( !(length(psi0) %in% c(n.psi, n.psi*n.ind))) {
    warning(sprintf(
      paste0("length(psi0) (%o) doesn't match n.psi (%s) nor n.ind * n.psi (%s)"),
      length(psi0), n.psi, n.psi*n.ind)
    )
    psi0.ok <- FALSE
  }
  # catching wrong psi0 initialization
  if(!psi0.ok) {
    # uniform initialization for psi0
    psi0 <- quantile(XX[[vars$segmented]], probs = (1:n.psi)/(n.psi+1))
    message(paste0(
      "psi0 incorrectly (not numeric / out of range of `", vars$segmented,
      "`) or not specified: \nUsing uniform initialization for all ",
      "individuals: psi0 = c(", paste0(format(psi0, digit=1), collapse = ", "), ")"
    ))
  }
  if (length(psi0)==n.psi) psi0 <- matrix(rep(psi0, n.ind), ncol = n.psi, byrow = T)
  # dim(psi0) = c(n.ind, n.psi)
  
  ## compute the U variables
  addU0var <- function(df, psi.index){
    for(i in seq_along(psi.index)){
      df[[paste0("U", i)]] <- pmax(0, df[[vars$segmented]] - psi0[df[[vars$group]],i])
    }
    
    df
  }
  # browser()
  XX <- addU0var(XX, 1:n.psi)
  vars$fixed <- c(vars$fixed, paste0("U", 1:n.psi)) # FORMULA USELESS, defaulted
  vars$random <- c(vars$random, paste0("U", 1:n.psi)) # FORMULA 
  plateau <- strsplit(pattern, split="")[[1]][2] == "0"
  if(plateau){  # PLATEAU
    vars$fixed <- setdiff(vars$fixed, vars$segmented) # keep only plateau block
    vars$random <- setdiff(vars$fixed, vars$segmented)
  }
  
  ## fit the initial lmm
  warn.list <- list()
  formula.init <- as.formula(paste(
    vars$response, "~", 
    paste(vars$fixed, collapse=" + "), "+ (",
    paste(paste(vars$random, collapse=" + "), vars$group, sep = " | "), ")")
  )
  tryCatch(
    mod.init <<- suppressMessages(lme4::lmer(formula.init, data = XX, REML = TRUE)),
    warning = function(w){
      warn.list$init <<- w$message
      mod.init <<- suppressWarnings(suppressMessages(
        lme4::lmer(formula.init, data = XX, REML = TRUE)
      ))
    },
    error = function(e){
      stop(e$message)
    }
  )
  if(dev.step==-1) return(mod.init)
  mod.work <- mod.init
  
  # HERE - weird behaviour with null variance of random effects in the initial mixed model
  delta.i <- coef(mod.init)[[vars$group]] %>% dplyr::select(starts_with("U"))
  XX[,paste0("delta.", 1:n.psi)] <- delta.i[XX[[vars$group]],]
  # FINAL - U1 name - sum of fixed & random coef
  
  # time-scale transformation for breakpoints (unconstrained)
  break.range <- quantile(data[[vars$segmented]][!is.na(data[[vars$response]])], c(1,9)/10)
  a1 <- unname(break.range[1])
  a2 <- unname(break.range[2])
  logit <- function(psi) return( log((psi-a1)/(a2-psi)) )
  expit <- function(eta) return( (a1 + a2 * exp(eta))/(1+exp(eta)) )
  eta.i <- logit(psi0)
  XX <- cbind(XX, eta.i[XX[[vars$group]],])
  names(XX)[(ncol(XX)-n.psi+1):ncol(XX)] <- paste0("eta.", 1:n.psi)
  XX <- cbind(XX, psi0[XX[[vars$group]],])
  names(XX)[(ncol(XX)-n.psi+1):ncol(XX)] <- paste0("psi.", 1:n.psi)
  
  # FOR LOOP
  it <- 0
  logLik.diff <- tol.logLik + 1
  psi.history <- array(NA_real_, dim = c(it.max+1, n.psi, n.ind))
  psi.history[1,,] <- t(psi0)
  
  vars$fixed <- c(vars$fixed, paste0("G", 1:n.psi)) # FORMULA
  vars$random <- c(vars$random, paste0("G", 1:n.psi))# FORMULA
  # FINAL - rename G-like variables to show breakpoint index:  break.1 ?
  formula.it <- as.formula(paste(
    vars$response, "~", 
    paste(vars$fixed, collapse=" + "), "+ (",
    paste(paste(vars$random, collapse=" + "), vars$group, sep = " | "), ")")
  )
  
  # 0. Update covariates U, G and O ====
  # assuming we have current breakpoints in XX
  # SLOW - FOR LOOP
  updateUGO <- function(df, delta, plateau){
    # browser() # HERE
    for(n in 1:n.psi){
      psi.n <- df[[paste0("psi.", n)]]
      df[[paste0("U", n)]] <- pmax(0, df[[vars$segmented]] - psi.n)
      
      df[[paste0("V", n)]] <- (-1) * (df[[vars$segmented]] > psi.n)
      
      eta.n <- df[[paste0("eta.", n)]]
      df[[paste0("D", n)]] <- (exp(eta.n) * (a2 - a1)) / (1 + exp(eta.n))^2
      
      df[[paste0("G", n)]] <- df[[paste0("V", n)]] * df[[paste0("D", n)]] * 
        df[[paste0("delta.", n)]]
      
      df[[paste0("Off", n)]] <- (-1) * df[[paste0("G", n)]] * df[[paste0("eta.", n)]]
    }
    
    # 10 / 101 plateau changes only for these patterns
    # FINAL - adapt to any kind of pattern
    if(plateau){
      df$U1 <- df[[vars$segmented]] - df$U1
      df$G1 <- -df$G1
    }
    
    df$OFF <- apply(df %>% select(paste0("Off", 1:n.psi)), 1, sum)
    
    df
  }
  
  while((logLik.diff > tol.logLik | logLik.diff < 0) & it < it.max){
    it <- it + 1
    
    # 1. Compute covariates U, G and Off ====
    XX <- updateUGO(XX, delta.i, plateau)
    XX[[vars$response]] <- data[[vars$response]] - XX$OFF
    # browser()
    
    # 2. Fit the working LMM ====
    mod.work.prev <- mod.work
    # browser()
    tryCatch(
      mod.work <- suppressMessages(lme4::lmer(formula.it, data = XX, REML = TRUE)),
      warning = function(w){
        warn.list[[it]] <<- w$message
        names(warn.list)[length(warn.list)] <<- it
        mod.work <<- suppressWarnings(suppressMessages(
          lme4::lmer(formula.it, data = XX, REML = TRUE)
        ))
      },
      error = function(e){
        stop("Error during LMM estimation: ", e$message)
      }
    ) 
    # message or warnings due to singular fits or non CV
    # (induces by plateau setting, cor(time,U~=-1), or null random effects
    REML.work <- lme4::REMLcrit(mod.work)
    REML.prev <- lme4::REMLcrit(mod.work.prev)
    logLik.diff <- (REML.work - REML.prev)/(REML.work+.1)
    
    # 3. Extract the breakpoint linear predictor (fixed+random) ===
    #    update breakpoint estimates
    #    extract the difference-in-slopes (fixed+random)
    # browser()
    eta.i <- coef(mod.work)[[vars$group]][,paste0("G", 1:n.psi)]
    if(n.psi==1) eta.i <- matrix(eta.i, nrow = n.ind)
    XX[,paste0("eta.", 1:n.psi)] <- eta.i[XX[[vars$group]],]
    psi.history[it+1,,] <- t(expit(eta.i))
    
    delta.i <- coef(mod.work)[[vars$group]][,paste0("U", 1:n.psi)]
    if(n.psi==1) delta.i <- matrix(delta.i, nrow = n.ind)
    XX[,paste0("delta.", 1:n.psi)] <- delta.i[XX[[vars$group]],]
    
    # browser()
  }
  # END FOR
  # print warning message at convergence if any
  if(it %in% names(warn.list)) {
    warning("Warning at convergence during last LMM estimation: \n    ", 
            warn.list[[it]])
  }
  # warning in case of non convergence (it.max iterations reached)
  if(it==it.max){
    warning(paste0("Maximum number of iterations allowed (`it.max` = ",it.max, 
                   ") reached: possible non-convergence. \n",
                   "  Consider increasing `it.max` at least once,",
                   " depending on final value of `logLik.diff`."))
  }
  
  # At convergence, does delta*(eta-tilde(eta)) approaches 0 ?
  # browser()
  XX[,paste0("VD", 1:n.psi)] <- XX[,paste0("V", 1:n.psi)] * XX[,paste0("D", 1:n.psi)]
  XX[[vars$response]] <- data[[vars$response]]
  formula.check <- stringr::str_replace_all(
    as.character(formula.it), 
    c("G" = "VD")
  )
  formula.check <- as.formula(paste(
    formula.check[2], formula.check[1], formula.check[3]
  ))
  tryCatch(
    mod.check <- suppressMessages(
      lme4::lmer(formula.check, data = XX, REML = TRUE,
                 control = lme4::lmerControl(optimizer = "nlminbwrap"))
    ),
    warning = function(w){
      warn.list$check <<- w$message
      mod.check <<- suppressWarnings(suppressMessages(
        lme4::lmer(formula.check, data = XX, REML = TRUE,
                   control = lme4::lmerControl(optimizer = "nlminbwrap"))
      ))
    },
    error = function(e){
      stop("Error in checking LMM estimation: ", e$message)
    }
  ) 
  # summary(mod.check) # significance of VD coefficients: VD1 yes (bad), VD2 no (good)
  
  
  # Return results ====
  # browser()
  summ <- summary(mod.work)
  if(plateau) { # PLATEAU
    rownames(summ$coefficients) <- stringr::str_replace_all(
      rownames(summ$coefficients), 
      c("U1" = paste0(vars$segmented, ".beta.1"),
        "U" = paste0(vars$segmented, ".delta."), 
        "G" = "logit.break.")
    )
    # PSILINK
  } else {
    # FORMULA - intercept or not
    rownames(summ$coefficients)[rownames(summ$coefficients) == vars$segmented] <- 
      paste0(vars$segmented, ".beta.1")
    rownames(summ$coefficients) <- stringr::str_replace_all(
      rownames(summ$coefficients), 
      c("U" = paste0(vars$segmented, ".delta."), 
        "G" = "logit.break.")
    )
    # PSILINK
  }
  # browser()
  # removing the offset term
  mod.work@frame[[vars$response]] <- data[[vars$response]][!is.na(data[[vars$response]])]
  z <- list(
    lme.fit = mod.work,
    fixed = summ$coefficients,
    random = coef(mod.work)[[vars$group]], # FORMULA
    psi.i = psi.history[it+1,,],
    fitted = mod.work@resp$mu,
    off = XX$OFF[!is.na(data[[vars$response]])],
    logLik = logLik(mod.work),
    REML = REML.work,
    psi.history = psi.history[1:(it+1),,],
    psi.range = break.range,
    n.iter = it,
    logLik.diff = logLik.diff,
    call = cl,
    var.name = vars,
    pattern = pattern,
    warn.list = warn.list,
    lme.fit.check = mod.check
  )
  if(n.psi > 1){
    z$psi.i <- data.frame(
      break. = t(psi.history[it+1,,]), 
      order = t(diff(psi.history[it+1,,]) > 0)
    )
    colnames(z$psi.i) <- c(
      paste0("break.", 1:n.psi), 
      paste0("order.", paste0(1:(n.psi-1), 2:n.psi))
    ) #FINAL - here is only for 2 breakpoints
  }
  
  if(ret.y) z$y <- XX[[vars$response]]
  
  model <- model.frame(mod.work)
  if(plateau){
    model[[vars$segmented]] <- data[[vars$segmented]][!is.na(data[[vars$response]])]
  }
  z$model <- model
  
  class(z) <- append(c("mixedBreak", paste0("mixedBreak", n.psi)), class(z)) 
  return(z)
}
