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
    formula, data, pattern = NA, psi0 = NA, psi.range = NULL,
    tol.logLik = 1e-4, it.max = 10L, 
    psi.history = TRUE, x = FALSE, y = FALSE, dev.step = 0)
{
  if(!(pattern %in% c("11", "10", "111", "101")))
    stop(paste(
      "'pattern' should either be one of four options:",
      "c(\"11\", \"10\", \"111\", \"101\").\n  See documentation for signification."))
  
  if(!is.null(psi.range) & !any(c("q.prob", "value") %in% names(psi.range)))
    stop(paste(
      "When specified, `psi.range` must have at least one of:",
      "c(\"q.prob\", \"value\") as name."
    ))
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
  # browser()
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
  n.obs <- nrow(XX)
  plateau <- (substr(pattern, 2L, 2L) == "0")

  # 0. Initialization ====
  
  ## Fix a starting value for the change points
  psi0.ok <- TRUE
  # should be numeric in the range of segmented variable (time)
  if (!is.numeric(psi0) | 
      any(min(XX[[vars$segmented]]) > psi0) | 
      any(max(XX[[vars$segmented]]) < psi0)) {
    psi0.ok <- FALSE
  } else if ( !(length(psi0) %in% c(n.psi ,n.ind, n.psi*n.ind))) {
    warning(sprintf(
      "`length(psi0)` doesn't match any of n.psi (%s), n.ind (%s), n.psi * n.ind (%s)", 
      n.psi ,n.ind, n.psi*n.ind)
    )
    psi0.ok <- FALSE
  }
  
  # uniform initialization for psi0
  if(!psi0.ok){
    psi0 <- unname(quantile(XX[[vars$segmented]], probs = 1:n.psi/(n.psi+1)))
    message(paste0(
      "psi0 incorrectly (not numeric / out of range of `", vars$segmented,
      "`) or not specified: \n  Using uniform initialization for all ",
      "individuals: psi0 = c(", paste0(format(psi0, digit=1), collapse = ", "), ")"
    ))
  }
  if (length(psi0)==n.psi) psi0 <- matrix(rep(psi0, n.ind), ncol = n.psi, byrow = T)
  
  
  ## compute the U variate
  # Here and in the rest of this function, variable(s) U are named only after
  # the letter, but may be declined into multiple sub-columns if n.psi > 1
  # Same approach is used ofr deriving U, V, D, G, O variables in main loop
  XX$U <- pmax(matrix(0, ncol = n.psi, nrow = n.obs), 
               XX[[vars$segmented]] - psi0[XX[[vars$group]],])
  vars$fixed <- c(vars$fixed, "U")
  vars$random <- c(vars$random, "U")
  
  if(plateau){
    XX$U[,1] <- XX[[vars$segmented]] - XX$U[,1]
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
  
  delta.i <- coef(mod.init)[[vars$group]]
  delta.i <- delta.i[,colnames(delta.i) %in% paste0("U", c("", 1:n.psi))]
  delta.i <- unname(as.matrix(delta.i))
  # sum of fixed & random coef
  
  # time-scale transformation for breakpoints
  # browser()
  if(is.null(psi.range)){
    # default 10-90%
    psi.range <- unname(quantile(XX[[vars$segmented]], probs = c(1, 9)/10))
  } else if ("q.prob" %in% names(psi.range)) {
    psi.range <- unname(quantile(XX[[vars$segmented]], probs = psi.range[["q.prob"]]))
  } else if ("value" %in% names(psi.range)) {
    psi.range <- psi.range[["value"]]
    stopifnot( length(psi.range)==2 & is.numeric(psi.range) )
    if(psi.range[1] < min(XX[[vars$segmented]]) | 
       psi.range[2] > max(XX[[vars$segmented]])){
      
      psi.range[1] <- max(psi.range[1], min(XX[[vars$segmented]]))
      psi.range[2] <- min(psi.range[2], max(XX[[vars$segmented]]))
      
      message(sprintf(paste(
        "`psi.range` specifies values out of observed segmented variable `%s` range.\n",
        "  Projecting `psi.range` onto observed range: `psi.range` = c(%s)"
        ), vars$segmented, paste(psi.range, collapse = ", "))
      )
    }
  }
  
  logit <- function(psi) return( log((psi-psi.range[1])/(psi.range[2]-psi)) )
  expit <- function(eta) return( (psi.range[1] + psi.range[2] * exp(eta))/(1+exp(eta)) )
  eta.i <- logit(psi0)
  XX$psi <- psi0[XX[[vars$group]],] # creates multiples col if n.psi > 2
  XX$eta <- eta.i[XX[[vars$group]],]
  XX$delta <- delta.i[XX[[vars$group]],]
  # browser()
  
  # FOR LOOP
  it <- 0
  logLik.diff <- tol.logLik + 1
  psi.history <- array(NA_real_, dim = c(it.max+1, n.psi, n.ind))
  psi.history[1,,] <- t(psi0)
  vars$fixed <- c(vars$fixed, "G") 
  vars$random <- c(vars$random, "G")
  formula.it <- as.formula(paste(
    vars$response, "~", 
    paste(vars$fixed, collapse=" + "), "+ (",
    paste(paste(vars$random, collapse=" + "), vars$group, sep = " | "), ")")
  )
  while((logLik.diff > tol.logLik | logLik.diff < 0) & it < it.max){
    it <- it + 1
    
    # 1. Compute covariates U, G and O ====
    XX$U <- pmax(matrix(0, ncol = n.psi, nrow = n.obs), XX[[vars$segmented]] - XX$psi)
    XX$V <- -1 * (XX[[vars$segmented]] > XX$psi)
    XX$D <- (exp(XX$eta) * diff(psi.range)) / (1 + exp(XX$eta))^2
    XX$G <- matrix(XX$delta * XX$V * XX$D, ncol = n.psi, nrow = n.obs)
    
    # mixed.break.1 code
    # XX$U1 <- pmax(0, XX[[vars$segmented]] - psi.history[it,1,XX[[vars$group]]])
    # XX$V1 <- -1 * (XX[[vars$segmented]] > psi.history[it,1,XX[[vars$group]]])
    # XX$D1 <- ((exp(eta.i) * (psi.range[2] - psi.range[1])) / (1 + exp(eta.i))^2)[XX[[vars$group]]]
    # XX$G1 <- delta.i[XX[[vars$group]]] * XX$V1 * XX$D1
    
    if(plateau) {
      XX$U[,1] <- XX[[vars$segmented]] - XX$U[,1]
      XX$G[,1] <- (-1)*XX$G[,1]
    }
    
    # XX$Off1 <- -eta.i[XX[[vars$group]]] * XX$G1 # offset term
    # XX[[vars$response]] <- data[[vars$response]] - XX$Off1
    
    XX$O <- -XX$eta * XX$G
    XX$OFF <- rowSums(XX$O)
    XX[[vars$response]] <- data[[vars$response]] - XX$OFF
    
    # 2. Fit the working LMM ====
    mod.work.prev <- mod.work
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
    eta.i <- coef(mod.work)[[vars$group]]
    eta.i <- as.matrix(
      eta.i[,colnames(eta.i) %in% paste0("G", c("", 1:n.psi))], ncol = n.psi
    )
    XX$eta <- eta.i[XX[[vars$group]],]
    XX$eta <- unname(XX$eta)
    
    XX$psi <- expit(XX$eta)
    psi.history[it+1,,] <- t(expit(eta.i))
    
    delta.i <- coef(mod.work)[[vars$group]]
    delta.i <- as.matrix(
      delta.i[,colnames(delta.i) %in% paste0("U", c("", 1:n.psi))], ncol = n.psi
    )
    XX$delta <- delta.i[XX[[vars$group]],]
    XX$delta <- unname(as.matrix(XX$delta))
    
  }
  # END FOR
  # print warning message at convergence if any
  if(it %in% names(warn.list)) {
    warning("Warning at convergence during last LMM estimation: \n  ", warn.list[[it]])
  }
  # warning in case of non convergence (it.max iterations reached)
  if(it==it.max){
    warning(paste0("Maximum number of iterations allowed (`it.max` = ",it.max, 
                   ") reached: possible non-convergence. \n",
                   "  Consider increasing the value of `it.max` at least once,",
                   " depending on value of `logLik.diff`."))
  }
  
  # At convergence, does delta*(eta-tilde(eta)) approaches 0 ?
  XX$VD <- XX$V * XX$D
  XX[[vars$response]] <- data[[vars$response]]
  formula.check <- gsub("G(\\d?)", "VD\\1", as.character(formula.it))
  formula.check <- as.formula(paste(formula.check[2], formula.check[1], formula.check[3]))
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
  
  
  # Return results ====
  # browser()
  summ <- summary(mod.work)
  if(plateau) { # PSILINK - both
    rownames(summ$coefficients) <- if(n.psi==1) {
      c(vars$segmented, "(logit.)break.") 
    } else {
      c(paste0(vars$segmented, c(".beta.1", paste0(".delta.", 3:(n.psi+1)))),
        paste0("(logit.)break.", 1:n.psi)) 
    }
    
  } else {
    rownames(summ$coefficients) <- if(n.psi==1) {
      c(paste0(vars$segmented, c(".beta.1", ".delta.2")), "(logit.)break.1") 
    } else {
      c(paste0(vars$segmented, c(".beta.1", paste0(".delta.", 2:(n.psi+1)))),
        paste0("(logit.)break.", 1:n.psi)) 
    } 
  }
  # removing the offset term
  mod.work@frame[[vars$response]] <- data[[vars$response]][!is.na(data[[vars$response]])]
  z <- list(
    lme.fit = mod.work,
    fixed = summ$coefficients,
    random = coef(mod.work)[[vars$group]][, 1:(1+n.psi)],
    psi.i = t(psi.history[it+1,,]),
    fitted = mod.work@resp$mu,
    off = XX$OFF[!is.na(data[[vars$response]])],
    logLik = logLik(mod.work),
    REML = REML.work,
    psi.history = psi.history[1:(it+1),,],
    psi.range = psi.range,
    n.iter = it,
    logLik.diff = logLik.diff,
    call = cl,
    var.name = vars,
    pattern = pattern,
    warn.list = warn.list,
    lme.fit.check = mod.check
  )
  if(ret.y) z$y <- XX[[vars$response]]
  
  model <- model.frame(mod.work)
  if(plateau){
    model[[vars$segmented]] <- data[[vars$segmented]][!is.na(data[[vars$response]])]
  }
  z$model <- model
  
  class(z) <- append(paste0("mixedBreak", c("", n.psi)), class(z)) 
  # z.it <<- append(z.it, list(z)); plot(z)
  return(z)
}
