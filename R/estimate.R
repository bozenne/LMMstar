### estimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 20 2021 (23:25) 
## Version: 
## Last-Updated: sep 20 2021 (11:00) 
##           By: Brice Ozenne
##     Update #: 157
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * estimate (documentation)
##' @title Optimizer for mixed models
##' @description Optimization procedure for mixed model (REML or ML).
##' Alternate between one step of gradient descent to update the variance parameters
##' and a GLS estimation of the mean parameters.
##'
##' @param design [list] information about the mean and variance structure. Obtained using \code{.model.matrix.lmm}.
##' @param time [list] information about the time variable (e.g. unique values)
##' @param method.fit [character] Should Restricted Maximum Likelihoood (\code{"REML"}) or Maximum Likelihoood (\code{"ML"}) be used to estimate the model parameters?
##' @param type.information [character] Should the expected information be computed  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param transform.sigma,transform.k,transform.rho possible transformations for the variance parameters.
##' @param precompute.moments [logical] have certain key quantities be pre-computed (e.g. \eqn{X'X}).
##' @param init [numeric vector] values used to initialize the mean and parameters.
##' @param n.iter [integer,>0] maximum number of iterations.
##' @param tol.score [double,>0] convergence is not reached unless each element of the score is smaller (in absolute value) than this value. 
##' @param tol.param [double,>0] convergence is not reached unless the change in parameter values between two iterations is smaller (in absolute value) than this value. 
##' @param trace [1, 2, or 3] should each iteration be displayed?
##' 
##' @examples
##' \dontrun{
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit Linear Mixed Model
##' LMMstar.options(optimizer = "gls")
##' eUN.gls <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##' 
##' LMMstar.options(optimizer = "FS")
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##' 
##' }
##'
##' @keywords internal

## * estimate (code)
.estimate <- function(design, time, method.fit, type.information,
                      transform.sigma, transform.k, transform.rho,
                      precompute.moments, init, n.iter, tol.score, tol.param, trace){

   
    if(!precompute.moments){
        stop("Only implemented when option \'precompute.moments\' is TRUE")
    }
    options <- LMMstar.options()
    if(is.null(n.iter)){
        n.iter <- options$param.optimizer["n.iter"]
    }
    if(is.null(tol.score)){
        tol.score <- options$param.optimizer["tol.score"]
    }
    if(is.null(tol.param)){
        tol.param <- options$param.optimizer["tol.param"]
    }
    if(is.null(trace)){
        trace <- FALSE
    }

    ## ** prepare
    param.mu <- design$param$mu
    param.sigma <- design$param$sigma
    param.k <- design$param$k
    param.rho <- design$param$rho
    param.Omega <- c(param.sigma,param.k,param.rho)
    param.name <- c(param.mu,param.Omega)
    n.param <- length(param.name)
    Upattern <- design$vcov$X$Upattern
    n.Upattern <- NROW(Upattern)
    precompute.XY <- design$precompute.XY
    precompute.XX <- design$precompute.XX$pattern
    key.XX <- design$precompute.XX$key

    ## ** intialization
    if(is.null(init)){

        param.value <- stats::setNames(rep(NA, n.param),param.name)

        ## mean value
        start.OmegaM1 <- stats::setNames(lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 2
            diag(1, nrow = length(Upattern[iPattern,"time"][[1]]), ncol = length(Upattern[iPattern,"time"][[1]]))
        }), Upattern$name)
        param.value[param.mu] <- .estimateGLS(OmegaM1 = start.OmegaM1, pattern = Upattern$name, precompute.XY = precompute.XY, precompute.XX = precompute.XX, key.XX = key.XX)
        
        ## vcov values
        iResiduals.long <- design$Y - design$mean %*% param.value[param.mu]
        outInit <- .initialize(design$vcov, residuals = iResiduals.long)
        param.value[names(outInit)] <- outInit
    }else{
        if(any(param.name %in% names(init) == FALSE)){
            stop("Initialization does not contain value for all parameters. \n",
                 "Missing parameters: \"",paste(param.name[param.name %in% names(init) == FALSE], collapse = "\" \""),"\". \n")
        }
        param.value <- init[param.name]
    }
    ## microbenchmark(test =     .estimateGLS(OmegaM1 = start.OmegaM1, pattern = pattern, precompute.XY = precompute.XY, precompute.XX = precompute.XX, key.XX = key.XX),
    ##                lm.fit = lm.fit(y = ncgsL$cholest[!is.na(ncgsL$cholest)], x = model.matrix(~time+highdose.time, data=ncgsL)[!is.na(ncgsL$cholest),]),
    ##                lm = lm(cholest~time+highdose.time, data=ncgsL)
    ##                )


    if(trace>1){
        cat("\nInitialization:\n")
        print(param.value)
    }
    
    ## ** loop
    cv <- FALSE
    param.valueM1 <- NULL
    if(trace>1){
        cat("\nLoop:\n")
    }

    for(iIter in 1:n.iter){ ## iIter <- 1
        param.valueM1 <- param.value
        
        ## *** moments
        outMoments <- .moments.lmm(value = param.value, design = design, time = time, method.fit = method.fit, type.information = "expected",
                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                   logLik = FALSE, score = TRUE, information = TRUE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = c("variance","correlation"), robust = FALSE,
                                   trace = FALSE, precompute.moments = TRUE, transform.names = FALSE)

        if(all(abs(outMoments$score)<tol.score) && (iIter==1 || all(abs(param.valueM1 - param.value)<tol.param))){
            cv <- TRUE
            break
        }

        ## *** variance estimate
        param.newvalue.trans <- outMoments$reparametrize$p + as.double(outMoments$score %*% solve(outMoments$information))
        param.value[param.Omega] <- .reparametrize(param.newvalue.trans,
                                                   type = design$param$type[names(param.newvalue.trans)],
                                                   strata = design$param$strata[names(param.newvalue.trans)],
                                                   time.levels = time$levels, time.k = design$param$time.k, time.rho = design$param$time.rho,
                                                   Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)$p
        

        ## *** mean estimate
        iOmega <- .calc_Omega(object = design$vcov, param = param.value, keep.interim = TRUE)
        param.value[param.mu] <- .estimateGLS(OmegaM1 = stats::setNames(lapply(iOmega, solve), names(iOmega)),
                                              pattern = Upattern$name, precompute.XY = precompute.XY, precompute.XX = precompute.XX,
                                              key.XX = key.XX)

        
        if(trace > 0 && trace < 3){
            cat("*")
        }else if(trace>2){
            M.print <- rbind(estimate = param.value,
                             diff = param.value - param.valueM1,
                             score = c(rep(0, length(param.mu)),outMoments$score))
            rownames(M.print) <- paste0(rownames(M.print),".",iIter)
            print(M.print)
        }
        
    }

    if(trace==1){
        cat("\n")
    }else if(trace>1){
        if(trace == 2){
            cat("\n")
            print(param.value)
        }
        if(cv){
            if(iIter==1){
                cat("Convergence after ",iIter," iteration: max score=",max(abs(outMoments$score)),"\n", sep = "")
            }else{
                cat("Convergence after ",iIter," iterations: max score=",max(abs(outMoments$score))," | max change in coefficient=",max(abs(param.valueM1 - param.value)),"\n", sep = "")
            }
        }else{
            if(iIter==1){
                cat("No convergence after ",iIter," iteration: max score=",max(abs(outMoments$score)),"\n")
            }else{
                cat("No convergence after ",iIter," iterations: max score=",max(abs(outMoments$score))," | max change in coefficient= ",max(abs(param.valueM1 - param.value)),"\n", sep = "")
            }
        }
    }

    ## ** export
    return(list(estimate = param.value,
                previous.estimate = param.valueM1,
                score = outMoments$score,
                n.iter = iIter,
                cv = cv))
}

.estimateGLS <- function(OmegaM1, pattern, precompute.XY, precompute.XX, key.XX){

    name.param <- colnames(key.XX)
    n.param <- length(name.param)
    max.key <- key.XX[n.param,n.param]
    numerator <- matrix(0, nrow = n.param, ncol = 1)
    denominator <- matrix(0, nrow = n.param, ncol = n.param)
    
    for(iPattern in pattern){ ## iPattern <- pattern[1]
        iVec.Omega <- as.double(OmegaM1[[iPattern]])
        iTime2 <- length(iVec.Omega)
        numerator  <- numerator + t(iVec.Omega %*%  matrix(precompute.XY[[iPattern]], nrow = iTime2, ncol = n.param))
        denominator  <- denominator + as.double(iVec.Omega %*%  matrix(precompute.XX[[iPattern]], nrow = iTime2, ncol = max.key))[key.XX]
    }

    out <- solve(denominator) %*% numerator
    return(stats::setNames(as.double(out), name.param))
}


##----------------------------------------------------------------------
### estimate.R ends here
