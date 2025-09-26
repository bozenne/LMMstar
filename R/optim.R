### optim.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 17 2025 (11:27) 
## Version: 
## Last-Updated: sep 26 2025 (16:45) 
##           By: Brice Ozenne
##     Update #: 35
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .optim (documentation)
##' @title Optimizer for mixed models
##' @description Optimization procedure for mixed model (REML or ML).
##' Alternate between one step of gradient descent to update the variance parameters
##' and a GLS estimation of the mean parameters.
##' @noRd
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
##' @param init.cor [1,2] method to initialize the correlation parameters.
##' @param trace [1, 2, or 3] should each iteration be displayed?
##' 
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit Linear Mixed Model
##' LMMstar.options(optimizer = "BFGS")
##' eUN.gls <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##' 
##' LMMstar.options(optimizer = "FS")
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##'
##' @keywords internal

## * .optim (code)
.optim <- function(design, time, method.fit, type.information, 
                   transform.sigma, transform.k, transform.rho,
                   precompute.moments, optimizer, init, n.iter, tol.score, tol.param, n.backtracking, init.cor, trace){

    ## ** apply default values
    if(is.null(trace)){
        trace <- FALSE
    }

    ## ** prepare
    index.cluster <- design$index.cluster

    param.mu <- design$param[design$param$type=="mu","name"]
    param.sigma <- design$param[design$param$type=="sigma","name"]
    param.k <- design$param[design$param$type=="k","name"]
    param.rho <- design$param[design$param$type=="rho","name"]
    param.Omega <- c(param.sigma,param.k,param.rho)
    param.name <- c(param.mu,param.Omega)

    param.fixed <- design$param[!is.na(design$param$constraint),"name"]
    param.Omega2 <- setdiff(param.Omega,param.fixed)
    param.mu2 <- setdiff(param.mu,param.fixed)
    param.fixed.mu <- setdiff(param.mu,param.mu2)
    design.param2 <- design$param[match(param.Omega, design$param$name),,drop=FALSE]
    
    n.param <- length(param.name)
    Upattern <- design$vcov$Upattern
    n.Upattern <- NROW(Upattern)

    if(length(param.fixed.mu)>0){
        partialHat <- design$mean[,param.fixed.mu,drop=FALSE] %*% init[param.fixed.mu]
        partialY <- design$Y - partialHat
        
        if(!is.null(design$precompute.XX) && !is.null(design$precompute.XY) && length(param.mu2)>0){
            ## remove XX involving fixed parameters
            key.XX <- design$precompute.XX$key[param.mu2,param.mu2,drop=FALSE]
            index.precompute <- unique(as.double(design$precompute.XX$key[param.mu2,param.mu2]))
            key.XX[] <- stats::setNames(1:length(index.precompute),index.precompute)[as.character(design$precompute.XX$key[param.mu2,param.mu2])]        
            precompute.XX <- lapply(design$precompute.XX$pattern, function(iX){iX[,index.precompute,drop=FALSE]})


            ## update XY with X(Y-Xb(fixed))
            if(attr(design$weights,"user-defined")){
                wX.mean <- sweep(design$mean[,param.mu2,drop=FALSE], FUN = "*", MARGIN = 1, STATS = sqrt(design$weights[attr(index.cluster,"vectorwise")]))
                wY <- cbind(partialY*sqrt(design$weights[attr(index.cluster,"vectorwise")]))
            }else{
                wX.mean <- design$mean[,param.mu2,drop=FALSE]
                wY <- partialY
            }
            precompute.XY <- .precomputeXR(X = wX.mean,
                                           residuals = wY,
                                           pattern = design$vcov$Upattern$name,
                                           pattern.ntime = stats::setNames(design$vcov$Upattern$n.time, design$vcov$Upattern$name),
                                           pattern.cluster = attr(design$vcov$pattern,"list"), index.cluster = design$index.cluster)

        }else{
            precompute.XY <- NULL
            precompute.XX <- NULL
            key.XX <- NULL
        }
        
    }else{
        partialY <- design$Y
        precompute.XY <- design$precompute.XY
        precompute.XX <- design$precompute.XX$pattern
        key.XX <- design$precompute.XX$key
    }
    
    effects <- c("mean","variance","correlation")
    
    ## ** intialization
    if(!is.null(init)){

        if(is.matrix(init)){
            if(NROW(init)!=time$n || NCOL(init)!=time$n){
                stop("When a matrix, initialization should be a square matrix with dimensions compatible with time (",time$n,"x",time$n,"). \n")
            }
            init.mu <- NULL
            init.Omega <- init
            init <- NULL
        }else if(!is.vector(init)){
            stop("Initialization should either be a vector containing value for all, or only mean parameters, \n",
                 "or be a full data variance-covariance matrix. \n")
        }else if(any(param.mu2 %in% names(init) == FALSE)){
            stop("Initialization does not contain value for all mean parameters. \n",
                 "Missing parameters: \"",paste(param.mu[param.mu %in% names(init) == FALSE], collapse = "\" \""),"\". \n")
        }else if(all(names(init) %in% param.mu)){
            init.mu <- init
            init.Omega <- NULL
            init <- NULL
        }else if(any(param.Omega %in% names(init) == FALSE)){
            stop("Initialization does not contain value for all variance-covariance parameters. \n",
                 "Missing parameters: \"",paste(param.Omega[param.Omega %in% names(init) == FALSE], collapse = "\" \""),"\". \n")
        }

    }else{
        init.mu <- NULL
        init.Omega <- NULL
    }

    if(is.null(init)){

        param.value <- stats::setNames(rep(NA, n.param),param.name)

        ## mean value
        if(!is.null(init.mu)){
            param.value[param.mu2] <- init.mu[param.mu2]
        }else if(length(param.mu2)>0){

            if(!is.null(init.Omega)){
                start.OmegaM1 <- stats::setNames(lapply(Upattern$name, function(iPattern){ ## iPattern <- 1
                    iCluster <- attr(design$vcov$pattern,"list")[[iPattern]][1]
                    iTime <- design$index.clusterTime[[iCluster]]
                    return(solve(init.Omega[iTime,iTime,drop=FALSE]))
                }), Upattern$name)
            }else{
                start.OmegaM1 <- stats::setNames(lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1
                    diag(1, nrow = Upattern[iPattern,"n.time"], ncol = Upattern[iPattern,"n.time"])
                }), Upattern$name)
            }
            param.value[param.mu2] <- .optimGLS(OmegaM1 = start.OmegaM1, pattern = Upattern$name, precompute.XY = precompute.XY, precompute.XX = precompute.XX, key.XX = key.XX,
                                                Y = partialY, design = design, param.mu = param.mu2)
        }
        ## vcov values
        iResiduals.long <- partialY - design$mean[,param.mu2,drop=FALSE] %*% param.value[param.mu2]
        if(length(param.Omega2)>0){
            if(is.null(init.Omega)){
                outInit <- .initialize(design$vcov, init.cor = init.cor, method.fit = method.fit, residuals = iResiduals.long, Xmean = design$mean, index.cluster = index.cluster)
            }else{
                outInit <- .initialize2(design$vcov, index.clusterTime = design$index.clusterTime, Omega = init.Omega)
            }
        }else{
            outInit <- NULL
        }

        ## check initialization leads to a positive definite matrix 
        initOmega <- .calc_Omega(object = design$vcov, param = outInit, simplify = FALSE)        
        test.npd <- sapply(initOmega,function(iOmega){any(eigen(iOmega, symmetric = TRUE)$values<0)})
        if(any(test.npd)){ ## otherwise initialize as compound symmetry
            param.value[setdiff(param.sigma,param.fixed)] <- outInit[setdiff(param.sigma,param.fixed)]
            param.value[setdiff(param.k,param.fixed)] <- outInit[setdiff(param.k,param.fixed)]
            param.value[setdiff(param.rho,param.fixed)] <- stats::median(outInit[setdiff(param.rho,param.fixed)])            
        }else{        
            param.value[names(outInit)] <- outInit
        }
    }else{
        param.value <- init[param.name]     
    }
    if(trace>1){
        cat("Initialization:\n")
        print(param.value)
    }

    ## ** loop
    update.param.Omega <- stats::setNames(rep(0, length(param.Omega)), param.Omega)

    if(n.iter==0 || length(param.Omega2)==0){
        cv <- as.numeric(length(param.Omega2)==0)
        param.valueM1 <- NULL
        logLik.value <- .moments.lmm(value = param.value, design = design, time = time, method.fit = method.fit, type.information = type.information,
                                     transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                     logLik = TRUE, score = FALSE, information = FALSE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = effects, robust = FALSE,
                                     trace = FALSE, precompute.moments = precompute.moments, transform.names = FALSE)$logLik
        logLik.valueM1 <- NULL
        score.value <- NULL
        iIter <- 0
    }else if(optimizer=="FS"){
        cv <- 0
        param.valueM1 <- NULL
        logLik.value <- -Inf
        score.value <- stats::setNames(rep(Inf, length(param.value)), names(param.value))
        information.value <- NULL
        type.information <- "expected"
        ## wolfe.c1 <- 1e-4
        ## wolfe.c2 <- 0.9
        iIter <- 0
        if(trace>1){
            cat("\nLoop:\n")
        }
        for(iiIter in 0:(n.iter-1)){ ## iIter <- 1
            logLik.valueM1 <- logLik.value
            score.valueM1 <- score.value
            information.valueM1 <- information.value

            ## *** estimate moments
            outMoments <- .moments.lmm(value = param.value, design = design, time = time, method.fit = method.fit, type.information = type.information,
                                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                       logLik = TRUE, score = TRUE, information = TRUE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = effects, robust = FALSE,
                                       trace = FALSE, precompute.moments = precompute.moments, transform.names = FALSE)

            logLik.value <- outMoments$logLik    
            score.value <- outMoments$score    
            information.value <- outMoments$information

            ## if(iiIter==0){
            ##     test.wolfe <- c(TRUE,TRUE)
            ## }else{
            ##     ## wolfe condition for full step
            ##     test.wolfe <- .wolfe(update = update.value,
            ##                          logLik.old = logLik.valueM1, logLik.new = logLik.value, 
            ##                          score.old = score.valueM1, score.new = score.value, alpha = 1, c1 = wolfe.c1, c2 = wolfe.c2)
            ## }
            if(length(param.Omega2)==0){ ## deal with special case of no variance coefficient (e.g. when performing profile likelihood on sigma)
                if(iIter==0){
                    ## continue (do one iteration to get GLS estimate)
                }else{
                    cv <- 1
                    break
                }
            }else if(all(!is.na(outMoments$score)) && all(abs(outMoments$score)<tol.score) && (iiIter==0 || all(abs(param.valueM1 - param.value)<tol.param))){
                if(iiIter==0){
                    param.valueM1 <- param.value * NA
                }
                cv <- 1
                break
            }else if(iiIter == 0 && is.na(logLik.value)){
                cv <- -2
                break
            }else if(is.na(logLik.value) || (logLik.value < logLik.valueM1)){ ## decrease in likelihood - try partial update
                outMoments <- .backtracking(valueM1 = param.valueM1, update = update.param.Omega, n.iter = n.backtracking,
                                            design = design, time = time, method.fit = method.fit, type.information = type.information,
                                            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                            logLikM1 = logLik.valueM1, scoreM1 = score.valueM1, informationM1 = information.valueM1, effects = effects, precompute.moments = precompute.moments,
                                            precompute.XY = precompute.XY, precompute.XX = precompute.XX, key.XX = key.XX, Y = partialY, param.mu = param.mu2, param.Omega = param.Omega2)
                
                if(attr(outMoments,"cv")==FALSE){
                    cv <- -1
                    param.value <- param.valueM1 ## revert back to previous iteration
                    break
                }else{
                    param.value <- attr(outMoments,"value")
                    logLik.value <- outMoments$logLik    
                    score.value <- outMoments$score    
                    information.value <- outMoments$information
                }
            }

            ## *** update variance-covariance estimate
            param.valueM1 <- param.value
            if(length(param.Omega2)>0){
                
                ## update variance-covariance parameters (transform scale)
                update.param.Omega[param.Omega2] <- stats::setNames(as.double(score.value[param.Omega2] %*% solve(information.value[param.Omega2,param.Omega2,drop=FALSE])), param.Omega2)
                param.newvalue.trans <- outMoments$reparametrize$p + update.param.Omega
                ## back to original (transform scale)
                param.value[param.Omega] <- .reparametrize(param.newvalue.trans,
                                                           type = design.param2$type,
                                                           sigma = design.param2$sigma,
                                                           k.x = design.param2$k.x,
                                                           k.y = design.param2$k.y,
                                                           Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                                           transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                                           transform.names = FALSE)$p
            }
            ## *** update mean estimate
            if(length(param.mu2)>0){
                iOmega <- .calc_Omega(object = design$vcov, param = param.value, simplify = FALSE)
                param.value[param.mu2] <- .optimGLS(OmegaM1 = stats::setNames(lapply(iOmega, solve), names(iOmega)),
                                                    pattern = Upattern$name, precompute.XY = precompute.XY, precompute.XX = precompute.XX, key.XX = key.XX,
                                                    Y = partialY, design = design,
                                                    param.mu = param.mu2)
            }
            
            ## *** display
            iIter <- iIter+1
            if(trace > 0 && trace < 3){
                cat("*")
            }else{
                if(!is.null(attr(outMoments,"n.backtrack"))){
                    txt.backtract <- paste0(" (",attr(outMoments,"n.backtrack")," backtrack)")
                }else{
                    txt.backtract <- ""
                }

                if(trace==3){
                    cat("iteration ",iIter,txt.backtract,": logLik=",formatC(outMoments$logLik, digits = 10),"\n",sep="")
                }else if(trace==4){
                    cat("iteration ",iIter,txt.backtract,": logLik=",formatC(outMoments$logLik, digits = 10),"\n",sep="")
                    print(param.value)
                }else if(trace > 4){
                    cat("iteration ",iIter,txt.backtract,": logLik=",formatC(outMoments$logLik, digits = 10),"\n",sep="")
                    M.print <- rbind(estimate = param.value,
                                     diff = c(param.value - param.valueM1),
                                     score = c(rep(NA, length(param.mu)),outMoments$score))
                    print(M.print)
                    cat("\n")
                    
                }
            }
            
        }

        if(cv>0){
            attr(cv,"message") <- "Convergence"
        }else if(cv==0){
            attr(cv,"message") <- "Stop optimization before convergence (maximum number of iterations reached)"
        }else if(cv==-1){
            attr(cv,"message") <- "Stop optimization before convergence (decreasing log-likelihood)"
        }else if(cv==-2){
            attr(cv,"message") <- "Stop optimization before convergence (log-likelihood=NA based on the initial values)"
        }
        if(trace>=1){
            if(trace %in% 2:3){
                cat("\n")
                print(param.value)
            }
            if(length(param.Omega2)==0){
                cat(attr(cv,"message")," after ",iIter," iteration. \n",sep="") ## only one iteration (GLS)
            }else if(cv==-2){
                cat(attr(cv,"message"),". \n",sep="") ## incorrect initialization
            }else if(iIter==0){
                cat(attr(cv,"message")," after ",iIter," iteration: max score=",max(abs(outMoments$score[param.Omega2])),"\n", sep = "")
            }else{
                cat(attr(cv,"message")," after ",iIter," iterations: max score=",max(abs(outMoments$score[param.Omega2]))," | max change in coefficient=",max(abs(param.valueM1 - param.value)),"\n", sep = "")
            }
        }
        score <- outMoments$score
    }else{
        warper_obj <- function(p){
            p.original <- .reparametrize(p,
                                         type = design$param[match(names(param.value), design$param$name), "type"],
                                         sigma = design$param[match(names(param.value), design$param$name), "sigma"],
                                         k.x = design$param[match(names(param.value), design$param$name), "k.x"],
                                         k.y = design$param[match(names(param.value), design$param$name), "k.y"],
                                         Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                         transform.names = FALSE)$p
            -.moments.lmm(value = p.original, design = design, time = time, method.fit = method.fit, type.information = "observed",
                          transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                          logLik = TRUE, score = FALSE, information = FALSE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = c("mean","variance","correlation"), robust = FALSE,
                          trace = FALSE, precompute.moments = precompute.moments, transform.names = FALSE)$logLik
        }
        warper_grad <- function(p){
            p.original <- .reparametrize(p,
                                         type = design$param[match(names(param.value), design$param$name), "type"],
                                         sigma = design$param[match(names(param.value), design$param$name), "sigma"],
                                         k.x = design$param[match(names(param.value), design$param$name), "k.x"],
                                         k.y = design$param[match(names(param.value), design$param$name), "k.y"],
                                         Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                         transform.names = FALSE)$p
            -.moments.lmm(value = p.original, design = design, time = time, method.fit = method.fit, type.information = "observed",
                          transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                          logLik = FALSE, score = TRUE, information = FALSE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = c("mean","variance","correlation"), robust = FALSE,
                          trace = FALSE, precompute.moments = precompute.moments, transform.names = FALSE)$score
        }
        warper_hess <- function(p){
            p.original <- .reparametrize(p,
                                         type = design$param[match(names(param.value), design$param$name), "type"],
                                         sigma = design$param[match(names(param.value), design$param$name), "sigma"],
                                         k.x = design$param[match(names(param.value), design$param$name), "k.x"],
                                         k.y = design$param[match(names(param.value), design$param$name), "k.y"],
                                         Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                         transform.names = FALSE)$p
            .moments.lmm(value = p.original, design = design, time = time, method.fit = method.fit, type.information = "observed",
                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                         logLik = FALSE, score = FALSE, information = TRUE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = c("mean","variance","correlation"), robust = FALSE,
                         trace = FALSE, precompute.moments = precompute.moments, transform.names = FALSE)$information
        }

        ## *** reparametrize (original -> unconstrain scale)
        param.value.trans <- .reparametrize(param.value,
                                            type = design$param[match(names(param.value), design$param$name), "type"],
                                            sigma = design$param[match(names(param.value), design$param$name), "sigma"],
                                            k.x = design$param[match(names(param.value), design$param$name), "k.x"],
                                            k.y = design$param[match(names(param.value), design$param$name), "k.y"],
                                            Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                            transform.names = FALSE)$p

        ## *** optimize
        ## warper_obj(param.value.trans)
        ## warper_obj(2*param.value.trans)
        ## numDeriv::jacobian(x = param.value.trans, func = warper_obj)-warper_grad(param.value.trans)
        if(trace<=0){trace <- 0}
        res.optim <- optimx::optimx(par = param.value.trans, fn = warper_obj, gr = warper_grad, hess = warper_hess,
                                    method = optimizer, itnmax = n.iter, control = list(trace = trace))
        ## solution <- stats::setNames(as.double(res.optim[1,1:length(param.value.trans)]), names(param.value.trans))
        ## warper_obj(solution)
        ## warper_grad(solution)

        ## *** reparametrize (unconstrain scale -> original)
        param.value[] <- .reparametrize(as.double(res.optim[1:length(param.value)]),
                                        type = design$param[match(names(param.value), design$param$name), "type"],
                                        sigma = design$param[match(names(param.value), design$param$name), "sigma"],
                                        k.x = design$param[match(names(param.value), design$param$name), "k.x"],
                                        k.y = design$param[match(names(param.value), design$param$name), "k.y"],
                                        Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                        transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                        transform.names = FALSE)$p

        param.valueM1 <- NULL
        score.value <- attr(res.optim,"details")[,"ngatend"][[1]]
        logLik.value <- NULL
        logLik.valueM1 <- NULL
        iIter <- res.optim$niter
        attr(iIter,"eval") <- c("logLik" = NA, "score" = NA)
        if("fevals" %in% names(res.optim)){
            attr(iIter,"eval")["logLik"] <- res.optim$fevals
        }
        if("gevals" %in% names(res.optim)){
            attr(iIter,"eval")["score"] <- res.optim$gevals
        }
        cv <- as.numeric(res.optim$convcode==0)
    }

    ## ** export
    return(list(estimate = param.value,
                previous.estimate = param.valueM1,
                logLik = logLik.value,
                previous.logLik = logLik.valueM1,
                score = score.value,
                n.iter = iIter,                
                cv = cv,
                control = c(n.iter = as.double(n.iter), tol.score = as.double(tol.score), tol.param = as.double(tol.param),
                            n.backtracking = as.double(n.backtracking), init.cor = as.double(init.cor))
                ))
}

## * .optimGLS
## Implement GLS estimator i.e. \beta = (tX \OmegaM1 X)^{1} tX \OmegaM1 Y
.optimGLS <- function(OmegaM1, pattern, precompute.XY, precompute.XX, key.XX, Y, design, param.mu){
    name.param <- param.mu
    n.param <- length(name.param)
    numerator <- matrix(0, nrow = n.param, ncol = 1)
    denominator <- matrix(0, nrow = n.param, ncol = n.param)

    if(!is.null(key.XX)){
        max.key <- key.XX[n.param,n.param]
    }

    for(iPattern in pattern){ ## iPattern <- pattern[1]
        if(!is.null(precompute.XX) && !is.null(precompute.XY)){
            iVec.Omega <- as.double(OmegaM1[[iPattern]])
            iTime2 <- length(iVec.Omega)
            numerator <- numerator + t(iVec.Omega %*%  precompute.XY[[iPattern]])
            denominator <- denominator + as.double(iVec.Omega %*%  precompute.XX[[iPattern]])[key.XX]
        }else{
            iIndexCluster <- design$index.cluster[design$vcov$pattern == which(pattern==iPattern)]
            for(iId in 1:length(iIndexCluster)){ ## iId <- 2
                iX <- design$mean[iIndexCluster[[iId]],param.mu,drop=FALSE]
                numerator  <- numerator + design$weights[iId] * (t(iX) %*% OmegaM1[[iPattern]] %*% Y[iIndexCluster[[iId]]])
                denominator  <- denominator + design$weights[iId] * (t(iX) %*% OmegaM1[[iPattern]] %*% iX)
            }

        }
    }

    out <- solve(denominator, numerator)    
    return(stats::setNames(as.double(out), name.param))
}

## * .wolfe
## Nocedal 2000 Numerical optimization page 34 (formula 3.
## (3.6a) f(x_k + \alpha_k p_k) <= f(x_k) + c_1 \alpha_k \nabla f_k p_k
##        -logLik_{k+1} <= -logLik_{k} + c_1 \alpha_k Score_k update_k
## MAKE SURE THAT THE CHANGE IN LOGLIK WORTH THE AMOUNT OF UPDATE OTHERWISE WE SHOULD HALF THE UPDATE
## (3.6b) \nable f(x_k + \alpha_k p_k) p_k >= c_2 \nabla f_k p_k
##        -Score_{k+1} update_k <=c2 Score_k update_k
## MAKE SURE THAT THE SCORE IS NOT TOO NEGATIVE OTHERWISE WE COULD GO FURTHER
.wolfe <- function(update, logLik.old, logLik.new, score.old, score.new, alpha, c1, c2){
    if(!is.na(c1)){
        test.a <- as.vector( (-logLik.new) <= (-logLik.old) + c1 * alpha * update %*% (-score.old) )
    }else{
        test.a <- TRUE
    }
    if(!is.na(c2)){
        test.b <- as.vector(update %*% (-score.new) <= - c2 * update %*% (-score.old))
    }else{
        test.b <- TRUE
    }
    return(c(test.a, test.b))
}

## * backtracking
.backtracking <- function(valueM1, update, n.iter,
                          design, time, method.fit, type.information,
                          transform.sigma, transform.k, transform.rho, 
                          logLikM1, scoreM1, informationM1, information, effects, precompute.moments,
                          precompute.XY, precompute.XX, key.XX, Y, param.mu, param.Omega){

    if(n.iter<=0){return(list(cv = FALSE))}
    
    alpha <- 1/2
    design.param <- design$param[match(param.Omega, design$param$name),,drop=FALSE]                                                
    valueNEW <- valueM1
    momentNEW <- NULL
    cv <- FALSE

    ## move param to transform scale
    valueM1.trans <- .reparametrize(valueM1[param.Omega],
                                    type = design.param$type,
                                    sigma = design.param$sigma,
                                    k.x = design.param$k.x,
                                    k.y = design.param$k.y,
                                    Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                    transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                    transform.names = FALSE)$p

    for(iIter in 1:n.iter){

        ## update variance-covariance parameters (transform scale)
        valueNEW.trans <- stats::setNames(valueM1.trans[param.Omega] + alpha * update[param.Omega], param.Omega)

        ## update variance-covariance parameters (original scale)
        valueNEW[param.Omega] <- .reparametrize(valueNEW.trans,
                                                type = design.param$type,
                                                sigma = design.param$sigma,
                                                k.x = design.param$k.x,
                                                k.y = design.param$k.y,
                                                Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                                transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                                transform.names = FALSE)$p

        ## estimate residual variance-covariance matrix
        iOmega <- .calc_Omega(object = design$vcov, param = valueNEW, simplify = FALSE)
        iDet <- sapply(iOmega,det) 
        if(any(is.na(iDet)) || any(iDet<0)){
            alpha <- alpha/2
            next
        }
        
        ## update mean parameters
        if(length(param.mu)>0){
            valueNEW[param.mu] <- .optimGLS(OmegaM1 = stats::setNames(lapply(iOmega, solve), names(iOmega)),
                                            pattern = design$vcov$Upattern$name,
                                            precompute.XY = precompute.XY,
                                            precompute.XX = precompute.XX,
                                            key.XX = key.XX,
                                            Y = Y,
                                            design = design,
                                            param.mu = param.mu)
        }

        ## compute moments
        momentNEW <- .moments.lmm(value = valueNEW, design = design, time = time, method.fit = method.fit, type.information = type.information,
                                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                  logLik = TRUE, score = TRUE, information = TRUE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = effects, robust = FALSE,
                                  trace = FALSE, precompute.moments = precompute.moments, transform.names = FALSE)

        ## ## check wolfe condition
        ## test.wolfe <- .wolfe(update,
        ##                      logLik.old = logLikM1, logLik.new = momentNEW$score,
        ##                      score.old = scoreM1, score.new = momentNEW$score,
        ##                      alpha = alpha, c1 = c1, c2 = c2)
        ## check convergence
        if(is.na(momentNEW$logLik) || momentNEW$logLik<logLikM1){
            alpha <- alpha/2
        }else{
            cv <- TRUE
            break
        }

    }

    ## ** export
    if(is.null(momentNEW)){
        momentNEW <- .moments.lmm(value = valueNEW, design = design, time = time, method.fit = method.fit, type.information = type.information,
                                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                  logLik = TRUE, score = TRUE, information = TRUE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = effects, robust = FALSE,
                                  trace = FALSE, precompute.moments = precompute.moments, transform.names = FALSE)
    }
    attr(momentNEW,"value") <- valueNEW
    attr(momentNEW,"cv") <- cv
    attr(momentNEW,"n.backtrack") <- iIter
    
    return(momentNEW)
}


##----------------------------------------------------------------------
### optim.R ends here
