### optim.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 17 2025 (11:27) 
## Version: 
## Last-Updated: okt 20 2026 (17:09) 
##           By: Brice Ozenne
##     Update #: 415
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

## * .optim.lmm (code)
.optim.lmm <- function(Y, design, n.time, method.fit, 
                       transform.sigma, transform.k, transform.rho, type.information = "expected",
                       optimizer, init, n.iter, tol.score, tol.param, n.backtracking, init.cor, trace){

    ## ** Set default option
    if(is.null(trace)){
        trace <- FALSE
    }

    ## ** Extract information
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

    n.param <- length(param.name)
    alpha.backtracking <- 1/2

    pattern2param <- design$structure$pattern2param
    index.clusterPattern <- design$index.clusterPattern

    ## ** Pre-compute
    index.cluster <- tapply(1:length(design$index.cluster),design$index.cluster,identity)
    precompute.weights <- .precomputeWeights(weights = design$weights, index.clusterPattern = index.clusterPattern, pattern2param = pattern2param)
    precompute.XXY <- .precomputeXXY(X = design$mean,
                                     Y = Y,
                                     index.cluster = index.cluster,
                                     index.clusterPattern = index.clusterPattern,
                                     pattern2param = pattern2param,
                                     weights = design$weights,
                                     precompute.weights = precompute.weights)
    pair.vcov <-  .precomputePairVcov(pattern2param = pattern2param,
                                      Xpattern.var = design$variance$Xpattern,
                                      Xpattern.cor = design$correlation$Xpattern,
                                      activeparam.vcov = param.Omega2,
                                      activeparam.all = c(param.mu2,param.Omega2))
    
    ## ** Estimate intialization values based on the empirical moments
    if(!is.null(init)){

        if(is.matrix(init)){
            if(NROW(init)!=n.time || NCOL(init)!=n.time){
                stop("When a matrix, initialization should be a square matrix with dimensions compatible with time (",n.time,"x",n.time,"). \n")
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

        ## Current residuals (possibly substract constant mean parameter effect)
        if(length(param.fixed.mu)>0){
            partialY <-  Y - design$mean[,param.fixed.mu,drop=FALSE] %*% design$param[!is.na(design$param$constraint),"constraint"]
            if(is.null(design$weights$multiplier)){
                pattern.ntime <- stats::setNames(pattern2param$n.time, pattern2param$name)
            }else{
                pattern.ntime <- stats::setNames(pattern2param$n.time, pattern2param$name)[precompute.weights$vcov]
            }
            pattern.cluster <- tapply(1:NROW(index.clusterPattern),factor(index.clusterPattern$all, levels = precompute.weights$pattern), FUN = identity)
            precompute.XpartialY <- .precomputeXR(X = design$mean, residuals = sweep(partialY, FUN = "*", MARGIN = 1, STATS = design$weights$XYR),
                                                  pattern = precompute.weights$pattern,
                                                  pattern.ntime = stats::setNames(pattern.ntime, precompute.weights$pattern),
                                                  pattern.cluster = pattern.cluster, index.cluster = index.cluster)
        }else{
            partialY <- Y
            precompute.XpartialY <- precompute.XXY$XY
        }

        ## initialize mean value
        if(!is.null(init.mu)){
            param.value[param.mu2] <- init.mu[param.mu2]
        }else if(length(param.mu2)>0){
            Upattern <- precompute.weights$pattern
            n.Upattern <- length(Upattern)
            if(!is.null(init.Omega)){
                start.OmegaM1 <- stats::setNames(lapply(1:n.Upattern, function(iP){ ## iP <- 1
                    iCluster <- which(index.clusterPattern == precompute.weights$vcov[iP])
                    iTime <- pattern2param[pattern2param$name == precompute.weights$vcov[iP],"n.time"]
                    return(solve(init.Omega[iTime,iTime,drop=FALSE]))
                }), Upattern)
            }else{
                start.OmegaM1 <- stats::setNames(lapply(1:n.Upattern, function(iP){ ## iP <- 1
                    iN.time <- pattern2param[pattern2param$name == precompute.weights$vcov[iP],"n.time"]
                    return(diag(1, nrow = iN.time, ncol = iN.time))
                }), Upattern)
            }
            param.value[param.mu2] <- .GLS(OmegaM1 = start.OmegaM1, precompute.XY = precompute.XpartialY, precompute.XX = precompute.XXY$XX, precompute.weights = precompute.weights)
        }                                                  
        ## initialize vcov values based on the residuals
        iResiduals.long <- partialY - design$mean[,param.mu2,drop=FALSE] %*% param.value[param.mu2]
        if(length(param.Omega2)>0){
            if(is.null(init.Omega)){
                browser()
                outInit <- .initialize(design$structure, init.cor = init.cor, method.fit = method.fit,
                                       residuals = iResiduals.long, Xmean = design$mean, index.cluster = design$index.cluster)
            }else{
                outInit <- .initialize2(design$structure,
                                        index.clusterTime = design$index.clusterTime, Omega = init.Omega)
            }
        }else{
            outInit <- NULL
        }
        
        ## check initialization leads to a positive definite matrix 
        initOmega <- .calc_Omega(object = design$vcov, param = outInit, 
                                 transform.sigma = "none", transform.k = "none", transform.rho = "none", simplify = FALSE)        
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

    ## ** Apply transformation to initialization values
    reparametrize <- .reparametrize(param.value,
                                    type = design$param$type,
                                    level = design$param$level,
                                    sigma = design$param$sigma,
                                    k.x = design$param$k.x,
                                    k.y = design$param$k.y,
                                    Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                    transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                    transform.names = TRUE)
    paramTrans.value <- reparametrize$p

    ## ** Gradient descent
    if(n.iter<1){ ## no iteration: n.iter = 0 (only return logLikelihood) or ]0;1[ also checks convergence
        logLik.valueOLD <- NULL
        paramTrans.valueOLD <- NULL
        
        outMoments <- .moments.lmm(value = paramTrans.value, Y = Y, design = design, method.fit = method.fit, type.information = type.information,
                                   precompute.weights = precompute.weights, precompute.XX = precompute.XXY$XX, precompute.XY = precompute.XXY$XY, pair.vcov = pair.vcov,
                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                   updateGLS = FALSE, logLik = TRUE, score = n.iter>0, information = FALSE, vcov = FALSE, df = FALSE, indiv = FALSE,
                                   effects = c("variance","correlation"), robust = FALSE,
                                   trace = FALSE, transform.names = FALSE, export.precompute = FALSE)

        logLik.value <- outMoments$logLik
        update <- 0
        if(n.iter>0){
            score.value <- outMoments$score
            cv <- all(!is.na(score.value)) && all(abs(score.value)<tol.score)
        }else{
            score.value <- NULL
            cv <- NA
        }
    }else if(length(param.Omega2)==0){ ## no variance-covariance paramater: GLS estimator
        logLik.valueOLD <- NULL
        paramTrans.valueOLD <- NULL

        outMoments <- .moments.lmm(value = paramTrans.value, Y = Y, design = design, method.fit = method.fit, type.information = type.information,
                                   precompute.weights = precompute.weights, precompute.XX = precompute.XXY$XX, precompute.XY = precompute.XXY$XY, pair.vcov = pair.vcov,
                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                   updateGLS = TRUE, logLik = TRUE, score = TRUE, information = FALSE, vcov = FALSE, df = FALSE, indiv = FALSE,
                                   effects = c("variance","correlation"), robust = FALSE,
                                   trace = FALSE, transform.names = FALSE, export.precompute = FALSE)

        paramTrans.value <- outMoments$value ## GLS update
        logLik.value <- outMoments$logLik
        score.value <- outMoments$score
        cv <- all(!is.na(score.value)) && all(abs(score.value)<tol.score)
        update <- 0        
    }else if(optimizer=="FS"){
        cv <- 0
        update <- 0
        backtrack <- 0
        paramTrans.valueOLD <- paramTrans.value
        logLik.value <- -Inf
        momentNEW <- NULL

        ## wolfe.c1 <- 1e-4
        ## wolfe.c2 <- 0.9
        if(trace>1){
            cat("\nLoop:\n")
        }

        for(iIter in 0:(n.iter-1)){ ## iIter <- 1
            ## *** is the current estimate good enough?
            ## vcov is not needed here but is needed at the next step to obtain the update
            moment <- .moments.lmm(value = paramTrans.value, Y = Y, design = design, method.fit = method.fit, type.information = type.information,
                                   precompute.weights = precompute.weights, precompute.XX = precompute.XXY$XX, precompute.XY = precompute.XXY$XY, pair.vcov = pair.vcov,
                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                   updateGLS = FALSE, logLik = (iIter==0), score = TRUE, information = FALSE, vcov = TRUE, df = FALSE, indiv = FALSE,
                                   effects = c("variance","correlation"), robust = FALSE,
                                   trace = FALSE, transform.names = FALSE)
            score.value <- moment$score

            if(trace > 0 && trace < 3){
                cat("*")
            }else{
                if(iIter>0 && backtrack>0){
                    txt.backtract <- paste0(" (",iBacktrack," backtrack)")
                }else{
                    txt.backtract <- ""
                }
                
                if(trace==3){
                    cat("iteration ",iIter,txt.backtract,": logLik=",formatC(moment$logLik, digits = 10),"\n",sep="")
                }else if(trace==4){
                    cat("iteration ",iIter,txt.backtract,": logLik=",formatC(moment$logLik, digits = 10),"\n",sep="")
                    print(paramTrans.value)
                }else if(trace > 4){
                    cat("iteration ",iIter,txt.backtract,": logLik=",formatC(moment$logLik, digits = 10),"\n",sep="")
                    if(iIter == 0){
                        M.print <- rbind(estimate = paramTrans.value,
                                         diff = NA,
                                         score = c(rep(NA, length(param.mu)),score.value))
                    }else{
                        M.print <- rbind(estimate = paramTrans.value,
                                         diff = c(paramTrans.value - paramTrans.valueM1),
                                         score = c(rep(NA, length(param.mu)),score.value))
                    }
                    print(M.print)
                    cat("\n")
                    
                }
            }

            if(iIter == 0 && is.na(logLik.value)){ #### initialization failure ####                
                logLik.value <- moment$logLik
                paramTrans.valueOLD <- NULL
                cv <- -2
                break
            }else if(iIter == 0 && all(!is.na(score.value)) && all(abs(score.value)<tol.score)){ #### succes: 0-score initialization ####
                logLik.value <- moment$logLik
                paramTrans.valueOLD <- NULL
                cv <- 1
                break
            }else if(all(!is.na(score.value)) && all(abs(score.value)<tol.score) && all(abs(paramTrans.value - paramTrans.valueOLD)<tol.param)){ #### success: 0-score and stable estimate ####
                cv <- 1
                break
            }else{ #### estimate should be improved ####
                
                ## *** possible update (fisher scoring)
                update.param.Omega <- as.double(score.value[param.Omega2] %*% moment$vcov[param.Omega2,param.Omega2,drop=FALSE])

                ## *** is the update acceptable? Otherwise backtrack
                valueNEW <- paramTrans.value
                cv.backtrack <- FALSE
                
                for(iBacktrack in 0:n.backtracking){ ## iBacktrack <- 0
                    valueNEW[param.Omega] <- paramTrans.value[param.Omega] + alpha.backtracking/2^iBacktrack * update.param.Omega

                    ## GLS estimator for mean parameters + log-likelihood
                    momentNEW <- .moments.lmm(value = valueNEW, Y = Y, design = design, Omega = NULL, method.fit = method.fit, type.information = type.information,
                                              precompute.weights = precompute.weights, precompute.XX = precompute.XXY$XX, precompute.XY = precompute.XXY$XY, pair.vcov = pair.vcov,
                                              transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                              updateGLS = TRUE, logLik = TRUE, score = FALSE, information = FALSE, vcov = FALSE, df = FALSE, indiv = FALSE,
                                              effects = c("variance","correlation"), robust = FALSE,
                                              trace = FALSE, transform.names = FALSE)

                    if(!is.na(momentNEW$logLik) && (momentNEW$logLik > logLik.value)){
                        cv.backtrack <- TRUE
                        backtrack <- iBacktrack
                        break
                    }
                }
                
                if(cv.backtrack){ ## update parameter and log-likelihood
                    update <- update+1 
                    paramTrans.valueOLD <- paramTrans.value
                    paramTrans.value <- momentNEW$value
                    logLik.valueOLD <- logLik.value
                    logLik.value <- momentNEW$logLik

                    ## ## check wolfe condition
                    ## test.wolfe <- .wolfe(update,
                    ##                      logLik.old = logLikM1, logLik.new = momentNEW$score,
                    ##                      score.old = scoreM1, score.new = momentNEW$score,
                    ##                      alpha = alpha, c1 = c1, c2 = c2)
                }else{  ## back-tracking failed
                    if(iIter == 0){
                        paramTrans.valueOLD <- NULL
                    }
                    backtrack <- n.backtracking
                    cv <- -1
                    break
                }
            }
        }
        
        if(cv>0){
            attr(cv,"message") <- "Convergence"
        }else if(cv==0){
            attr(cv,"message") <- "Stop optimization before convergence (maximum number of iterations reached)"
        }else if(cv==-1){
            if(is.na(logLik.value)){
                attr(cv,"message") <- "Stop optimization before convergence (log-likelihood = NA)"
            }else{
                attr(cv,"message") <- "Stop optimization before convergence (decreasing log-likelihood)"
            }
        }else if(cv==-2){
            attr(cv,"message") <- "Stop optimization before convergence (log-likelihood=NA based on the initial values)"
        }
        if(trace>=1){
            if(trace %in% 2:3){
                cat("\n")
                print(paramTrans.value)
            }
            if(length(param.Omega2)==0 || cv==-2){
                cat(attr(cv,"message"),". \n",sep="")
            }else if(update==0){
                cat(attr(cv,"message")," after ",iIter," iteration(s): max score=",max(abs(score.value)),"\n", sep = "")
            }else{
                cat(attr(cv,"message")," after ",iIter," iteration(s): max score=",max(abs(score.value))," | max change in coefficient=",max(abs(paramTrans.valueOLD - paramTrans.value)),"\n", sep = "")
            }
        }
    }else{
        warper_obj <- function(p){
            -.moments.lmm(value = p, Y = Y, design = design, method.fit = method.fit, type.information = type.information,
                          precompute.weights = precompute.weights, precompute.XX = precompute.XXY$XX, precompute.XY = precompute.XXY$XY, pair.vcov = pair.vcov,
                          transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                          updateGLS = TRUE, logLik = TRUE, score = FALSE, information = FALSE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = c("variance","correlation"), robust = FALSE,
                          trace = FALSE, transform.names = FALSE)$logLik
        }
        warper_grad <- function(p){
            -.moments.lmm(value = p, Y = Y, design = design, method.fit = method.fit, type.information = type.information,
                          precompute.weights = precompute.weights, precompute.XX = precompute.XXY$XX, precompute.XY = precompute.XXY$XY, pair.vcov = pair.vcov,
                          transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                          updateGLS = TRUE, logLik = FALSE, score = TRUE, information = FALSE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = c("variance","correlation"), robust = FALSE,
                          trace = FALSE, transform.names = FALSE)$score
        }
        warper_hess <- function(p){
            .moments.lmm(value = p, Y = Y, design = design, method.fit = method.fit, type.information = type.information,
                         precompute.weights = precompute.weights, precompute.XX = precompute.XXY$XX, precompute.XY = precompute.XXY$XY, pair.vcov = pair.vcov,
                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                         updateGLS = TRUE, logLik = FALSE, score = FALSE, information = TRUE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = c("variance","correlation"), robust = FALSE,
                         trace = FALSE, transform.names = FALSE)$information
        }

        ## *** optimize
        ## warper_obj(paramTrans.value)
        ## warper_obj(2*paramTrans.value)
        ## numDeriv::jacobian(x = paramTrans.value, func = warper_obj)-warper_grad(paramTrans.value)
        if(trace<=0){trace <- 0}
        res.optim <- optimx::optimx(par = paramTrans.value[param.Omega2], fn = warper_obj, gr = warper_grad, hess = warper_hess,
                                    method = optimizer, itnmax = n.iter, control = list(trace = trace))
        ## solution <- stats::setNames(as.double(res.optim[1,1:length(paramTrans.value)]), names(paramTrans.value))
        ## warper_obj(solution)
        ## warper_grad(solution)

        paramTrans.valueM1 <- NULL
        score.value <- attr(res.optim,"details")[,"ngatend"][[1]]
        logLik.value <- NULL
        logLik.valueOLD <- NULL
        update <- res.optim$niter
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
    param.value <- .reparametrize(paramTrans.value,
                                  type = design$param$type,
                                  sigma = design$param$sigma,
                                  k.x = design$param$k.x,
                                  k.y = design$param$k.y,
                                  Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                  transform.names = FALSE)$p

    return(list(estimate.originalScale = param.value,
                estimate = stats::setNames(paramTrans.value, reparametrize$newname),
                previous.estimate = paramTrans.valueOLD,
                logLik = logLik.value,
                previous.logLik = logLik.valueOLD,
                score = score.value,
                n.iter = update,                
                cv = cv,
                precompute.weights = precompute.weights,
                precompute.XX = precompute.XXY$XX,
                precompute.XY = precompute.XXY$XY,
                pair.vcov = pair.vcov))
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

##----------------------------------------------------------------------
### optim.R ends here
