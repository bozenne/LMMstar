### estimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 20 2021 (23:25) 
## Version: 
## Last-Updated: jun  1 2022 (14:51) 
##           By: Brice Ozenne
##     Update #: 686
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * estimate.lmm
##' @title Delta Method for Mixed Models
##' @description Perform a first order delta method
##'
##' @param x  a \code{lmm} object.
##' @param f [function] function of the model coefficient computing the parameter(s) of interest. Can accept extra-arguments.
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors. 
##' @param df [logical] Should degree of freedom, computed using Satterthwaite approximation, for the parameter of interest be output.
##' @param type.information [character] Should the expected information be used  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param ... extra arguments passed to \code{f}.
##'
##' @examples
##' if(require(lava)){
##' 
##' #### Random effect ####
##' set.seed(10)
##' dL <- sampleRem(1e2, n.times = 3, format = "long")
##' e.lmm1 <- lmm(Y ~ X1+X2+X3, repetition = ~visit|id, structure = "CS", data = dL)
##' coef(e.lmm1, effects = "ranef")
##' e.ranef <- estimate(e.lmm1, f  = function(p){coef(e.lmm1, p = p, effects = "ranef")})
##' e.ranef
##'
##' if(require(ggplot2)){
##' df.gg <- cbind(index = 1:NROW(e.ranef), e.ranef)
##' gg.ranef <- ggplot(df.gg, aes(x = index, y=estimate, ymin=lower, ymax = upper))
##' gg.ranef + geom_point() + geom_errorbar() + ylab("estimated random effect") + xlab("id")
##' }
##' 
##' #### ANCOVA via mixed model ####
##' set.seed(10)
##' d <- sampleRem(1e2, n.time = 2)
##' e.ANCOVA1 <- lm(Y2~Y1+X1, data = d)
##'
##' if(require(reshape2)){
##'    dL2 <- melt(d, id.vars = c("id","Y1","X1"),  measure.vars = c("Y1","Y2"))
##'    e.lmm <- lmm(value ~ variable + variable:X1, data = dL2, repetition = ~variable|id)
##' 
##'    e.delta <- estimate(e.lmm, function(p){
##'        c(Y1 = p["rho(Y1,Y2)"]*p["k.Y2"],
##'          X1 = p["variableY2:X1"]-p["k.Y2"]*p["rho(Y1,Y2)"]*p["variableY1:X1"])
##' })
##'    ## same estimate and similar standard errors. 
##'    e.delta
##'    summary(e.ANCOVA1)$coef
##'    ## Degrees of freedom are a bit off though
##' }
##'
##' }
##' @export
estimate.lmm <- function(x, f, df = TRUE, robust = FALSE, type.information = NULL, level = 0.95,
                         transform.sigma = "none", transform.k = "none", transform.rho = "none", ...){

    
    ## estimate
    beta <- coef(x, effects = "all", transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)

    ## partial derivative
    f.formals <- names(formals(f))
    if(length(f.formals)==1){
        fbeta <- f(beta)
        grad <- numDeriv::jacobian(func = f, x = beta)
    }else{
        fbeta <- f(beta, ...)
        grad <- numDeriv::jacobian(func = f, x = beta, ...)
    }

    ## extract variance-covariance
    Sigma <- vcov(x, df = 2*(df>0), effects = "all", robust = robust, type.information = type.information, ## 2*df is needed to return dVcov
                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)

    ## ** delta-method
    C.Sigma.C <- grad %*% Sigma %*% t(grad)
    C.sigma.C <- sqrt(diag(C.Sigma.C))

    ## second order?
    ## g(\thetahat) = g(\theta) + (\thetahat-\theta)grad + 0.5(\thetahat-\theta)lap(\thetahat-\theta) + ...
    ## Var[g(\thetahat)] = grad\Var[\thetahat-\theta]grad + 0.25\Var[(\thetahat-\theta)lap(\thetahat-\theta)] + \Cov((\thetahat-\theta)grad,(\thetahat-\theta)lap(\thetahat-\theta)) + ...
    ## https://stats.stackexchange.com/questions/427332/variance-of-quadratic-form-for-multivariate-normal-distribution
    ## Var[g(\thetahat)] = grad\Var[\thetahat-\theta]grad + 0.25*2*tr((lap\Var[\thetahat-\theta])^2) + 0 + ...
    
    ## lap <- numDeriv::hessian(func = f, x = beta) ## laplacian
    ## 2 * sum(diag(Sigma %*% lap %*% Sigma %*% lap))
    
    ## df 
    if(!is.null(attr(Sigma, "dVcov"))){
        keep.param <- dimnames(attr(Sigma, "dVcov"))[[3]]
        C.dVcov.C <- sapply(keep.param, function(iM){ ##  iName  <- dimnames(attr(Sigma, "dVcov"))[[3]][1]
            rowSums(grad %*% attr(Sigma, "dVcov")[,,iM] * grad)
        })
        numerator <- 2 *diag(C.Sigma.C)^2
        denom <- rowSums(C.dVcov.C %*% Sigma[keep.param,keep.param,drop=FALSE] * C.dVcov.C)
        df <- numerator/denom
    }else{
        df <- rep(Inf, NROW(grad))
    }

    ## ** export
    alpha <- 1-level
    out <- data.frame(estimate = as.double(fbeta),
                      se = as.double(C.sigma.C),
                      df = as.double(df),
                      lower = as.double(fbeta + stats::qt(alpha/2, df = df) * C.sigma.C),
                      upper = as.double(fbeta + stats::qt(1-alpha/2, df = df) * C.sigma.C),
                      p.value = as.double(2*(1-stats::pt(abs(fbeta/C.sigma.C), df = df))))
    attr(out,"gradient") <- grad
    if(!is.null(names(fbeta))){
        rownames(out) <- names(fbeta)
        rownames(attr(out,"gradient")) <- names(fbeta)
    }
    colnames(attr(out,"gradient")) <- names(beta)
    return(out)
}


## * .estimate (documentation)
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
##' @param trace [1, 2, or 3] should each iteration be displayed?
##' 
##' @examples
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
##' @keywords internal

## * .estimate (code)
.estimate <- function(design, time, method.fit, type.information,
                      transform.sigma, transform.k, transform.rho,
                      precompute.moments, optimizer, init, n.iter, tol.score, tol.param, trace){

   
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
    index.cluster <- design$index.cluster
    
    param.mu <- design$param[design$param$type=="mu","name"]
    param.sigma <- design$param[design$param$type=="sigma","name"]
    param.k <- design$param[design$param$type=="k","name"]
    param.rho <- design$param[design$param$type=="rho","name"]
    param.Omega <- c(param.sigma,param.k,param.rho)
    param.name <- c(param.mu,param.Omega)
    n.param <- length(param.name)
    Upattern <- design$vcov$X$Upattern
    n.Upattern <- NROW(Upattern)
    precompute.XY <- design$precompute.XY
    precompute.XX <- design$precompute.XX$pattern
    key.XX <- design$precompute.XX$key
    effects <- c("variance","correlation")

    ## ** intialization
    if(is.null(init)){

        param.value <- stats::setNames(rep(NA, n.param),param.name)
        if(trace>1){
            cat("\nInitialization:\n")
        }

        ## mean value
        start.OmegaM1 <- stats::setNames(lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1
            diag(1, nrow = Upattern[iPattern,"n.time"], ncol = Upattern[iPattern,"n.time"])
        }), Upattern$name)
        param.value[param.mu] <- .estimateGLS(OmegaM1 = start.OmegaM1, pattern = Upattern$name, precompute.XY = precompute.XY, precompute.XX = precompute.XX, key.XX = key.XX,
                                              design = design)

        ## vcov values
        iResiduals.long <- design$Y - design$mean %*% param.value[param.mu]
        outInit <- .initialize(design$vcov, residuals = iResiduals.long, Xmean = design$mean, index.cluster = index.cluster)
        param.value[names(outInit)] <- outInit
        if(trace>1){
            print(param.value)
        }
    }else{
        if(any(param.name %in% names(init) == FALSE)){
            stop("Initialization does not contain value for all parameters. \n",
                 "Missing parameters: \"",paste(param.name[param.name %in% names(init) == FALSE], collapse = "\" \""),"\". \n")
        }
        param.value <- init[param.name]

        if(trace>1){
            cat("\nInitialization:\n")
            print(param.value)
        }
    }
    ## microbenchmark(test =     .estimateGLS(OmegaM1 = start.OmegaM1, pattern = pattern, precompute.XY = precompute.XY, precompute.XX = precompute.XX, key.XX = key.XX),
    ##                lm.fit = lm.fit(y = ncgsL$cholest[!is.na(ncgsL$cholest)], x = model.matrix(~time+highdose.time, data=ncgsL)[!is.na(ncgsL$cholest),]),
    ##                lm = lm(cholest~time+highdose.time, data=ncgsL)
    ##                )


    
    ## ** loop
    if(n.iter==0){
        cv <- FALSE
        param.valueM1 <- NULL
        score <- NULL
        iIter <- 0
    }else if(optimizer=="FS"){
        cv <- FALSE
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

            if(all(abs(outMoments$score)<tol.score) && (iiIter==0 || all(abs(param.valueM1 - param.value)<tol.param))){
                if(iiIter==0){param.valueM1 <- param.value * NA}
                cv <- TRUE
                break
            }else if(iiIter == 0 && is.na(logLik.value)){
                cv <- -1
                break
            }else if(is.na(logLik.value) || (logLik.value < logLik.valueM1)){ ## decrease in likelihood - try partial update
                outMoments <- .backtracking(valueM1 = param.valueM1, update = update.value, n.iter = options$param.optimizer["n.backtracking"],
                                            design = design, time = time, method.fit = method.fit, type.information = type.information,
                                            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                            logLikM1 = logLik.valueM1, scoreM1 = score.valueM1, informationM1 = information.valueM1, effects = effects, precompute.moments = precompute.moments,
                                            precompute.XY = precompute.XY, precompute.XX = precompute.XX, key.XX = key.XX, param.mu = param.mu)
                
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
            ## update variance-covariance parameters (transform scale)
            param.valueM1 <- param.value
            update.value <- stats::setNames(as.double(score.value %*% solve(information.value)), names(outMoments$reparametrize$p))
            param.newvalue.trans <- outMoments$reparametrize$p + update.value
            ## back to original (transform scale)
            param.value[names(param.newvalue.trans)] <- .reparametrize(param.newvalue.trans,
                                                                       type = design$param[match(names(param.newvalue.trans), design$param$name), "type"],
                                                                       sigma = design$param[match(names(param.newvalue.trans), design$param$name), "sigma"],
                                                                       k.x = design$param[match(names(param.newvalue.trans), design$param$name), "k.x"],
                                                                       k.y = design$param[match(names(param.newvalue.trans), design$param$name), "k.y"],
                                                                       Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                                                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                                                       transform.names = FALSE)$p
            
            ## *** update mean estimate
            iOmega <- .calc_Omega(object = design$vcov, param = param.value, keep.interim = TRUE)
            ## eigen(iOmega[[1]])
            param.value[param.mu] <- .estimateGLS(OmegaM1 = stats::setNames(lapply(iOmega, solve), names(iOmega)),
                                                  pattern = Upattern$name, precompute.XY = precompute.XY, precompute.XX = precompute.XX,
                                                  key.XX = key.XX,
                                                  design = design)

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

        if(trace==1){
            cat("\n")
        }else if(trace>1){
            if(trace %in% 2:3){
                cat("\n")
                print(param.value)
            }
            if(cv>0){
                if(iIter==0){
                    cat("Convergence after ",iIter," iteration: max score=",max(abs(outMoments$score)),"\n", sep = "")
                }else{
                    cat("Convergence after ",iIter," iterations: max score=",max(abs(outMoments$score))," | max change in coefficient=",max(abs(param.valueM1 - param.value)),"\n", sep = "")
                }
            }else if(cv==0){
                if(iIter==0){
                    cat("No convergence after ",iIter," iteration: max score=",max(abs(outMoments$score)),"\n")
                }else if(iIter==n.iter){
                    cat("No convergence after ",iIter," iterations: max score=",max(abs(outMoments$score))," | max change in coefficient= ",max(abs(param.valueM1 - param.value)),"\n", sep = "")
                }
            }else if(cv==-1){
                cat("Stop optimization after ",iIter," iterations as the log-likelihood was increasing. Max score=",max(abs(outMoments$score)),"\n", sep = "")
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
            -.moments.lmm(value = p.original, design = design, time = time, method.fit = method.fit, 
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
            -.moments.lmm(value = p.original, design = design, time = time, method.fit = method.fit, 
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
            .moments.lmm(value = p.original, design = design, time = time, method.fit = method.fit,
                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, type.information = "observed",
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
        res.optim <- optimx::optimx(par = param.value.trans, fn = warper_obj, gr = warper_grad, hess = warper_hess,
                                    method = optimizer, itnmax = n.iter, control = list(trace = trace))
        ## solution <- setNames(as.double(res.optim[1,1:length(param.value.trans)]), names(param.value.trans))
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
        score <- attr(res.optim,"details")[,"ngatend"][[1]]
        iIter <- res.optim$niter
        attr(iIter,"eval") <- c("logLik" = NA, "score" = NA)
        if("fevals" %in% names(res.optim)){
            attr(iIter,"eval")["logLik"] <- res.optim$fevals
        }
        if("gevals" %in% names(res.optim)){
            attr(iIter,"eval")["score"] <- res.optim$gevals
        }
        cv <- (res.optim$convcode==0)
    }

    ## ** export
    return(list(estimate = param.value,
                previous.estimate = param.valueM1,
                score = score,
                n.iter = iIter,                
                cv = cv>0))
}

## * .estimateGLS
## Implement GLS estimator i.e. \beta = (tX \OmegaM1 X)^{1} tX \OmegaM1 Y
.estimateGLS <- function(OmegaM1, pattern, precompute.XY, precompute.XX, key.XX, design){

    name.param <- design$param[design$param$type=="mu","name"] ## colnames(key.XX)
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
            numerator  <- numerator + t(iVec.Omega %*%  matrix(precompute.XY[[iPattern]], nrow = iTime2, ncol = n.param))
            denominator  <- denominator + as.double(iVec.Omega %*%  matrix(precompute.XX[[iPattern]], nrow = iTime2, ncol = max.key))[key.XX]
        }else{
            iOmegaM1 <- OmegaM1[[iPattern]]
            iIndexCluster <- design$index.cluster[which(design$vcov$X$pattern.cluster$pattern==iPattern)]
            for(iId in 1:length(iIndexCluster)){ ## iId <- 2
                iX <- design$mean[iIndexCluster[[iId]],,drop=FALSE]
                if(is.null(design$weight)){
                    iWeight <- 1
                }else{
                    iWeight <- design$weights[iId]
                }
                if(!is.null(design$scale.Omega)){
                    iWeight <- iWeight * design$scale.Omega[iId]
                }
                numerator  <- numerator + iWeight * (t(iX) %*% iOmegaM1 %*% design$Y[iIndexCluster[[iId]]])
                denominator  <- denominator + iWeight * (t(iX) %*% iOmegaM1 %*% iX)
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
                          precompute.XY, precompute.XX, key.XX, param.mu){

    if(n.iter<=0){return(list(cv = FALSE))}
    
    alpha <- 1/2
    name.update <- names(update)
    design.param <- design$param[match(name.update, design$param$name),,drop=FALSE]
    valueNEW <- stats::setNames(rep(NA, length(valueM1)), names(valueM1))
    cv <- FALSE

    ## move param to transform scale
    valueM1.trans <- .reparametrize(valueM1[name.update],
                                    type = design.param$type,
                                    sigma = design.param$sigma,
                                    k.x = design.param$k.x,
                                    k.y = design.param$k.y,
                                    Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                    transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                    transform.names = FALSE)$p

    for(iIter in 1:n.iter){

        ## update variance-covariance parameters (transform scale)
        valueNEW.trans <- valueM1.trans[name.update] + alpha * update

        ## update variance-covariance parameters (original scale)
        valueNEW[name.update] <- .reparametrize(valueNEW.trans,
                                                type = design.param$type,
                                                sigma = design.param$sigma,
                                                k.x = design.param$k.x,
                                                k.y = design.param$k.y,
                                                Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                                transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                                transform.names = FALSE)$p

        ## estimate residual variance-covariance matrix
        iOmega <- .calc_Omega(object = design$vcov, param = valueNEW, keep.interim = TRUE)
        if(any(sapply(iOmega,det)<0)){
            alpha <- alpha/2
            next
        }
        
        ## update mean parameters
        valueNEW[param.mu] <- .estimateGLS(OmegaM1 = stats::setNames(lapply(iOmega, solve), names(iOmega)),
                                           pattern = design$vcov$X$Upattern$name, precompute.XY = precompute.XY, precompute.XX = precompute.XX,
                                           key.XX = key.XX,
                                           design = design)

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
    attr(momentNEW,"value") <- valueNEW
    attr(momentNEW,"cv") <- cv
    attr(momentNEW,"n.backtrack") <- iIter
    return(momentNEW)
}
                
##----------------------------------------------------------------------
### estimate.R ends here
