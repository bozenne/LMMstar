### constrain.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 17 2022 (05:36) 
## Version: 
## Last-Updated: aug  3 2023 (15:41) 
##           By: Brice Ozenne
##     Update #: 86
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

.constrain.lmm <- function(x, effects, trace = FALSE, init = NULL, ...){

    ## ** normalize user input
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    name.effects <- names(effects)
    if(is.null(init)){
        init <- coef(x, effects = "all")
    }
    if(any(duplicated(name.effects))){
        stop("Incorrect argument \'effects\': contain duplicated names \"",paste(unique(name.effects[duplicated(name.effects)]), collapse = "\" \""),"\".\n")
    }
    if(any(name.effects %in% names(init) == FALSE)){
        stop("Incorrect argument \'effects\': unknown parameter(s) \"",paste(name.effects[name.effects %in% names(init) == FALSE], collapse = "\" \""),"\".\n")
    }
    init[name.effects] <- effects

    ## ** update model
    ## *** add constrain
    x$design$param[match(name.effects,x$design$param$name),"constraint"] <- effects
    ## no need to update the pairs of mean-vcov parameters as .estimate is only here to provide a point estimate

    ## *** update design matrix
    browser()
    if(meanconstraint){
    }

    ## *** update precomputation 
    browser()
    if(meanconstraint && x$args$precompute.moments && NCOL(x$design$mean)>0){
        if(is.na(var.weights[1])){
            wX.mean <- X.mean
            wY <- cbind(data[[var.outcome]])
        }else{
            wX.mean <- sweep(X.mean, FUN = "*", MARGIN = 1, STATS = sqrt(data[[var.weights[1]]]))
            wY <- cbind(data[[var.outcome]]*sqrt(data[[var.weights[1]]]))
        }

        precompute.wXX <-  .precomputeXX(X = wX.mean, pattern = structure$Upattern$name, 
                                        pattern.ntime = stats::setNames(structure$Upattern$n.time, structure$Upattern$name),
                                        pattern.cluster = attr(structure$pattern,"list"), index.cluster = index.cluster)

        precompute.wXY <-  .precomputeXR(X = precompute.wXX$Xpattern, residuals = wY, pattern = structure$Upattern$name,
                                        pattern.ntime = stats::setNames(structure$Upattern$n.time, structure$Upattern$name),
                                        pattern.cluster = attr(structure$pattern,"list"), index.cluster = index.cluster)
        ## update meanvcov pair - even if not useful since we use the expected information matrix
    }else if(vcovconstraint){
        ## update vcov pair, meanvcov pair
    }
    

    ## *** refit
    eee <- .estimate(design = x$design, time = x$time, method.fit = x$args$method.fit, 
                     transform.sigma = x$reparametrize$transform.sigma, transform.k = x$reparametrize$transform.k, transform.rho = x$reparametrize$transform.rho,
                     precompute.moments = "precompute.wXX" %in% names(x$design),
                     optimizer = "FS", init = init, n.iter = x$opt$control[["n.iter"]], tol.score = x$opt$control[["tol.score"]], tol.param = x$opt$control[["tol.param"]],
                     n.backtracking = x$opt$control[["n.backtracking"]], type.information = x$opt$control[["type.information"]],
                     options = options, trace = trace)

    ## ** export
    x$opt[c("cv","n.iter","score","previous.estimate")] <- eee[c("cv","n.iter","score","previous.estimate")]
    x$param <- eee$estimate
    x$logLik <- eee$logLik
    x[c("reparametrize","fitted","residuals","Omega","score","information","vcov","df","dVcov")] <- NULL

    class(x) <- append("clmm",class(x))
    return(x)
}

##----------------------------------------------------------------------
### constrain.R ends here
