### constrain.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 17 2022 (05:36) 
## Version: 
## Last-Updated: jul 21 2023 (15:56) 
##           By: Brice Ozenne
##     Update #: 71
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
    ## add constrain
    x$design$param[match(name.effects,x$design$param$name),"constraint"] <- effects
    ## refit
    eee <- .estimate(design = x$design, time = x$time, method.fit = x$args$method.fit, type.information = x$args$type.information,
                     transform.sigma = x$reparametrize$transform.sigma, transform.k = x$reparametrize$transform.k, transform.rho = x$reparametrize$transform.rho,
                     precompute.moments = "precompute.XX" %in% names(x$design),
                     optimizer = "FS", init = init, n.iter = x$opt$control["n.iter"], tol.score = x$opt$control["tol.score"], tol.param = x$opt$control["tol.param"], trace = trace)

    x$opt[c("cv","n.iter","score","previous.estimate")] <- eee[c("cv","n.iter","score","previous.estimate")]
    x$param <- eee$estimate
    x$logLik <- eee$logLik
    x[c("reparametrize","fitted","residuals","Omega","OmegaM1","dOmega","d2Omega","score","information","vcov","df","dVcov")] <- NULL

    class(x) <- append("clmm",class(x))
    return(x)
}

##----------------------------------------------------------------------
### constrain.R ends here
