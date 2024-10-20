### constrain.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 17 2022 (05:36) 
## Version: 
## Last-Updated: okt 20 2024 (16:43) 
##           By: Brice Ozenne
##     Update #: 145
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .constrain.lmm
##' @description Fit a linear mixed model where some parameters are set to pre-defined value.
##'
##' @param x lmm object
##' @param effects [numeric vector] value at which parameters should be set.
##' The name of the vector indicate the parameter.
##' @param transform.sigma,transform.k,transform.rho [character] transformation after which the constrained should be applied.
##' @param init [numeric vector] values of the (untransformed) parameters used to initialize the optimization procedure.
##' Does not need to include the constrained defined in \code{effect}: they will be included by the function.
##' @param trace [logical]
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details The argument \code{transform.sigma}, \code{transform.k}, \code{transform.rho} are used to make sense of the parameter name in argument \code{effects}.
##' For instance if \code{effects=c("sigma^2"=1)} then \code{transform.sigma} should be set to \code{"square"}.
##' 
##' @noRd
.constrain.lmm <- function(x, effects, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL,
                           init = NULL, trace = FALSE, ...){

    
    ## ** normalize user input

    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** transform
    if(is.null(transform.rho)){
        transform.rho <- "none"
    }
    if(is.null(transform.k) && (!is.null(transform.rho) && transform.rho == "none")){
        transform.k <- "none"
    }
    if(is.null(transform.sigma) && (!is.null(transform.k) && transform.k == "none") && (!is.null(transform.rho) && transform.rho == "none")){
        transform.sigma <- "none"
    }
    transform.init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = x$reparametrize$transform.sigma, x.transform.k = x$reparametrize$transform.k, x.transform.rho = x$reparametrize$transform.rho)
    transform.sigma <- transform.init$transform.sigma
    transform.k <- transform.init$transform.k
    transform.rho <- transform.init$transform.rho
    
    table.param <- stats::model.tables(x,
                                       effects = "param",
                                       transform.sigma = transform.sigma,
                                       transform.k = transform.k,
                                       transform.rho = transform.rho)

    ## *** init
    if(is.null(init)){
        init <- stats::coef(x, effects = "all", options = options)
    }else if(any(table.param$name %in% names(init) == FALSE)){
        stop("Incorrect argument \'init\': missing parameters ",paste(setdiff(table.param$name,names(init)), collapse = "\", \""),"\".\n")
    }else{
        init <- init[table.param$name]
    }

    ## *** effects
    name.effects <- names(effects)
    if(any(duplicated(name.effects))){
        stop("Incorrect argument \'effects\': contain duplicated names \"",paste(unique(name.effects[duplicated(name.effects)]), collapse = "\", \""),"\".\n")
    }
    if(all(name.effects %in% table.param$name)){
        init[name.effects] <- effects        
        constraint.transform <- FALSE
    }else if(all(name.effects %in% table.param$trans.name)){
        
        name.effectsMu <- intersect(name.effects, table.param[table.param$type=="mu","name"])
        if(length(name.effectsMu)>0){
            init[name.effectsMu] <- effects[name.effectsMu]
            constraint.transform <- FALSE
        }
        name.Omega <- stats::setNames(table.param[table.param$type!="mu","trans.name"], table.param[table.param$type!="mu","name"])
        name.effectsOmega <- intersect(name.effects, name.Omega)
        if(length(name.effectsOmega)>0){
            init2 <- coef(x, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, options = options)[name.Omega]
            init2[name.effectsOmega] <- effects[name.effectsOmega]
            init[names(name.Omega)] <- .reparametrize(stats::setNames(init2,names(name.Omega)),
                                                      type = stats::setNames(table.param$type,table.param$trans.name)[name.Omega],
                                                      sigma = stats::setNames(table.param$sigma,table.param$trans.name)[name.Omega],
                                                      k.x = stats::setNames(table.param$k.x,table.param$trans.name)[name.Omega],
                                                      k.y = stats::setNames(table.param$k.y,table.param$trans.name)[name.Omega],
                                                      Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                                      transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                                      transform.names = FALSE)$p
            constraint.transform <- TRUE
        }
    }else{
        stop("Incorrect argument \'effects\': parameter(s) \"",paste(setdiff(name.effects,table.param$trans.name), collapse = "\", \""),"\" not in the lmm.\n")
    }
    
    ## ** update model

    ## add constrain
    if(constraint.transform){
        x$design$param[match(name.effects,table.param$trans.name),"constraint"] <- effects
        
        test.k <- any(table.param[table.param$trans.name %in% name.effects,"type"]=="k" && transform.k %in% c("sd","logsd","var","logvar"))
        test.rho <- any(table.param[table.param$trans.name %in% name.effects,"type"]=="rho" && transform.rho %in% "cov")
        if(test.k || test.rho){
            x$reparametrize$transform.sigma <- transform.sigma
            x$reparametrize$transform.k <- transform.k
            x$reparametrize$transform.rho <- transform.rho
        }
    }else{
        x$design$param[match(name.effects,table.param$name),"constraint"] <- effects
    }

    ## refit
    eee <- .estimate(design = x$design, time = x$time, method.fit = x$args$method.fit, type.information = x$args$type.information,
                     transform.sigma = x$reparametrize$transform.sigma, transform.k = x$reparametrize$transform.k, transform.rho = x$reparametrize$transform.rho,
                     precompute.moments = "precompute.XX" %in% names(x$design),
                     optimizer = "FS", init = init,
                     n.iter = x$opt$control[["n.iter"]],
                     tol.score = x$opt$control[["tol.score"]],
                     tol.param = x$opt$control[["tol.param"]],
                     n.backtracking = x$opt$control[["n.backtracking"]],
                     init.cor = x$opt$control[["init.cor"]],
                     trace = trace)

    x$opt[c("cv","n.iter","score","previous.estimate")] <- eee[c("cv","n.iter","score","previous.estimate")]
    x$param <- eee$estimate
    x$logLik <- eee$logLik
    x[c("reparametrize","fitted","residuals","Omega","OmegaM1","dOmega","d2Omega","score","information","vcov","df","dVcov")] <- NULL

    class(x) <- append("clmm",class(x))
    return(x)
}

##----------------------------------------------------------------------
### constrain.R ends here
