### ranef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 26 2022 (11:18) 
## Version: 
## Last-Updated: mar 15 2024 (17:20) 
##           By: Brice Ozenne
##     Update #: 528
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * ranef.lmm (documentation)
##' @title Estimate Random Effect From a Linear Mixed Model
##' @description Recover the random effects from the variance-covariance parameter of a linear mixed model.
##' @param object a \code{lmm} object.
##' @param effects [character] should the estimated random effects (\code{"mean"}) or the estimated variance/standard deviation of the random effects (\code{"variance"},\code{"std"}) be output?
##' @param scale [character] should the total variance, variance relative to each random effect, and residual variance be output (\code{"absolute"}).
##' Or the ratio of these variances relative to the total variance (\code{"relative"}).
##' @param p [numeric vector] value of the model coefficients to be used. Only relevant if differs from the fitted values.
##' @param se [logical] should standard error and confidence intervals be evaluated using a delta method?
##' Will slow down the execution of the function.
##' @param format [character] should each type of random effect be output in a data.frame (\code{format="long"})
##' @param transform [logical] should confidence intervals for the variance estimates (resp. relative variance estimates) be evaluated using a log-transform (resp. atanh transformation)?
##' @param simplify [logical] when relevant will convert list with a single element to vectors and omit unessential output.
##' @param ... for internal use.
##'
##' @details Consider the following mixed model:
##' \deqn{Y = X\beta + \epsilon = X\beta + Z\eta + \xi}
##' where the variance of \eqn{\epsilon} is denoted \eqn{\Omega},
##' the variance of \eqn{\eta} is denoted \eqn{\Omega_{\eta}},
##' and the variance of \eqn{\xi} is \eqn{\sigma^2 I} with \eqn{I} is the identity matrix. \cr
##' The random effets are estimating according to:
##' \deqn{E[Y|\eta] = \Omega_{\eta} Z^{t} \Omega^{-1} (Y-X\beta)}
##' 
##' @keywords methods
##' 
##' @return A data.frame or a list depending on the argument \code{format}.
##' 
##' @examples
##' if(require(nlme)){
##' data(gastricbypassL, package = "LMMstar")
##' 
##' ## random intercept
##' e.RI <- lmm(weight ~ time + (1|id), data = gastricbypassL)
##' ranef(e.RI, effects = "mean")
##' ranef(e.RI, effects = "variance")
##'
##' }

## * ranef.lmm (code)
##' @export
ranef.lmm <- function(object, effects = "mean", scale = "absolute", se = (format=="long"), transform = (effects %in% c("std","variance")),
                      p = NULL, newdata = NULL, format = "long", simplify = TRUE, ...){



    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    mycall <- match.call()
    if(!inherits(object$design$vcov,"RE")){
        stop("Cannot estimate random effects linear mixed models defined by covariance structure (argument \'structure\'). \n",
             "Consider adding random effects in the argument \'formula\' instead. \n")
    }
    effects <- match.arg(effects, c("mean","std","variance"))
    if(length(scale)==1 && scale == "all"){
        scale <- c("absolute","relative")
    }
    scale <- match.arg(scale, c("absolute","relative"), several.ok = TRUE)
    format <- match.arg(format, c("wide","long"))
    if(length(transform)!=1 || (!is.numeric(transform) && !is.logical(transform))){
        stop("Argument \'transform\' must be a logical value. \n")
    }
    if(length(se)!=1 || (!is.numeric(se) && !is.logical(se))){
        stop("Argument \'transform\' must be a logical value. \n")
    }
    if(!is.null(mycall$se) && se == TRUE && format == "wide"){
        se <- FALSE
        message("Argument \'se\' ignored when using the wide format. \n",
                "Considering setting argument \'format\' to \"long\". \n")
    }
    if(!is.null(mycall$scale) && effects == "mean"){
        message("Argument \'scale\' ignored when argument \'effects\' is set to \"mean\". \n")
    }

    table.param <- object$design$param
    param.name <- table.param$name
    if(!is.null(p)){
        if(any(duplicated(names(p)))){
            stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
        }
        if(any(param.name %in% names(p) == FALSE)){
            stop("Incorrect argument \'p\': missing parameter(s) \"",paste(param.name[param.name %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
        }
        new.p <- TRUE
        p <- p[param.name]
    }else{
        new.p <- FALSE
        p <- object$param
    }

    if(transform && se && effects == "mean"){
        transform <- FALSE
        message("Argument \'transform\' is ignore when evaluating confidence intervals for the random effects. \n")
    }

    if(any(object$design$vcov$name$cor[[1]] %in% c("total","relative"))){
        stop("The ranef function cannot be used when the variable(s) w.r.t. which the random are defined are called \"total\" or \"relative\". \n")
    }
    
    ## ** extract from object    
    ## param
    param.type <- stats::setNames(table.param$type,param.name)
    param.rho <- param.name[param.type=="rho"]
    param.sigma <- param.name[param.type=="sigma"]

    ## cluster 
    var.cluster <- object$cluster$var
    n.cluster <- object$design$cluster$n
    index.cluster <- object$design$index.cluster
    Vindex.cluster <- attr(index.cluster, "vectorwise")
    
    ## strata
    var.strata <- object$strata$var
    index.clusterStrata <- object$design$index.clusterStrata
    U.strata <- object$strata$levels
    n.strata <- length(U.strata)

    ## design
    X.cor <- object$design$vcov$cor$X
    Xpattern.cor <- object$design$vcov$cor$Xpattern
    pattern.cluster <- object$design$vcov$pattern
    Upattern <- object$design$vcov$Upattern
    infoRanef <- object$design$vcov$ranef
    name.hierarchy <- unlist(infoRanef$hierarchy, use.names = FALSE)
    index.hierarchy <- unlist(lapply(1:length(infoRanef$hierarchy), function(iH){rep(iH,length(infoRanef$hierarchy[[iH]]))}))

    ## missing values
    index.na <- object$index.na

    ## all vars
    var.all <- unique(lava::manifest(object))
    
    ## ** converting correlation parameters into random effect variance
    ## convert coefficients into variance
    theta <- coef(object, p = p, effects = c("variance","correlation"), transform.sigma = "square", transform.rho = "cov", transform.names = FALSE)
    if(se){
        theta.var <- switch(effects,
                            "mean" = vcov(object, p = p, effects = c("all"), transform.sigma = "square", transform.rho = "cov", transform.names = FALSE),
                            vcov(object, p = p, effects = c("variance","correlation"), transform.sigma = "square", transform.rho = "cov", transform.names = FALSE)
                            )
    }

    ## *** prepare contrast matrix
    n.hierarchy <- length(infoRanef$param)
    index.hierarchy <- unlist(lapply(1:n.hierarchy, function(iH){rep(iH, length(infoRanef$param[[iH]]))}))
    name.RE <- unname(unlist(lapply(infoRanef$param, rownames)))
    n.RE <- length(name.RE)
    contrastRE <- array(0, dim = c(n.RE+2, length(theta), n.strata),
                        dimnames = list(c("total",name.RE,"residual"), names(theta), U.strata))

    ## total variance
    tableSigma.param <- table.param[match(param.sigma,table.param$name),,drop=FALSE]
    if(n.strata==1){
        contrastRE["total",param.sigma,U.strata] <- 1
    }else{
        diag(contrastRE["total",tableSigma.param[order(tableSigma.param$index.strata),"name"],U.strata]) <- 1
    }
    ## change in variance within hierarchy
    for(iH in 1:n.hierarchy){ ## iH <- 1
        iHierarchy <- infoRanef$param[[iH]]
        iContrast <- diag(1, nrow = NROW(iHierarchy), ncol = NROW(iHierarchy))        
        iContrast[lower.tri(iContrast, diag = FALSE)] <- -1
        for(iS in 1:n.strata){ ## iS <- 1
            contrastRE[rownames(iHierarchy),iHierarchy[,U.strata[iS]],U.strata[iS]] <- iContrast
        }        
    }
    ## residual variance
    if(length(name.RE)==1){
        contrastRE["residual",,] <- contrastRE["total",,] - contrastRE[name.RE,,]
    }else{
        for(iS in 1:n.strata){ ## iS <- 1
            contrastRE["residual",,iS] <- contrastRE["total",,iS] - colSums(contrastRE[name.RE,,iS])
        }
    }

    ## *** variance decomposition
    varDecomp <- apply(contrastRE, MARGIN = 3, FUN = `%*%`, theta)
    rownames(varDecomp) <- c("total",name.RE,"residual")
    if(any(varDecomp<=0)){
        stop("Variance for the random effects is found to be negative - cannot estimate the random effects. \n")
    }

    if(effects %in% c("std","variance")){

        ## absolute
        if(se){
            varDecomp.vcov <- apply(contrastRE, MARGIN = 3, FUN = function(iM){iM %*% theta.var %*% t(iM)}, simplify = FALSE)
            varDecomp.var <- do.call(cbind,lapply(varDecomp.vcov,diag))
        }else{
            varDecomp.var <- stats::setNames(vector(mode =  "numeric", length = n.strata), U.strata)
        }

        ## relative
        RvarDecomp <- sweep(varDecomp, FUN = "/", MARGIN = 2, STATS = varDecomp["total",])
        if(se){
            RvarDecomp.var <- do.call(cbind,lapply(1:n.strata, function(iS){ ## iS <- 1
                A <- varDecomp[-1,iS]
                B <- varDecomp["total",iS]
                c(total = 0,diag(varDecomp.vcov[[iS]])[-1]/B^2 + varDecomp.vcov[[iS]]["total","total"]*A^2/B^4 - 2*varDecomp.vcov[[iS]][-1,"total"]*A/B^3)
            }))
            colnames(RvarDecomp.var) <- U.strata
        }else{
            RvarDecomp.var <- stats::setNames(vector(mode =  "numeric", length = n.strata), U.strata)
        }

        ## take square root
        if(effects=="std"){
            varDecomp <- sqrt(varDecomp)
            RvarDecomp <- sqrt(RvarDecomp)            

            if(se){
                varDecomp.var <- varDecomp.var/(2*varDecomp)^2
                RvarDecomp.var <- RvarDecomp.var/(2*RvarDecomp)^2
            }
        }

        ## gather results into a single data.frame
        out <- do.call(rbind,mapply(u = as.data.frame(varDecomp), v = as.data.frame(varDecomp.var),
                                    x = as.data.frame(RvarDecomp), y = as.data.frame(RvarDecomp.var), s = colnames(varDecomp),
                                    FUN = function(u,v,x,y,s){ ## x <-
                                        iDF <- data.frame(strata = s,
                                                          type = rep(c("total",name.RE,"residual"),2),
                                                          scale = c(rep("absolute",n.RE+2),rep("relative",n.RE+2)),
                                                          estimate = c(u,x))
                                        if(se){iDF$se <- sqrt(c(v,y))}
                                        return(iDF)
                                    }, SIMPLIFY = FALSE))
        if(se & transform){
            out$lower = exp(log(out$estimate) + stats::qnorm(0.025)*out$se/out$estimate)
            out$upper = exp(log(out$estimate) + stats::qnorm(0.975)*out$se/out$estimate)
        }else if(se){
            out$lower <- out$estimate + stats::qnorm(0.025)*out$se
            out$upper <- out$estimate + stats::qnorm(0.975)*out$se
        }
        
        if(format=="long"){
            out <- out[out$scale %in% scale,,drop=FALSE]
            rownames(out) <- NULL
            if(simplify){
                if(length(scale)==1){
                    out$scale <- NULL
                }
                if(n.strata==1){
                    out$strata <- NULL
                }
                if(se == FALSE && n.strata==1 && length(scale)==1){
                    out <- stats::setNames(out$estimate,out$type)
                }
            }
        }else if(format=="wide"){
            if(!simplify || n.strata>1){
                out$scale <- paste(out$scale,out$strata,sep=".")
            }
            out$strata <- NULL
            out <- stats::reshape(out, direction  = "wide",
                                  idvar = "type", timevar = "scale",
                                  varying = unique(out$scale))
            rownames(out) <- NULL
        
        }
        return(out)

    }else if(effects == "mean"){

        ## gradient w.r.t. variance/covariance parameters
        contrastRE2 <- array(0, dim = c(length(name.RE),length(param.name),n.strata), dimnames = list(name.RE,param.name,U.strata))
        contrastRE2[,names(theta),] <- contrastRE[name.RE,,]

        df.epsilon <- stats::residuals(object, p = p, keep.data = TRUE, type = "normalized2", format = "long", simplify = FALSE)
            
        if(se){
            ## derivative of chol(Omega^(-1))(Y-XB)
            ## = chol(Omega^(-1))(Y-XB) - chol(Omega^(-1))X
            grad <- matrix(0, nrow = NROW(df.epsilon), ncol = length(p), dimnames = list(NULL, names(p)))
            
            if(!is.null(newdata) || new.p){
                design <- model.matrix(object, data = newdata, effects = "variance")
                Omega <- .calc_Omega(object = structure, param = p)
                precision <- lapply(Omega, solve)
            }else{
                precision <- object.OmegaM1
            }
        }
        

        
        sqrtPrecision$normalized <- lapply(Omega,function(iP){solve(chol(iP))})
        precision <- lapply(Omega, solve)
        resnorm <- as.double(iResidual %*% precision[[pattern[iId]]])
                    


        
        grid.ranef <- unlist(lapply(1:n.hierarchy, function(iH){ ## iH <- 1
            stats::setNames(lapply(1:length(infoRanef$hierarchy[[iH]]), function(iP){
                infoRanef$hierarchy[[iH]][1:iP]
            }), infoRanef$hierarchy[[iH]])
        }), recursive = FALSE)

        ls.out <- lapply(1:length(grid.ranef), function(iR){ ## iR <- 1
            do.call(rbind,by(df.epsilon, df.epsilon[grid.ranef[[iR]]], function(iDF){ ## iDF <- df.epsilon[df.epsilon$school == "I",]
                iNames <- grid.ranef[[iR]]
                iTau <- varDecomp[names(grid.ranef)[iR],which(iDF[[var.strata]][1]==U.strata)]
                iOut <- data.frame(variable = NA,
                                   strata = iDF[[var.strata]][1],
                                   level = NA,
                                   estimate = unname(sum(iDF$r.normalized)*iTau))
                iOut$variable <- list(iNames)
                iOut$level <- list(unlist(iDF[1,iNames,drop=FALSE]))
                return(iOut)
            }))
        })
        browser()
    }

    ## ** extract normalized residuals
    ## head(stats::residuals(object, p = p, keep.data = TRUE, type = "response", format = "long"))
    
    
    ## ** export
    out <- do.call(rbind,ls.out)
    rownames(out) <- NULL
    if(format == "wide"){
        out <- stats::setNames(lapply(1:n.hierarchy, function(iH){ ## iH <- 1
            iOut <- out[sapply(out$variable, utils::tail,1) %in% infoRanef$hierarchy[[iH]],,drop=FALSE]
            iOut$col <- sapply(iOut$level, function(iLevel){paste0(iLevel[-1], collapse = ":")})
            iOut$level <- sapply(iOut$level, utils::head, 1)
            iOutW <- stats::reshape(iOut, direction = "wide", 
                                    idvar = "level", timevar = "col", times = unique(iOut$col), drop = c("strata","variable"))
            colnames(iOutW)[2] <- "estimate"
            return(iOutW)
        }), sapply(infoRanef$hierarchy,"[",1))
        if(simplify && n.hierarchy == 1){
            names(out[[1]])[1] <- names(out)
            out <- out[[1]]            
        }
        
    }else if(format == "long"){
        if(simplify && n.strata == 1){
            out$strata <- NULL
        }
        if(simplify && all(lengths(infoRanef$hierarchy)==1)){
            out$variable <- unlist(out$variable)
            out$level <- unlist(out$level)
        }
    }
    return(out)
}

## * ranef.mlmm (code)
##' @export
ranef.mlmm <- function(object, ...){

    return(lapply(object$model, ranef, ...))
}

##----------------------------------------------------------------------
### ranef.R ends here
