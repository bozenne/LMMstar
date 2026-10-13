### information.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (22:13) 
## Version: 
## Last-Updated: okt 13 2026 (17:45) 
##           By: Brice Ozenne
##     Update #: 1348
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * information.lmm (documentation)
##' @title Extract The Information From a Linear Mixed Model
##' @description Extract or compute minus the second derivative of the log-likelihood of a linear mixed model.
##' 
##' @param x a \code{lmm} object.
##' @param newdata [data.frame] dataset relative to which the information should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param indiv [logical] Should the contribution of each cluster to the information be output? Otherwise output the sum of all clusters of the derivatives.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the information. Only relevant if differs from the fitted values.
##' @param effects [character] Should the information relative to all coefficients be output (\code{"all"} or \code{"fixed"}),
##' or only coefficients relative to the mean (\code{"mean"}),
##' or only coefficients relative to the variance and correlation structure (\code{"variance"} or \code{"correlation"}).
##' @param type.information [character] Should the expected information be computed (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details For details about the arguments \bold{transform.sigma}, \bold{transform.k}, \bold{transform.rho}, see the documentation of the \link[LMMstar]{coef.lmm} function.
##'
##' @return
##' When argument indiv is \code{FALSE}, a matrix with the value of the infroamtion relative to each pair of coefficient (in rows and columns) and each cluster (in rows).
##' When argument indiv is \code{TRUE}, a 3-dimensional array with the value of the information relative to each pair of coefficient (dimension 2 and 3) and each cluster (dimension 1).
##'
##' @keywords methods

## * information.lmm (code)
##' @export
information.lmm <- function(x, effects = NULL, newdata = NULL, p = NULL, indiv = FALSE, type.information = NULL,
                            transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

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

    ## *** effects
    if(is.null(effects)){
        if((is.null(transform.sigma) || identical(transform.sigma,"none")) && (is.null(transform.k) || identical(transform.k,"none")) && (is.null(transform.rho) || identical(transform.rho,"none"))){
            effects <- options$effects
        }else{
            effects <- c("mean","variance","correlation")
        }
    }else{
        if(!is.character(effects) || !is.vector(effects)){
            stop("Argument \'effects\' must be a character vector. \n")
        }
        valid.effects <- c("mean","fixed","variance","correlation","all")
        if(any(effects %in% valid.effects == FALSE)){
            stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
                 "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
        }
        if(all("all" %in% effects)){
            if(length(effects)>1){
                stop("Argument \'effects\' must have length 1 when containing the element \"all\". \n")
            }else{
                effects <- c("mean","variance","correlation")
            }
        }else{
            effects[effects == "fixed"] <- "mean"
        }
    }

    x.param <- stats::model.tables(x, effects = c("param",effects), transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)

    ## *** type.information
    if(is.null(type.information)){
        type.information <- x$args$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    ## *** transformation & p
    init <- .init_transform(p = p, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = x$reparametrize$transform.sigma, x.transform.k = x$reparametrize$transform.k, x.transform.rho = x$reparametrize$transform.rho,
                            table.param = x$design$param)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform
    if(is.null(p)){
        theta <- x$param
    }else{
        theta <- init$p
    }
    
    ## ** extract or recompute information
    if(is.null(newdata) && is.null(p) && (indiv == FALSE) && test.notransform && x$args$type.information==type.information){
        keep.name <- stats::setNames(x.param$name, x.param$trans.name)    

        design <- x$design ## useful in case of NA
        out <- x$information[keep.name,keep.name,drop=FALSE]
        attr(out,"type.information") <- attr(x$information,"type.information")
        if(transform.names){
            dimnames(out) <- list(names(keep.name),names(keep.name))
        }
    }else{
         
        if(!is.null(newdata)){
            design <- stats::model.matrix(x, newdata = newdata, effects = "all", simplify = FALSE)
        }else{
            design <- x$design
        }

        out <- .moments.lmm(value = theta, design = design, time = x$time, method.fit = x$args$method.fit, type.information = type.information,
                            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                            logLik = FALSE, score = FALSE, information = TRUE, vcov = FALSE, df = FALSE, indiv = indiv, effects = effects, robust = FALSE,
                            trace = FALSE, transform.names = transform.names)$information
    }

    ## ** restaure NAs and name
    if(indiv){

        if(!is.numeric(x$cluster$levels)){
            dimnames(out)[[1]] <- x$cluster$levels[match(1:dim(out)[[1]],x$cluster$index)]

        } 
        out <- restaureNA(out, index.na = x$index.na,
                          level = "cluster", cluster = x$cluster)
        
    }

    ## ** re-order values when converting to sd with strata (avoid sd0:0 sd0:1 sd1:0 sd1:1 sd2:0 sd2:1 ...)
    if("variance" %in% effects && transform.k %in% c("sd","var","logsd","logvar") && x$strata$n>1 && transform.names){
        if(indiv){
            out <- out[,x.param$trans.name,x.param$trans.name,drop=FALSE]
        }else{
            out <- out[x.param$trans.name,x.param$trans.name,drop=FALSE]
        }
    }

    ## ** export
    return(out)
}

## * .information
## REML term
## d 0.5 tr[(X \OmegaM1 X)^{-1} (X \OmegaM1 d\Omega \OmegaM1 X)] = 0.5 tr[ (X \OmegaM1 d'\Omega \OmegaM1 X) (X \OmegaM1 X)^{-2} (X \OmegaM1 d\Omega \OmegaM1 X) ]
##                                                                 - 0.5 tr[ (X \OmegaM1 X)^{-1} (X \OmegaM1 d'\Omega \OmegaM1 d\Omega \OmegaM1 X) + (X \OmegaM1 X)^{-1} (X \OmegaM1 d\Omega \OmegaM1 d'\Omega \OmegaM1 X) ]
##                                                                 + 0.5 tr[ (X \OmegaM1 X)^{-1} (X \OmegaM1 d2\Omega \OmegaM1 X) ]
.information <- function(X, residuals, weights, 
                         pattern, index.cluster, name.allcoef, pair.vcov,
                         indiv, REML, type.information, effects, precompute){

    ## ** extract information
    n.cluster <- length(pattern)
    n.mucoef <- NCOL(X)
    
    name.mucoef <- colnames(X)
    name.vcovcoef <- intersect(setdiff(name.allcoef,name.mucoef), unlist(lapply(precompute$Omega$dOmegaM1,function(iA){dimnames(iA)[[2]]})))
    n.vcovcoef <- length(name.vcovcoef)
    npair.vcov <- NROW(pair.vcov) ## number of pairs of variance parameters
    
    ## ** prepare output
    name.effects <- c(name.mucoef, name.vcovcoef)
    n.effects <- length(name.effects)
    if(indiv){
        info <- array(0, dim = c(n.cluster, n.effects, n.effects),
                      dimnames = list(NULL, name.effects, name.effects))
    }else{
        info <- matrix(0, nrow = n.effects, ncol = n.effects,
                       dimnames = list(name.effects, name.effects)
                       )
    }    

    ## ** restrict to relevant parameters
    test.mean <- (n.mucoef>0) && ("mean" %in% effects)
    test.vcov <- (n.vcovcoef>0) && (("variance" %in% effects) || ("correlation" %in% effects))

    if(test.vcov){
        if(REML && indiv){
            message <- "approximate individual REML information"
        }else{
            message <- NULL
        }

        if("variance" %in% effects == FALSE || "correlation" %in% effects == FALSE){            
            if("variance" %in% effects == FALSE){            
                name.vcovcoef <- setdiff(name.vcovcoef, name.allcoef[attr(name.allcoef, "type") %in% c("sigma","k")])
            }else if("correlation" %in% effects == FALSE){
                name.vcovcoef <- setdiff(name.vcovcoef, name.allcoef[attr(name.allcoef, "type") %in% c("rho")])
            }
            precompute$Omega$d2Trace <- lapply(precompute$Omega$d2Trace, function(iO){iO[intersect(names(iO),name.effects),,drop=FALSE]})
            if(type.information=="observed"){
                precompute$Omega$dOmegaM1 <- lapply(precompute$Omega$dOmegaM1, function(iO){iO[,intersect(colnames(iO),name.effects),,drop=FALSE]})
                precompute$Omega$d2OmegaM1 <- lapply(precompute$Omega$d2OmegaM1, function(iO){iO[,intersect(colnames(iO),name.effects),,drop=FALSE]})
            }
            n.vcovcoef <- length(name.vcovcoef)
            test.vcov <- n.vcovcoef>0
        }        
    }else{
        message <- NULL
    }

    ## ** compute information
    if(indiv){
        
        ## *** looping over individuals
        for(iId in 1:n.cluster){ ## iId <- 7
            iIndex <- index.cluster[[iId]]
            iPattern <- pattern[iId]
            if(type.information == "observed"){
                iResidual <- residuals[iIndex,,drop=FALSE]
            }
            iX <- t(X[iIndex,,drop=FALSE])
            iOmegaM1 <- precision[[iPattern]] * weights[iId,"Omega"]
            
            ## **** mean,mean
            if(test.mean){
                info[iId,name.mucoef,name.mucoef] <- weights[iId,"likelihood"] * (iX %*% iOmegaM1 %*% t(iX))
            }
            
            ## **** var,var
            if(test.vcov){
                ## subset in case effects is only variance or correlation
                iname.vcovcoef <- intersect(names(dOmega[[iPattern]]), name.vcovcoef)
                iPair.vcov <- pair.vcov[which(pair.vcov[[iPattern]]),,drop=FALSE]
                
                ## compute and store contribution
                if(type.information == "expected"){
                    iValue <- 0.5 * weights[iId,"likelihood"] * precompute$Omega$tr.OmegaM1.dOmega.OmegaM1.dOmega[[iPattern]]
                }else if(type.information == "observed"){
                    iValue <- - 0.5 * weights[iId,"likelihood"] * (precompute$Omega$tr.OmegaM1.dOmega.OmegaM1.dOmega[[iPattern]] - precompute$Omega$tr.OmegaM1.d2Omega[[iPattern]])
                    iValue <- iValue - 0.5 * prod(weights[iId,]) * (as.vector(tcrossprod(iResidual)) %*% precompute$Omega$OmegaM1.d2OmegaAndCo.OmegaM1[[iPattern]])[1,]
                    ## same as iResidual[,1] %*% matrix(precompute$Omega$OmegaM1.d2OmegaAndCo.OmegaM1[[iPattern]][,2], nrow = length(iIndex), ncol = length(iIndex)) %*% iResidual 
                }
                iKey <- c(iPair.vcov$key1,stats::na.omit(iPair.vcov$key2)) + (iId-1)*prod(dim(info)[2:3])
                info[iKey]  <- info[iKey] + c(iValue,iValue[!is.na(iPair.vcov$key2)])
                
                if(REML){
                    ## APPROXIMATION: the REML score w.r.t. variance parameter is not linear in the individual contribution
                    for(iP2 in 1:NROW(iPair.vcov)){ ## iP2 <- 1
                        iParam2.1 <- iPair.vcov[iP2,"param1"]
                        iParam2.2 <- iPair.vcov[iP2,"param2"]
                        iName <- iPair.vcov[iP2,"param"]
                        
                        iREML.num1.1 <- prod(weights[iId,]) * iX %*% OmegaM1.dOmega.OmegaM1[[iPattern]][,,iParam2.1] %*% t(iX)
                        ## shortcut for iOmegaM1 %*% dOmega[[iPattern]][[iParam2.1]] %*% iOmegaM1
                        iREML.num1.2 <- prod(weights[iId,]) * iX %*% OmegaM1.dOmega.OmegaM1[[iPattern]][,,iParam2.2] %*% t(iX)
                        ## shortcut for iOmegaM1 %*% dOmega[[iPattern]][[iParam2.2]] %*% iOmegaM1
                        iREML.num2 <- prod(weights[iId,]) * iX %*% OmegaM1.d2OmegaAndCo.OmegaM1[[iPattern]][,,iName] %*% t(iX)
                        ## shortcut for iOmegaM1 %*% (d2Omega[[iPattern]][[iParam2]] - 2 * dOmega[[iPattern]][[iParam2.1]] %*% iOmegaM1 %*% dOmega[[iPattern]][[iParam2.2]]) %*% iOmegaM1
                        
                        iValue <- 0.5 * sum(REML.denom * (iREML.num2 + 0.5 * iREML.num1.1 %*% REML.denom %*% REML.num1[[iParam2.2]] + 0.5 * REML.num1[[iParam2.1]] %*% REML.denom %*% iREML.num1.2))

                        info[iId,iParam2.1,iParam2.2] <- info[iId,iParam2.1,iParam2.2] - iValue
                        if(iParam2.1 != iParam2.2){
                            info[iId,iParam2.1,iParam2.2] <- info[iId,iParam2.2,iParam2.1] - iValue
                        }
                    }
                }
            }

            ## **** mean,var
            if(type.information == "observed" && test.mean && test.vcov){

                ## compute
                iValue <- prod(weights[iId,]) * t(iResidual %x% t(iX)) %*% precompute$Omega$OmegaM1.dOmega.OmegaM1[[iPattern]][,iname.vcovcoef,drop=FALSE]
                ## kronecker product is the same as apply(iX, MARGIN = 1, FUN = function(iCol){as.double(tcrossprod(iCol,iResidual))})
                
                ## store
                info[iId,name.mucoef,iname.vcovcoef] <- iValue
                info[iId,iname.vcovcoef,name.mucoef] <- t(iValue)

            }
        }
    }else{
        ## *** looping over covariance patterns
        for (iG in 1:NROW(precompute$weights)) { ## iG <- 1

            iPattern <- precompute$weights[iG,"pattern"]
            browser()

            ## **** mean,mean
            if(test.mean){
                iOmegaM1 <- precompute$Omega$OmegaM1[[iPattern]]
                iXX <- precompute$XX$pattern[[iPattern]]## has already been weigthed
                info[name.mucoef,name.mucoef] <- info[name.mucoef,name.mucoef] + (as.vector(iOmegaM1) %*% iXX)[as.double(precompute$XX$key)]                
            }

            ## **** var,var
            if(test.vcov){
                ## compute and store contribution
                if(type.information == "expected"){
                    iValue <- 0.5 * precompute$weights[iPattern] * precompute$Omega$tr.OmegaM1.dOmega.OmegaM1.dOmega[[iPattern]]
                }else if(type.information == "observed"){                    
                    iValue <- - 0.5 * precompute$weights[iPattern] * (precompute$Omega$tr.OmegaM1.dOmega.OmegaM1.dOmega[[iPattern]] - precompute$Omega$tr.OmegaM1.d2Omega[[iPattern]])
                    iValue <- iValue - 0.5 * (precompute$RR[[iPattern]] %*% precompute$Omega$OmegaM1.d2OmegaAndCo.OmegaM1[[iPattern]])[1,]
                }
                iKey <- c(pair.vcov[pair.vcov[[iPattern]],"key1"],stats::na.omit(pair.vcov[pair.vcov[[iPattern]],"key2"]))
                info[iKey]  <- info[iKey] + c(iValue,iValue[!is.na(pair.vcov[pair.vcov[[iPattern]],"key2"])])
                ## any(duplicated(sort(c(pair.vcov$key1,stats::na.omit(pair.vcov$key2)))))
            }

            ## **** mean,var
            if(type.information == "observed" && test.mean && test.vcov){
                ## compute
                iValue <- t(precompute$Omega$OmegaM1.dOmega.OmegaM1[[iPattern]]) %*% precompute$XR[[iPattern]]

                ## store
                info[iName.vcovcoef,name.mucoef] <- info[iName.vcovcoef,name.mucoef] + iValue
                info[name.mucoef,iName.vcovcoef] <- info[name.mucoef,iName.vcovcoef] + t(iValue)
            }
        }

        ## *** REML contribution
        if(REML && test.vcov){

            for(iP2 in 1:npair.vcov){ ## iP2 <- 1                
                iParam2.1 <- pair.vcov[iP2,"param1"]
                iParam2.2 <- pair.vcov[iP2,"param2"]
                iName <- pair.vcov[iP2,"name"]
                ## same at 0.5 * tr(REML.denom %*% (REML.num2[[iParam2]] + REML.num1[[iParam2.1]] %*% REML.denom %*% REML.num1[[iParam2.2]]))
                iValue <- 0.5 * sum(REML.denom * (REML.num2[[iName]] + REML.num1[[iParam2.1]] %*% REML.denom %*% REML.num1[[iParam2.2]]))
                info[iParam2.1,iParam2.2] <- info[iParam2.1,iParam2.2] - iValue
                if(iParam2.1 != iParam2.2){
                    info[iParam2.2,iParam2.1] <- info[iParam2.2,iParam2.1] - iValue
                }
            }
        }

    }
    ## ** export
    attr(info,"message") <- message
    return(info)
}

##----------------------------------------------------------------------
### information.R ends here

