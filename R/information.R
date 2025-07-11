### information.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (22:13) 
## Version: 
## Last-Updated: jul  9 2025 (14:00) 
##           By: Brice Ozenne
##     Update #: 1260
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
                            trace = FALSE, precompute.moments = !is.null(x$design$precompute.XX), transform.names = transform.names)$information
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
.information <- function(X, residuals, precision, dOmega, d2Omega, weights, 
                         pattern, index.cluster, name.allcoef,
                         pair.meanvcov, pair.vcov, indiv, REML, type.information, effects, 
                         precompute){

    ## ** extract information
    n.cluster <- length(pattern)
    name.mucoef <- colnames(X)
    n.mucoef <- length(name.mucoef)
    name.varcoef <- unique(unlist(lapply(dOmega,names)))
    n.varcoef <- length(name.varcoef)    
    name.varcoef2 <- colnames(attr(pair.vcov,"global"))
    n.varcoef2 <- length(name.varcoef2)
    U.pattern <- names(dOmega)

    ## ** prepare output
    compute.indiv <- indiv || is.null(precompute$weights) || is.null(precompute$XR) || is.null(precompute$RR)
    name.effects <- attr(effects,"original.names")
    n.effects <- length(name.effects)
    if(compute.indiv){
        info <- array(0, dim = c(n.cluster, n.effects, n.effects),
                      dimnames = list(NULL, name.effects, name.effects))
    }else if(any(sapply(precision, inherits, "try-error"))){ ## when evaluating score at parameter values where the residual variance-covariance matrix is singular
        return(matrix(NA, nrow = n.effects, ncol = n.effects,
                      dimnames = list(name.effects, name.effects)
                      ))
    }else{
        info <- matrix(0, nrow = n.effects, ncol = n.effects,
                       dimnames = list(name.effects, name.effects)
                       )
    }    

    ## ** restrict to relevant parameters
    test.mean <- (n.mucoef>0) && ("mean" %in% effects)
    test.vcov <- (n.varcoef>0) && (("variance" %in% effects) || ("correlation" %in% effects))

    if(test.vcov){
        if(REML && indiv){
            message <- "approximate individual REML information"
        }else{
            message <- NULL
        }

        if(("variance" %in% effects == FALSE) || ("correlation" %in% effects == FALSE)){## restrict to requested coefficients
            name.varcoef <- intersect(name.varcoef,name.effects)
            n.varcoef <- length(name.varcoef)
            precompute$Omega$tr.OmegaM1.dOmega <- lapply(precompute$Omega$tr.OmegaM1.dOmega, function(iO){iO[intersect(names(iO),name.varcoef)]})
            precompute$Omega$OmegaM1.dOmega.OmegaM1 <- lapply(precompute$Omega$OmegaM1.dOmega.OmegaM1, function(iO){iO[,intersect(colnames(iO),name.varcoef),drop=FALSE]})
            test.vcov <- n.varcoef>0

            cbindpair.vcov <- attr(pair.vcov,"global")
            cbindpair.vcov[] <- cbindpair.vcov %in% name.varcoef
            name.varcoef2 <- names(which(colSums(cbindpair.vcov==TRUE)==2))
            n.varcoef2 <- length(name.varcoef2)
            precompute$Omega$tr.OmegaM1.dOmega.OmegaM1.dOmega <- lapply(precompute$Omega$tr.OmegaM1.dOmega.OmegaM1.dOmega, function(iO){iO[,intersect(colnames(iO),name.varcoef2),drop=FALSE]})
            if(type.information == "observed"){
                precompute$Omega$OmegaM1.d2OmegaAndCo.OmegaM1 <- lapply(precompute$Omega$OmegaM1.d2OmegaAndCo.OmegaM1, function(iO){iO[,intersect(colnames(iO),name.varcoef2),drop=FALSE]})
                precompute$Omega$Omega$tr.OmegaM1.d2Omega <- lapply(precompute$Omega$Omega$tr.OmegaM1.d2Omega, function(iO){iO[,intersect(colnames(iO),name.varcoef2),drop=FALSE]})
            }
        }
    }else{
        message <- NULL
    }

    ## ** prepare REML term
    if(REML && test.vcov){
        if(is.null(precompute$REML)){
            X.OmegaM1.X <- matrix(0, nrow = n.mucoef, ncol = n.mucoef)
            REML.num1 <- stats::setNames(replicate(n.varcoef, matrix(0, nrow = n.mucoef, ncol = n.mucoef), simplify = FALSE), name.varcoef)
            REML.num2 <- stats::setNames(replicate(n.varcoef2, matrix(0, nrow = n.mucoef, ncol = n.mucoef), simplify = FALSE), name.varcoef2)
            for(iId in 1:n.cluster){ ## iId <- 1
                iX <- X[index.cluster[[iId]],,drop=FALSE]
                iOmegaM1 <- precision[[pattern[iId]]]
                iOmegaM1.dOmega.OmegaM1 <- precompute$Omega$OmegaM1.dOmega.OmegaM1[[pattern[iId]]]
                iOmegaM1.d2OmegaAndCo.OmegaM1 <- precompute$Omega$OmegaM1.d2OmegaAndCo.OmegaM1[[pattern[iId]]]
                iWeights <- weights[iId]
                X.OmegaM1.X <- X.OmegaM1.X + iWeights * t(iX) %*% iOmegaM1 %*% iX
                for(iParam in intersect(names(dOmega[[pattern[iId]]]), name.varcoef)){ ## intersect to handle when argument effects is only "variance" or "correlation"
                    REML.num1[[iParam]] <- REML.num1[[iParam]] + iWeights * t(iX) %*% iOmegaM1.dOmega.OmegaM1[[iParam]] %*% iX
                }
                for(iParam2 in intersect(names(d2Omega[[pattern[iId]]]), name.varcoef2)){ ## intersect to handle when argument effects is only "variance" or "correlation"
                    REML.num2[[iParam2]] <- REML.num2[[iParam2]] + iWeights * t(iX) %*% iOmegaM1.d2OmegaAndCo.OmegaM1[[iParam2]] %*% iX
                }
            }
            REML.denom <- solve(X.OmegaM1.X)
        }else{
            REML.num2 <- precompute$REML$X.OmegaM1.d2OmegaAndCo.OmegaM1.X
            REML.num1 <- precompute$REML$X.OmegaM1.dOmega.OmegaM1.X
            REML.denom <- precompute$REML$X.OmegaM1.X_M1
        }
    }

    ## ** compute information
    if(compute.indiv){
        
        ## *** looping over individuals
        if(REML && test.vcov){
            OmegaM1.dOmega.OmegaM1 <- lapply(precompute$Omega$OmegaM1.dOmega.OmegaM1, function(iO){
                array(iO, dim = c(sqrt(NROW(iO)),sqrt(NROW(iO)),NCOL(iO)), dimnames = list(NULL,NULL,colnames(iO)))
            })
            OmegaM1.d2OmegaAndCo.OmegaM1 <- lapply(precompute$Omega$OmegaM1.d2OmegaAndCo.OmegaM1, function(iO){
                array(iO, dim = c(sqrt(NROW(iO)),sqrt(NROW(iO)),NCOL(iO)), dimnames = list(NULL,NULL,colnames(iO)))
            })
        }
        
        for(iId in 1:n.cluster){ ## iId <- 7
            iIndex <- index.cluster[[iId]]
            iPattern <- pattern[iId]
            iWeights <- weights[iId]
            if(type.information == "observed"){
                iResidual <- residuals[iIndex,,drop=FALSE]
            }
            iX <- t(X[iIndex,,drop=FALSE])
            iOmegaM1 <- precision[[iPattern]]
        
            ## **** mean,mean
            if(test.mean){
                info[iId,name.mucoef,name.mucoef] <- iWeights * (iX %*% iOmegaM1 %*% t(iX))
            }
            
            ## **** var,var
            if(test.vcov){
                ## subset in case effects is only variance or correlation
                iName.varcoef <- intersect(names(dOmega[[iPattern]]), name.varcoef)

                ## compute and store contribution
                if(type.information == "expected"){
                    iValue <- 0.5 * iWeights * precompute$Omega$tr.OmegaM1.dOmega.OmegaM1.dOmega[[iPattern]]
                }else if(type.information == "observed"){
                    iValue <- - 0.5 * iWeights * (precompute$Omega$tr.OmegaM1.dOmega.OmegaM1.dOmega[[iPattern]] - precompute$Omega$tr.OmegaM1.d2Omega[[iPattern]])
                    iValue <- iValue - 0.5 * (as.vector(tcrossprod(iResidual)) %*% precompute$Omega$OmegaM1.d2OmegaAndCo.OmegaM1[[iPattern]])[1,]
                    ## same as iResidual[,1] %*% matrix(precompute$Omega$OmegaM1.d2OmegaAndCo.OmegaM1[[iPattern]][,2], nrow = length(iIndex), ncol = length(iIndex)) %*% iResidual 
                }
                info[iId,iName.varcoef,iName.varcoef]  <- info[iId,iName.varcoef,iName.varcoef] + iValue[as.vector(attr(pair.vcov[[iPattern]],"key"))]

                if(REML){
                    ## APPROXIMATION: the REML score w.r.t. variance parameter is not linear in the individual contribution
                    for(iParam2 in name.varcoef2){ ## iParam2 <- name.varcoef2[1]
                        iParam2.1 <- attr(pair.vcov,"global")[1,iParam2]
                        iParam2.2 <- attr(pair.vcov,"global")[2,iParam2]
                        
                        iREML.num1.1 <- iWeights * iX %*% OmegaM1.dOmega.OmegaM1[[iPattern]][,,iParam2.1] %*% t(iX)
                        ## shortcut for iOmegaM1 %*% dOmega[[iPattern]][[iParam2.1]] %*% iOmegaM1
                        iREML.num1.2 <- iWeights * iX %*% OmegaM1.dOmega.OmegaM1[[iPattern]][,,iParam2.2] %*% t(iX)
                        ## shortcut for iOmegaM1 %*% dOmega[[iPattern]][[iParam2.2]] %*% iOmegaM1
                        iREML.num2 <- iWeights * iX %*% OmegaM1.d2OmegaAndCo.OmegaM1[[iPattern]][,,iParam2] %*% t(iX)
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
                iValue <- t(iResidual %x% t(iX)) %*% precompute$Omega$OmegaM1.dOmega.OmegaM1[[iPattern]][,iName.varcoef,drop=FALSE]
                ## kronecker product is the same as apply(iX, MARGIN = 1, FUN = function(iCol){as.double(tcrossprod(iCol,iResidual))})
                
                ## store
                info[iId,name.mucoef,iName.varcoef] <- iValue
                info[iId,iName.varcoef,name.mucoef] <- t(iValue)

            }
        }
    }else{
        ## *** looping over covariance patterns
        for (iPattern in U.pattern) { ## iPattern <- name.pattern[1]
            
            iOmegaM1 <- precision[[iPattern]]
            iTime2 <- length(iOmegaM1)
            iName.varcoef <- names(dOmega[[iPattern]])

            ## **** mean,mean
            if(test.mean){
                info[name.mucoef,name.mucoef] <- info[name.mucoef,name.mucoef] + (attr(iOmegaM1,"vectorize") %*% precompute$XX$pattern[[iPattern]])[as.double(precompute$XX$key)]                
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
                info[iName.varcoef,iName.varcoef]  <- info[iName.varcoef,iName.varcoef] + iValue[as.vector(attr(pair.vcov[[iPattern]],"key"))]

            }

            ## **** mean,var
            if(type.information == "observed" && test.mean && test.vcov){
                ## compute
                iValue <- t(precompute$Omega$OmegaM1.dOmega.OmegaM1[[iPattern]]) %*% precompute$XR[[iPattern]]

                ## store
                info[iName.varcoef,name.mucoef] <- info[iName.varcoef,name.mucoef] + iValue
                info[name.mucoef,iName.varcoef] <- info[name.mucoef,iName.varcoef] + t(iValue)
            }
        }

        ## *** REML contribution
        if(REML && test.vcov){

            for(iParam2 in name.varcoef2){ ## iParam2 <- name.varcoef2[4]
                iParam2.1 <- attr(pair.vcov,"global")[1,iParam2]
                iParam2.2 <- attr(pair.vcov,"global")[2,iParam2]
                ## same at 0.5 * tr(REML.denom %*% (REML.num2[[iParam2]] + REML.num1[[iParam2.1]] %*% REML.denom %*% REML.num1[[iParam2.2]]))
                iValue <- 0.5 * sum(REML.denom * (REML.num2[[iParam2]] + REML.num1[[iParam2.1]] %*% REML.denom %*% REML.num1[[iParam2.2]]))
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

