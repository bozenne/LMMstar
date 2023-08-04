### information.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (22:13) 
## Version: 
## Last-Updated: aug  3 2023 (18:08) 
##           By: Brice Ozenne
##     Update #: 1167
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
##' @description Extract or compute the (expected) second derivative of the log-likelihood of a linear mixed model.
##' 
##' @param x a \code{lmm} object.
##' @param data [data.frame] dataset relative to which the information should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param indiv [logical] Should the contribution of each cluster to the information be output? Otherwise output the sum of all clusters of the derivatives.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the information. Only relevant if differs from the fitted values.
##' @param effects [character] Should the information relative to all coefficients be output (\code{"all"} or \code{"fixed"}),
##' or only coefficients relative to the mean (\code{"mean"}),
##' or only coefficients relative to the variance and correlation structure (\code{"variance"} or \code{"correlation"}).
##' @param type.information [character] Should the expected information be computed  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
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
information.lmm <- function(x, effects = NULL, data = NULL, p = NULL, indiv = FALSE, type.information = NULL,
                            transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

    ## ** normalize user input
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(is.null(type.information)){
        robust <- FALSE
        type.information <- x$args$type.information
    }else{
        robust <- identical(attr(type.information,"robust"),TRUE)
        type.information <- match.arg(type.information, c("expected","observed"))
    }    
    if(is.null(effects)){
        effects <- options$effects
    }else if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","fixed","variance","correlation"), several.ok = TRUE)
    effects[effects== "fixed"] <- "mean"

    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = x$reparametrize$transform.sigma, x.transform.k = x$reparametrize$transform.k, x.transform.rho = x$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform
    
    ## ** extract or recompute information
    if(is.null(data) && is.null(p) && (indiv == FALSE) && test.notransform && (robust==FALSE) && x$args$type.information==type.information){
        keep.name <- stats::setNames(names(coef(x, effects = effects, transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE)),
                                     names(coef(x, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)))    

        design <- x$design ## useful in case of NA
        out <- x$information[keep.name,keep.name,drop=FALSE]
        if(transform.names){
            dimnames(out) <- list(names(keep.name),names(keep.name))
        }
    }else{
        test.precompute <- !is.null(x$design$precompute.XX) && !indiv
         
        if(!is.null(data)){
            design <- stats::model.matrix(x, data = data, effects = "all", simplify = FALSE)
        }else{
            design <- x$design
        }

        if(!is.null(p)){
            if(any(duplicated(names(p)))){
                stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
            }
            if(any(names(x$param) %in% names(p) == FALSE)){
                stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(x$param)[names(x$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
            }
            p <- p[names(x$param)]
        }else{
            p <- x$param
        }
        out <- .moments.lmm(value = p, design = design, time = x$time, method.fit = x$args$method.fit, type.information = type.information,
                            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                            logLik = FALSE, score = FALSE, information = TRUE, vcov = FALSE, df = FALSE, indiv = indiv, effects = effects, robust = robust,
                            trace = FALSE, precompute.moments = test.precompute, transform.names = transform.names)$information
    }

    ## ** restaure NAs and name
    if(indiv){

        if(!is.numeric(x$cluster$levels)){
            dimnames(out)[[1]] <- x$cluster$levels[match(1:dim(out)[[1]],x$cluster$index)]

        } 
        out <- addNA(out, index.na = x$index.na,
                     level = "cluster", cluster = x$cluster)
        
    }


    ## ** re-order values when converting to sd with strata (avoid sd0:0 sd0:1 sd1:0 sd1:1 sd2:0 sd2:1 ...)
    if("variance" %in% effects && transform.k %in% c("sd","var","logsd","logvar") && x$strata$n>1 && transform.names){
        out.name <- names(stats::coef(x, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = TRUE))
        if(indiv){
            out <- out[,out.name,out.name,drop=FALSE]
        }else{
            out <- out[out.name,out.name,drop=FALSE]
        }
    }

    ## ** export
    return(out)
}

## * .information
.information <- function(X, residuals, precompute, pair.vcov,
                         Upattern.ncluster, weights, scale.Omega,
                         pattern, index.cluster, 
                         indiv, REML, type.information, effects, robust){

    ## ** extract information
    if(indiv && REML && c("variance","correlation") %in% effects){
        stop("Cannot compute cluster-specific information contributions relative to the variance-covariance parameters when using REML. \n")
    }
    test.loopIndiv <- indiv || !attr(precompute,"moments")

    ## all parameters requested by the user (maybe no variance but correlation)
    name.effects <- attr(effects,"original.names")
    n.effects <- length(name.effects)

    if(!indiv && any(is.na(precompute$Omega.logdet))){ ## non positive definite residual variance covariance
        return(matrix(NA, nrow = n.effects, ncol = n.effects,
                      dimnames = list(name.effects, name.effects)))
    }

    ## all free parameters
    name.mucoef <- intersect(colnames(X), name.effects)
    n.mucoef <- length(name.mucoef)

    name.varcoef <- lapply(precompute$dOmega_OmegaM1, function(iO){intersect(names(iO),name.effects)})
    name.allvarcoef <- unique(unlist(name.varcoef))
    n.allvarcoef <- length(name.allvarcoef)

    name.allcoef <- c(name.mucoef, name.allvarcoef)
    n.allcoef <- length(name.allcoef)
    
    ## ** compute information
    ## *** looping over individuals
    if(test.loopIndiv){
        
        n.cluster <- length(pattern)
        info <- array(0, dim = c(n.cluster, n.allcoef, n.allcoef),
                      dimnames = list(NULL, name.allcoef, name.allcoef))

        ## loop
        for(iId in 1:n.cluster){ ## iId <- 7
            iPattern <- pattern[iId]
            iIndex <- index.cluster[[iId]]
            iWeight <- weights[iId]

            iX <- X[iIndex,,drop=FALSE]
            tiX <- t(iX)
            iOmegaM1 <- precision[[iPattern]] * scale.Omega[iId]
            if(type.information == "observed"){
                iResidual <- residuals[iIndex,,drop=FALSE]
            }
        
            ## **** mean,mean
            iValue <-  iWeight * (tiX %*% iOmegaM1 %*% iX)
            if(test.mean){
                if(indiv){
                    info[iId,name.mucoef,name.mucoef] <- iValue
                }else{
                    info[name.mucoef,name.mucoef] <- info[name.mucoef,name.mucoef] + iValue
                }
            }
            if(REML && test.vcov){
                REML.denom <- REML.denom + iValue
                for(iVarcoef in name.varcoef[[iPattern]]){ ## iVarcoef <- 1
                    REML.numerator1[,,iVarcoef] <- REML.numerator1[,,iVarcoef] + iWeight * (tiX %*% ls.OmegaM1_dOmega_OmegaM1[[iPattern]][[iVarcoef]] %*% iX) * scale.Omega[iId]
                }
            }

            ## **** var,var
            if(test.vcov){

                for(iPair in 1:npair.vcov[[iPattern]]){ ## iPair <- 1
                    iCoef1 <- pair.vcov[[iPattern]][1,iPair]
                    iCoef2 <- pair.vcov[[iPattern]][2,iPair]

                    iValue <- 0.5 * iWeight * tr_OmegaM1_d2OmegaAndCo[[iPattern]][iPair]
                    ## 0.5 * ntr(iOmegaM1 %*% idOmega$sigma %*% iOmegaM1 %*% idOmega$sigma)

                    if(type.information == "observed"){
                        iValue <- iValue - 0.5 * iWeight * (t(iResidual) %*% OmegaM1_d2OmegaAndCo_OmegaM1[[iPattern]][,,iPair] %*% iResidual) * scale.Omega[iId]
                    }
                    if(indiv){
                        info[iId,iCoef1,iCoef2] <- iValue
                        if(iCoef1 != iCoef2){
                            info[iId,iCoef2,iCoef1] <- iValue
                        }
                    }else{
                        info[iCoef1,iCoef2] <- info[iCoef1,iCoef2] + iValue
                        if(iCoef1 != iCoef2){
                            info[iCoef2,iCoef1] <- info[iCoef2,iCoef1] + iValue
                        }
                    }

                    if(REML){
                        iKey <- REML.key[iCoef1,iCoef2]
                        REML.numerator2[,,iKey] <- REML.numerator2[,,iKey] + iWeight * (tiX %*% OmegaM1_d2OmegaAndCo_OmegaM1[[iPattern]][,,iPair] %*% iX) * scale.Omega[iId]
                    }
                }
            }

            ## **** mean,var
            if(type.information == "observed" && test.mean && test.vcov){

                for(iPair in 1:npair.meanvcov[[iPattern]]){ ## iPair <- 1
                    iCoef1 <- pair.meanvcov[[iPattern]][1,iPair]
                    iCoef2 <- pair.meanvcov[[iPattern]][2,iPair]

                    iValue <- iWeight * (tiX[iCoef1,,drop=FALSE] %*% ls.OmegaM1_dOmega_OmegaM1[[iPattern]][[iCoef2]] %*% iResidual) * scale.Omega[iId]

                    if(indiv){
                        info[iId,iCoef1,iCoef2] <- iValue
                        info[iId,iCoef2,iCoef1] <- iValue
                    }else{
                        info[iCoef1,iCoef2] <- info[iCoef1,iCoef2] + iValue
                        info[iCoef2,iCoef1] <- info[iCoef2,iCoef1] + iValue
                    }
                }

            }
        }
    }

    ## *** looping over covariance patterns
    if(!test.loopIndiv){
        U.pattern <- names(Upattern.ncluster)
        n.pattern <- length(U.pattern)
        vec.wXX.key <- as.double(precompute$wXX$key)
        
        info <- matrix(0, nrow = n.allcoef, ncol = n.allcoef,
                       dimnames = list(name.allcoef, name.allcoef))

        ## loop
        for(iPattern in U.pattern) { ## iPattern <- U.pattern[1]
            iName.varcoef <- name.varcoef[[iPattern]]
            iVec.vcov.key <- as.double(attr(pair.vcov[[iPattern]],"key"))
        
            ## **** mean,mean
            if(n.mucoef>0){
                info[name.mucoef,name.mucoef] <- info[name.mucoef,name.mucoef] + (attr(precompute$OmegaM1[[iPattern]],"vectorwise") %*% precompute$wXX$pattern[[iPattern]])[vec.wXX.key]                        }

            ## **** mean,var
            if(type.information == "observed" && n.allvarcoef>0 && n.mucoef>0){
                iTerm <- attr(precompute$OmegaM1_dOmega_OmegaM1[[iPattern]],"vectorwise") %*% t(attr(precompute$wXR[[iPattern]],"vectorwise"))
                info[iName.varcoef,name.mucoef] <- info[iName.varcoef,name.mucoef] + iTerm
                info[name.mucoef,iName.varcoef] <- info[name.mucoef,iName.varcoef] + t(iTerm)
            }

            ## **** var,var
            if(n.allvarcoef>0){                
                if(type.information == "expected"){
                    iValue <- 0.5 * Upattern.ncluster[iPattern] * attr(precompute$OmegaM1_dOmega_OmegaM1_dOmega[[iPattern]],"tr")
                }else if(type.information == "observed"){
                    iValue1 <- Upattern.ncluster[iPattern] * attr(precompute$d2Omega_dOmega_OmegaM1_dOmega[[iPattern]],"trModified")
                    iValue2 <- attr(precompute$d2Omega_dOmega_OmegaM1_dOmega[[iPattern]],"vectorwise") %*% attr(precompute$wRR[[iPattern]],"vectorwise")
                    iValue <- - 0.5 * (iValue1 + iValue2[,1])
                }
                info[iName.varcoef,iName.varcoef]  <- info[iName.varcoef,iName.varcoef] +  iValue[iVec.vcov.key]
            }
        }
    }

    ## *** add REML term
    if(REML && n.allvarcoef>0){
        browser()
        npair.vcov <- NCOL(attr(pair.vcov,"global"))


        term1 <- mapply(x = precompute$REML$wXOmegaM1X.M1, y = precompute$REML$wXOmegaM1dOmegaOmegaM1X, z = precompute$REML$wXOmegaM1.d2OmegaAndCo.OmegaM1X, function(x,y,z){
            xy <- x %*% y 
            sum(xy*xy) + sum(xy*xy)
        })
        tr(precompute$REML$wXOmegaM1X.M1 %*% precompute$REML$wXOmegaM1X.M1)
        sum(precompute$REML$wXOmegaM1X.M1 * precompute$REML$wXOmegaM1X.M1)
        

        REML.all <- as.double(REML.denomM1) %*% matrix(REML.numerator2.bis + REML.numerator2, nrow = prod(dim(REML.numerator2)[1:2]), ncol = dim(REML.numerator2)[3], byrow = FALSE)
        info[name.allvarcoef,name.allvarcoef] <- info[name.allvarcoef,name.allvarcoef] - 0.5 * REML.all[as.double(REML.key)]
        
    }


    ## ** export

    if(robust){
        if(REML){
            effects2 <- "mean"
            attr(effects2,"original.names") <- attr(effects,"original.names")
            attr(effects2,"reparametrize.names") <- attr(effects,"reparametrize.names")
        }else{
            effects2 <- effects
        }
        if(is.null(weights)){
            weights <- rep(1, length(pattern))
        }
        if(is.null(scale.Omega)){
            scale.Omega <- rep(1, length(pattern))
        }
        attr.info <- info
        attr.bread <- crossprod(.score(X = X, residuals = residuals, precision = precision, dOmega = dOmega,
                                       weights = weights, scale.Omega = scale.Omega, pattern = pattern, 
                                       index.cluster = index.cluster, name.allcoef = name.allcoef, indiv = TRUE, REML = REML, effects = effects2,
                                       precompute = precompute) )

        
        
        if(any(c("mean","variance","correlation") %in% effects2 == FALSE)){
            keep.cols <- intersect(names(which(rowSums(abs(attr.bread))!=0)),names(which(rowSums(abs(attr.bread))!=0)))
            info <- NA*attr.info

            attr.infoM1 <- solve(attr.info)
            info[keep.cols,keep.cols] <- solve(attr.infoM1[keep.cols,keep.cols] %*% attr.bread[keep.cols,keep.cols,drop=FALSE] %*% attr.infoM1[keep.cols,keep.cols])
        }else{
            attr.infoM1 <- solve(attr.info)
            info <- solve(attr.infoM1 %*% attr.bread %*% attr.infoM1)
        }
    }
    return(info)
}

##----------------------------------------------------------------------
### information.R ends here

