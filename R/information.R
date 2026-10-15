### information.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (22:13) 
## Version: 
## Last-Updated: okt 15 2026 (14:21) 
##           By: Brice Ozenne
##     Update #: 1411
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
    x.transform <- transform(x)
    init <- .init_transform(p = p, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = x.transform$sigma, x.transform.k = x.transform$k, x.transform.rho = x.transform$rho,
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
                            updateGLS = FALSE, logLik = FALSE, score = FALSE, information = TRUE, vcov = FALSE, df = FALSE, indiv = indiv, effects = effects, robust = FALSE,
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

    ## ** restrict to relevant parameters
    test.mean <- (n.mucoef>0) && ("mean" %in% effects)
    test.vcov <- (n.vcovcoef>0) && (("variance" %in% effects) || ("correlation" %in% effects))
    if(test.vcov && REML && indiv){
        message <- "approximate individual REML information"
    }else{
        message <- NULL
    }
    

    ## ** prepare output
    name.effects <- NULL
    if(test.mean){
        name.effects <- c(name.effects,name.mucoef)
    }
    if(test.vcov){
        name.effects <- c(name.effects,name.vcovcoef)
        if(test.mean){
            pair.vcov$key1 <- pair.vcov$key1.all
            pair.vcov$key2 <- pair.vcov$key2.all
        }else{
            pair.vcov$key1 <- pair.vcov$key1.vcov
            pair.vcov$key2 <- pair.vcov$key2.vcov
        }
    }
    n.effects <- length(name.effects)
    if(indiv){
        info <- array(0, dim = c(n.cluster, n.effects, n.effects),
                      dimnames = list(NULL, name.effects, name.effects))
    }else{
        info <- matrix(0, nrow = n.effects, ncol = n.effects,
                       dimnames = list(name.effects, name.effects)
                       )
    }    

    ## ** compute information
    if(REML && test.vcov){
        n.time <- sqrt(length(precompute$REML$X.OmegaM1.X_M1))
        M.REML.denom <- matrix(precompute$REML$X.OmegaM1.X_M1, ncol = n.time, nrow = n.time)
        REML.num1 <- do.call(cbind,apply(precompute$REML$X.dOmegaM1.X, MARGIN = 2, function(iVec){as.vector(M.REML.denom %*% matrix(iVec, ncol = n.time, nrow = n.time))}, simplify = FALSE))
    }

    if(indiv){ ## *** looping over individuals

        if(test.vcov){
            iStore.vcov <- matrix(0, nrow = dim(info)[2], ncol = dim(info)[3], dimnames = dimnames(info)[2:3])
        }
        
        for(iId in 1:n.cluster){ ## iId <- 1
            
            iIndex <- index.cluster[[iId]]
            iPattern <- pattern[iId,"all"]
            iWeight <- weights$likelihood[iId]
            iReplicate <- weights$Omega[iId]
            if(test.mean || (REML && test.vcov)){
                iX <- t(X[iIndex,,drop=FALSE])
            }
            if(test.vcov && (REML || (type.information == "observed" && test.mean))){
                idOmegaM1 <- precompute$Omega$dOmegaM1[[iPattern]]
            }
            if(test.vcov && (REML || (type.information == "observed"))){
                id2OmegaM1 <- precompute$Omega$d2OmegaM1[[iPattern]]
            }
            
            if(test.mean){ ## mean mean
                iOmegaM1 <- precompute$Omega$OmegaM1[[iPattern]]
                info[iId,name.mucoef,name.mucoef] <- (iWeight * iReplicate) * (iX %*% iOmegaM1 %*% t(iX))
            }
            
            if(test.vcov){ ## var var
                iStore.vcov[] <- 0

                ## position in the information matrix
                iPattern.vcov <- pattern[iId,"vcov"]
                iKey <- c(stats::setNames(pair.vcov$key1,pair.vcov$name)[which(pair.vcov[[iPattern.vcov]])],
                          stats::setNames(pair.vcov$key2,pair.vcov$name)[intersect(which(pair.vcov[[iPattern.vcov]]),which(!is.na(pair.vcov$key2)))])
                ## any(duplicated(iKey))

                ## trace term
                id2Trace <- precompute$Omega$d2Trace[[iPattern]]
                iStore.vcov[iKey] <- - (0.5 * iWeight * id2Trace)[names(iKey)]

                ## residual term
                if(type.information == "observed"){
                    iResidual <- residuals[iIndex,,drop=FALSE]
                    iStore.vcov[iKey]  <- iStore.vcov[iKey] - 0.5 * (as.vector(tcrossprod(iResidual)) %*% id2OmegaM1)[1,names(iKey)]
                }

                ## store
                info[iId,,] <- iStore.vcov
            }
            
            if(type.information == "observed" && test.mean && test.vcov){ ## mean var
                
                iValue <- (iWeight * iReplicate) * kronecker(t(iResidual),iX) %*% idOmegaM1
                ## vec(AMB^T) = (B %x% A) vec(M)
                ## iX %*% matrix(idOmegaM1[,1], NROW(iOmegaM1), NCOL(iOmegaM1)) %*% iResidual - kronecker(t(iResidual),iX) %*% idOmegaM1[,1]
                info[iId,name.mucoef,colnames(idOmegaM1)] <- iValue
                info[iId,colnames(idOmegaM1),name.mucoef] <- t(iValue)

            }

            if(REML && test.vcov){
                ## APPROXIMATION: the REML score w.r.t. variance parameter is not linear in the individual contribution
                REML.denom <- precompute$REML$X.OmegaM1.X_M1
                Key <- c(stats::setNames(pair.vcov$key1,pair.vcov$name),
                         stats::setNames(pair.vcov$key2,pair.vcov$name)[which(!is.na(pair.vcov$key2))])
                
                iREML.trace1a <- colSums((kronecker(iX,iX) %*% idOmegaM1[,pair.vcov$param1,drop=FALSE]) * REML.num1[,pair.vcov$param2,drop=FALSE])
                iREML.trace1b <- colSums(REML.num1[,pair.vcov$param1,drop=FALSE] * (kronecker(iX,iX) %*% idOmegaM1[,pair.vcov$param2,drop=FALSE]))
                iREML.trace2 <- colSums(sweep(kronecker(iX,iX) %*% id2OmegaM1, FUN = "*", MARGIN = 1, STATS = REML.denom))
                
                iStore.vcov[] <- 0
                iStore.vcov[iKey] <- - 0.5 * stats::setNames(iREML.trace1a + 0.5*iREML.trace1b + 0.5*iREML.trace2, pair.vcov$name)[names(Key)]
                info[iId,,] <- info[iId,,] + iStore.vcov

            }
        }

    }else{ ## *** looping over covariance patterns

        for (iG in 1:NROW(precompute$weights)) { ## iG <- 1

            iPattern <- precompute$weights[iG,"pattern"]

            if(test.mean){ ## mean mean
                iOmegaM1 <- precompute$Omega$OmegaM1[[iPattern]]
                iXX <- precompute$XX$pattern[[iPattern]]## has already been weigthed
                info[name.mucoef,name.mucoef] <- info[name.mucoef,name.mucoef] + (as.vector(iOmegaM1) %*% iXX)[as.double(precompute$XX$key)]                
            }

            if(test.vcov){ ## var var
                ## position in the information matrix
                iPattern.vcov <- precompute$weights[iG,"vcov"]
                iKey <- c(stats::setNames(pair.vcov$key1,pair.vcov$name)[which(pair.vcov[[iPattern.vcov]])],
                          stats::setNames(pair.vcov$key2,pair.vcov$name)[intersect(which(pair.vcov[[iPattern.vcov]]),which(!is.na(pair.vcov$key2)))])
                ## any(duplicated(iKey))

                ## trace term
                iWeight <- precompute$weights[iG,"sum_w"]
                id2Trace <- precompute$Omega$d2Trace[[iPattern]]
                info[iKey]  <- info[iKey] - (0.5 * iWeight * id2Trace)[names(iKey)]

                ## residual term
                if(type.information == "observed"){
                    iRR <- precompute$RR[[iPattern]]
                    id2OmegaM1 <- precompute$Omega$d2OmegaM1[[iPattern]]
                    info[iKey]  <- info[iKey] - 0.5 * (iRR %*% id2OmegaM1)[1,names(iKey)]
                }
            }

            if(type.information == "observed" && test.mean && test.vcov){ ## mean var
                iXR <- precompute$XR[[iPattern]]
                idOmegaM1 <- precompute$Omega$dOmegaM1[[iPattern]]
                info[colnames(idOmegaM1),name.mucoef] <- info[colnames(idOmegaM1),name.mucoef] + t(idOmegaM1) %*% iXR
                info[name.mucoef,colnames(idOmegaM1)] <- info[name.mucoef,colnames(idOmegaM1)] + t(iXR) %*% idOmegaM1
            }
        }

        if(REML && test.vcov){ ## REML

            REML.denom <- precompute$REML$X.OmegaM1.X_M1
            REML.num2 <- precompute$REML$X.d2OmegaM1.X[,pair.vcov$name,drop=FALSE]
            Key <- c(stats::setNames(pair.vcov$key1,pair.vcov$name),
                     stats::setNames(pair.vcov$key2,pair.vcov$name)[which(!is.na(pair.vcov$key2))])

            REML.trace1 <- colSums(REML.num1[,pair.vcov$param1,drop=FALSE] * REML.num1[,pair.vcov$param2,drop=FALSE])
            REML.trace2 <- colSums(sweep(REML.num2, FUN = "*", MARGIN = 1, STATS = REML.denom))
            info[iKey] <- info[iKey] - 0.5 * stats::setNames(REML.trace1 + REML.trace2, pair.vcov$name)[names(Key)]
            
            ## tr(M.REML.denom %*% matrix(precompute$REML$X.dOmegaM1.X[,1],n.time,n.time) %*% M.REML.denom %*% matrix(precompute$REML$X.dOmegaM1.X[,1],n.time,n.time))
            ## colSums(REML.num1[,pair.vcov$param1,drop=FALSE] * REML.num1[,pair.vcov$param2,drop=FALSE])[1]
        }
    }

    ## ** export
    attr(info,"message") <- message
    return(info)
}

##----------------------------------------------------------------------
### information.R ends here

