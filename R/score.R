### score.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:59) 
## Version: 
## Last-Updated: aug  2 2024 (15:58) 
##           By: Brice Ozenne
##     Update #: 928
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * score.lmm (documentation)
##' @title Extract The Score From a Linear Mixed Model
##' @description Extract or compute the first derivative of the log-likelihood of a linear mixed model.
##' 
##' @param x a \code{lmm} object.
##' @param effects [character] Should the score relative to all coefficients be output (\code{"all"}),
##' @param indiv [logical] Should the contribution of each cluster to the score be output? Otherwise output the sum of all clusters of the derivatives.
##' or only coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only coefficients relative to the variance and correlation structure (\code{"variance"} or \code{"correlation"}).
##' @param newdata [data.frame] dataset relative to which the score should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the score. Only relevant if differs from the fitted values.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details For details about the arguments \bold{transform.sigma}, \bold{transform.k}, \bold{transform.rho}, see the documentation of the \link[LMMstar]{coef.lmm} function.
##'
##' @return
##' When argument indiv is \code{FALSE}, a vector with the value of the score relative to each coefficient.
##' When argument indiv is \code{TRUE}, a matrix with the value of the score relative to each coefficient (in columns) and each cluster (in rows).
##'
##' @keywords methods

## * score.lmm (code)
##' @export
score.lmm <- function(x, effects = "mean", indiv = FALSE, newdata = NULL, p = NULL,
                      transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

    options <- LMMstar.options()

    ## ** normalize user input

    ## *** dots
    dots <- list(...)
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

    ## *** p and transform
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

    ## ** extract or recompute score
    if(is.null(newdata) && is.null(p) && (indiv == FALSE) && test.notransform){
        keep.name <- stats::setNames(names(coef(x, effects = effects, transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE)),
                                     names(coef(x, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)))    

        design <- x$design ## useful in case of NA
        out <- x$score[keep.name]
        if(transform.names){
            names(out) <- names(keep.name)
        }

    }else{
        if(!is.null(newdata)){
            design <- stats::model.matrix(x, newdata = newdata, effects = "all", simplify = FALSE)
        }else{
            design <- x$design
        }
        out <- .moments.lmm(value = theta, design = design, time = x$time, method.fit = x$args$method.fit, type.information = x$args$type.information,
                            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                            logLik = FALSE, score = TRUE, information = FALSE, vcov = FALSE, df = FALSE, indiv = indiv, effects = effects, robust = FALSE,
                            trace = FALSE, precompute.moments = !is.null(x$design$precompute.XX), transform.names = transform.names)$score
    }

    ## ** name and restaure NAs
    if(indiv){

        if(!is.numeric(x$cluster$levels)){
            rownames(out) <- x$cluster$levels[match(1:NROW(out),x$cluster$index)]
        } 
        out <- restaureNA(out, index.na = x$index.na,
                          level = "cluster", cluster = x$cluster)        
        
    }

    ## ** export
    return(out)
}

## * score.mlmm (documentation)
##' @title Extract The Score From Multiple Linear Mixed Models
##' @description Extract or compute the first derivative of the log-likelihood of each linear mixed model.
##'
##' @param x a \code{mlmm} object.
##' @param effects [character] By default will output the estimates relative to the hypotheses being tested (\code{"contrast"}).
##' But can also output all model coefficients (\code{"all"}),
##' or only coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only coefficients relative to the variance structure (\code{"variance"}),
##' or only coefficients relative to the correlation structure (\code{"correlation"}).
##' @param indiv [logical] Should the contribution of each cluster to the score be output? Otherwise output the sum of all clusters of the derivatives.
##' @param p [list of numeric vector] list of model coefficients to be used. Only relevant if differs from the fitted values.
##' @param newdata [NULL] Not used. For compatibility with the generic method.
##' @param ordering [character] should the output be ordered by type of parameter (\code{parameter}) or by model (\code{by}).
##' @param simplify [logical] should the score be combined across models into a single vector (\code{indiv=FALSE}) or matrix (\code{indiv=TRUE})?
##' @param ... passed to \code{score.lmm}.
##' 
##' @return
##' When argument indiv is \code{FALSE}, a vector with the value of the score relative to each coefficient.
##' When argument indiv is \code{TRUE}, a matrix with the value of the score relative to each coefficient (in columns) and each cluster (in rows).
##' When \code{effects} differs from \code{"contrast"} and \code{simplify=FALSE}, it will store the score in a list with an element relative to each parameter or model (argument \code{ordering}).

## * score.mlmm (code)
##' @export
score.mlmm <- function(x, effects = "contrast", indiv = FALSE, p = NULL, newdata = NULL, ordering = "by", simplify = TRUE, ...){

    level.by <- names(x$model)

    ## ** normalize use input

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character vector")
    }
    valid.effects <- c("contrast","mean","fixed","variance","correlation","all")
    if(any(effects %in% valid.effects == FALSE)){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }    
    if("contrast" %in% effects){
        if(length(effects)>1){
            stop("Argument \'effects\' must have length 1 when containing the element \'effects\'. \n")
        }
        if(!is.null(x$univariate) & all(x$univariate$type=="mu")){
            effects2 <- "mean"
        }else{
            effects2 <- "all"
        }
    }else{
        effects2 <- effects
    }
    if("all" %in% effects && length(effects)>1){
        stop("Argument \'effects\' must have length 1 when containing the element \'all\'. \n")
    }

    ## *** p
    if(!is.null(p)){
        if(!is.list(p)){
            stop("Argument \'p\' should either be NULL or a list. \n")
        }
        if(is.null(names(p))){
            stop("Argument \'p\' should either be NULL or a named list. \n")
        }
        if(any(names(p) %in% level.by == FALSE)){
            stop("Incorrect names for argument \'p\': \"",paste(setdiff(names(p),level.by), collapse = "\", \""),"\". \n", 
                 "Should be among \"",paste(level.by, collapse = "\", \""),"\". \n")
        }
    }

    ## *** newdata
    if(!is.null(newdata)){
        message("Argument \'newdata\' is being ignored. \n")
    }

    ## *** ordering    
    ordering <- match.arg(ordering, c("by","parameter"))

    ## *** simplify
    if(!is.numeric(simplify) && !is.logical(simplify)){
        stop("Argument \'simplify\' must be numeric or logical. \n")
    }
    if(length(simplify)!=1){
        stop("Argument \'simplify\' must have length 1. \n")
    }
    if(simplify %in% c(0,1) == FALSE){
        stop("Argument \'simplify\' must be TRUE/1 or FALSE/0. \n")
    }

    ## ** extract score
    if(indiv){
        cluster <- x$object$cluster
        n.cluster <- length(cluster)
    }

    if(all(effects=="contrast")){
        ls.contrast <- stats::coef(x, type = "ls.contrast") ## extract linear combination from each model
        contrast <- coef(x, type = "contrast") ## contrast combinations across models
        ls.score <- lapply(level.by, function(iBy){
            score(x$model[[iBy]], effects = effects2, indiv = indiv, p = p[[iBy]], transform.names = FALSE, ...)            
        })

        ls.out <- mapply(iScore = ls.score, iC = ls.contrast, FUN = function(iScore,iC){ ## iScore <- ls.score[[1]] ; iC <- ls.contrast[[1]]
            if(indiv){
                iOut <- matrix(0, nrow = n.cluster, ncol = NROW(iC), dimnames = list(cluster, rownames(iC)))
                iOut[rownames(iScore),] <- iScore %*% t(iC[,colnames(iScore), drop = FALSE])
                return(iOut)
            }else{
                iOut <- iScore %*% t(iC[,names(iScore), drop = FALSE])
                return(iOut[1,])
            }
        }, SIMPLIFY = FALSE)

        if(indiv){
            out <- (do.call("cbind",ls.out) %*% t(contrast))
        }else{
            out <- (do.call("c",ls.out) %*% t(contrast))
        }
        attr(out,"message") <- unique(unlist(lapply(ls.out,attr,"message")))

    }else{

        ls.out <- lapply(level.by, function(iBy){
            iO <- x$model[[iBy]]

            iScore <- score(iO, effects = effects, indiv = indiv, p = p[[iBy]], transform.names = FALSE, ...)
            if(indiv){
                if(simplify){
                    iOut <- matrix(0, nrow = n.cluster, ncol = NCOL(iScore), dimnames = list(cluster, paste0(iBy,": ", colnames(iScore))))
                    attr(iOut,"name") <- colnames(iScore)
                    attr(iOut,"by") <- rep(iBy, NCOL(iScore))                
                }else{
                    iOut <- matrix(0, nrow = n.cluster, ncol = NCOL(iScore), dimnames = list(cluster, colnames(iScore)))
                }
                iOut[rownames(iScore),] <- iScore
            }else{
                iOut <- iScore                
                if(simplify){
                    names(iOut) <- paste0(iBy,": ", names(iScore))
                    attr(iOut,"name") <- names(iScore)
                    attr(iOut,"by") <- rep(iBy, length(iScore))
                }
            }
            attr(iOut,"message") <- attr(iScore,"message")
            return(iOut)
        })
        names(ls.out) <- level.by
    }

    ## ** re-order
    if(all(effects=="contrast")){
        name.order <- names(coef(x, effects = "contrast", ordering = ordering))        
        if(indiv){
            out <- out[,name.order,drop=FALSE]
        }else{
            out <- out[,name.order]
        }
    }else if(simplify){

        if(indiv){
            out <- do.call("cbind",unname(ls.out))
        }else{
            out <- do.call("c",unname(ls.out))
        }
        attr(out,"original.name") <- unname(unlist(lapply(ls.out,attr,"name")))
        attr(out,"by") <- unname(unlist(lapply(ls.out,attr,"by")))
        attr(out,"message") <- unname(unique(unlist(lapply(ls.out,attr,"message"))))

        if(ordering == "parameter"){
            reorder <- order(factor(attr(out,"original.name"),levels = unique(attr(out,"original.name"))))
            out.save <- out
            if(indiv){
                out <- out.save[,reorder,drop=FALSE]
            }else{
                out <- out.save[reorder]
            }
            attr(out, "original.name") <- attr(out.save, "original.name")[reorder]
            attr(out, "by") <- attr(out.save, "by")[reorder]
            attr(out, "message") <- attr(out.save, "message")
        }

    }else if(ordering == "by"){

        out <- ls.out

    }else if(ordering == "parameter"){

        ## unique parameters
        Uname <- unique(unlist(lapply(ls.out,function(iOut){
            if(indiv){colnames(iOut)}else{names(iOut)}
        })))

        ## combine score
        if(indiv){
            out <- stats::setNames(lapply(Uname, function(iName){ ## iName <- "X1"
                iOut <- do.call(cbind,lapply(ls.out, function(iOut){ iOut[,iName]}))
                colnames(iOut) <- names(ls.out)
                attr(iOut,"message") <- unname(unique(unlist(lapply(ls.out, attr, "message"))))
                return(iOut)
            }), Uname)
        }else{
            out <- stats::setNames(lapply(Uname, function(iName){ ## iName <- "X1"
                iOut <- sapply(ls.out,"[",iName)
                names(iOut) <- names(ls.out)
                attr(iOut,"message") <- unname(unique(unlist(lapply(ls.out, attr, "message"))))
                return(iOut)
            }), Uname)
        }        
    }

    ## ** export
    return(out)
}

## * .score
.score <- function(X, residuals, precision, dOmega, weights, 
                   pattern, index.cluster, name.allcoef, indiv, REML, effects, precompute){

    ## ** extract information
    n.cluster <- length(pattern)
    name.mucoef <- colnames(X)
    n.mucoef <- length(name.mucoef)
    name.varcoef <- unique(unlist(lapply(dOmega,names)))
    n.varcoef <- length(name.varcoef)
    U.pattern <- names(dOmega)

    ## ** prepare output
    compute.indiv <- indiv || is.null(precompute$weights) || is.null(precompute$XR) || is.null(precompute$RR)
    name.effects <- attr(effects,"original.names")
    n.effects <- length(name.effects)
    if(compute.indiv){
        Score <- matrix(NA, nrow = n.cluster, ncol = n.effects,
                        dimnames = list(NULL, name.effects))
    }else{
        Score <- stats::setNames(rep(0, n.effects), name.effects)
    }

    ## ** restrict to relevant parameters
    test.mean <- (n.mucoef>0) && ("mean" %in% effects)
    test.vcov <- (n.varcoef>0) && (("variance" %in% effects) || ("correlation" %in% effects))

    if(test.vcov){
        if(REML && indiv){
            message <- "approximate individual REML score"
        }else{
            message <- NULL
        }

        if(("variance" %in% effects == FALSE) || ("correlation" %in% effects == FALSE)){## restrict to requested coefficients
            name.varcoef <- intersect(name.varcoef,name.effects)
            n.varcoef <- length(name.varcoef)
            precompute$Omega$tr.OmegaM1.dOmega <- lapply(precompute$Omega$tr.OmegaM1.dOmega, function(iO){iO[intersect(names(iO),name.effects)]})
            precompute$Omega$OmegaM1.dOmega.OmegaM1 <- lapply(precompute$Omega$OmegaM1.dOmega.OmegaM1, function(iO){iO[,intersect(colnames(iO),name.effects),drop=FALSE]})
            test.vcov <- n.varcoef>0
        }
    }else{
        message <- NULL
    }

    ## ** prepare REML term
    if(REML && test.vcov){
        if(is.null(precompute$REML)){
            X.OmegaM1.X <- matrix(0, nrow = n.mucoef, ncol = n.mucoef)
            REML.num <- stats::setNames(replicate(n.varcoef, matrix(0, nrow = n.mucoef, ncol = n.mucoef), simplify = FALSE), name.varcoef)
            for(iId in 1:n.cluster){ ## iId <- 1
                iX <- X[index.cluster[[iId]],,drop=FALSE]
                iOmegaM1 <- precision[[pattern[iId]]]
                iOmegaM1.dOmega.OmegaM1 <- precompute$Omega$OmegaM1.dOmega.OmegaM1[[pattern[iId]]]
                iWeights <- weights[iId]
                X.OmegaM1.X <- X.OmegaM1.X + iWeights * t(iX) %*% iOmegaM1 %*% iX
                for(iParam in intersect(names(iOmegaM1.dOmega.OmegaM1), name.varcoef)){ ## intersect to handle when argument effects is only "variance" or "correlation"
                    REML.num[[iParam]] <- REML.num[[iParam]] + iWeights * t(iX) %*% iOmegaM1.dOmega.OmegaM1[[iParam]] %*% iX
                }
            }
            REML.denom <- solve(X.OmegaM1.X)
        }else{
            REML.num <- precompute$REML$X.OmegaM1.dOmega.OmegaM1.X
            REML.denom <- precompute$REML$X.OmegaM1.X_M1
        }
    }

    ## ** compute score
    if(compute.indiv){
        ## *** looping over individuals
        if(REML && test.vcov){
            OmegaM1.dOmega.OmegaM1 <- lapply(precompute$Omega$OmegaM1.dOmega.OmegaM1, function(iO){
                array(iO, dim = c(sqrt(NROW(iO)),sqrt(NROW(iO)),NCOL(iO)), dimnames = list(NULL,NULL,colnames(iO)))
            })
        }
                        
        for(iId in 1:n.cluster){ ## iId <- 7
            iIndex <- index.cluster[[iId]]
            iPattern <- pattern[iId]
            iWeights <- weights[iId]
            iResidual <- residuals[iIndex,,drop=FALSE]
            iX <- t(X[iIndex,,drop=FALSE])
            iOmegaM1 <- precision[[iPattern]]

            if(test.mean){
                Score[iId,name.mucoef] <- weights[iId] * (iX %*% iOmegaM1 %*% iResidual)
            }
            if(test.vcov){
                iVarcoef <- names(precompute$Omega$tr.OmegaM1.dOmega[[pattern[iId]]])
                Score[iId,iVarcoef] <- 0.5 * iWeights * (-precompute$Omega$tr.OmegaM1.dOmega[[iPattern]] + as.vector(tcrossprod(iResidual)) %*% precompute$Omega$OmegaM1.dOmega.OmegaM1[[iPattern]])
                ## second term is the same as
                ## iResidual[,1] %*% matrix(precompute$Omega$OmegaM1.dOmega.OmegaM1[[pattern[iId]]][,2], nrow = length(iIndex), ncol = length(iIndex)) %*% iResidual 

                if(REML){
                    ## APPROXIMATION: the REML score w.r.t. variance parameter is not linear in the individual contribution tr(\sum_i a_i / \sum_i b_i)/2
                    ## first order expansion:  tr(a_j / \sum_i b_i)/2

                    for(iParam in intersect(dimnames(OmegaM1.dOmega.OmegaM1[[iPattern]])[[3]], name.varcoef)){ ## intersect to handle when argument effects is only "variance" or "correlation"
                        iREML.num <- iWeights * iX %*% OmegaM1.dOmega.OmegaM1[[iPattern]][,,iParam] %*% t(iX)
                        ## iREML.num <- iWeights * iX %*% iOmegaM1 %*% dOmega[[iPattern]][[iParam]] %*% iOmegaM1 %*% t(iX)
                        ## iREML.denom <- iWeights * iX %*% dOmega[[iPattern]][[iParam]] %*% t(iX)
                        Score[iId,iParam] <-  Score[iId,iParam] + 0.5 * sum(iREML.num * REML.denom)
                    }
                }

            }
        }

        if(!indiv){
            Score <- colSums(Score)
        }

    }else{
        ## *** looping over covariance patterns
        for (iPattern in U.pattern) { ## iPattern <- U.pattern[1]
            iOmegaM1 <- precision[[iPattern]]
            iTime2 <- length(iOmegaM1)

            if(test.mean){
                Score[name.mucoef] <- Score[name.mucoef] + attr(iOmegaM1,"vectorize") %*% precompute$XR[[iPattern]]
            }

            if(test.vcov){
                iVarcoef <- names(precompute$Omega$tr.OmegaM1.dOmega[[iPattern]])
                ## precompute$RR has already been weigthed
                Score[iVarcoef] <- Score[iVarcoef] + 0.5 * (-precompute$weights[iPattern] * precompute$Omega$tr.OmegaM1.dOmega[[iPattern]] + precompute$RR[[iPattern]] %*% precompute$Omega$OmegaM1.dOmega.OmegaM1[[iPattern]])
            }
        }

        ## *** REML contribution
        if(REML && test.vcov){
            Score[name.varcoef] <-  Score[name.varcoef] + 0.5 * sapply(REML.num[name.varcoef], function(iO){sum(iO * REML.denom)})
        }

    }

    ## ** export
    attr(Score,"message") <- message
    return(Score)
}


##----------------------------------------------------------------------
### score.R ends here
