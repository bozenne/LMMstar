### score.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:59) 
## Version: 
## Last-Updated: jul 16 2024 (17:20) 
##           By: Brice Ozenne
##     Update #: 765
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
        effects <- options$effects
    }else if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","fixed","variance","correlation"), several.ok = TRUE)
    effects[effects== "fixed"] <- "mean"

    ## *** indiv
    if(identical(indiv,"REML2ML") && x$args$method.fit == "REML"){
        indiv <- TRUE
        attr(indiv,"REML2ML") <- TRUE
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
        test.precompute <- !is.null(x$design$precompute.XX) && !indiv

        if(!is.null(newdata)){
            design <- stats::model.matrix(x, newdata = newdata, effects = "all", simplify = FALSE)
        }else{
            design <- x$design
        }        
        out <- .moments.lmm(value = theta, design = design, time = x$time, method.fit = x$args$method.fit,
                            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                            logLik = FALSE, score = TRUE, information = FALSE, vcov = FALSE, df = FALSE, indiv = indiv, effects = effects,
                            trace = FALSE, precompute.moments = test.precompute, transform.names = transform.names)$score
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

    ## *** indiv
    if(identical(indiv,"REML2ML") && x$object$method.fit == "REML"){
        indiv <- TRUE
        attr(indiv,"REML2ML") <- TRUE
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
.score <- function(X, residuals, precision, dOmega,
                   Upattern.ncluster, weights, scale.Omega,
                   pattern, index.cluster, name.allcoef,
                   indiv, REML, effects,
                   precompute){


    ## ** extract information
    test.loopIndiv <- indiv || is.null(precompute)
    n.obs <- length(index.cluster)
    n.cluster <- length(pattern)
    name.mucoef <- colnames(X)
    n.mucoef <- length(name.mucoef)
    name.varcoef <- lapply(dOmega,names)
    n.varcoef <- lengths(name.varcoef)
    name.allvarcoef <- unique(unlist(name.varcoef))
    U.pattern <- names(dOmega)
    n.pattern <- length(U.pattern)

    ## ** prepare output
    name.effects <- attr(effects,"original.names")
    n.effects <- length(name.effects)
    if(test.loopIndiv){
        Score <- matrix(0, nrow = n.cluster, ncol = n.effects,
                        dimnames = list(NULL, name.effects))
    }else{
        Score <- stats::setNames(rep(0, n.effects), name.effects)
    }
    if(any(is.na(attr(precision, "logdet")))){ ## non positive definite residual variance covariance
        return(Score*NA)
    }

    ## restrict to relevant parameters
    if(("variance" %in% effects == FALSE) && ("correlation" %in% effects == FALSE)){ ## compute score only for mean parameters
        test.vcov <- FALSE
        test.mean <- n.mucoef>0
        message <- NULL
    }else{
        if(REML && indiv){
            if(identical(attr(indiv,"REML2ML"),TRUE)){
                REML <- FALSE
                message <- "score:REML2ML"
            }else{
                stop("Not possible to compute individual REML score for variance and/or correlation coefficients.\n")
                ## "Consider setting the attribute score:REML2ML to TRUE in argument \'indiv\' to use individual ML score with REML estimates as an approximation."
            }
        }else{
            message <- NULL
        }
        if(("variance" %in% effects == FALSE) || ("correlation" %in% effects == FALSE)){
            name.varcoef <- stats::setNames(lapply(U.pattern,function(iPattern){intersect(name.effects,name.varcoef[[iPattern]])}),
                                            U.pattern)
            n.varcoef <- lengths(name.varcoef)
            name.allvarcoef <- unique(unlist(name.varcoef))
        }
        if("mean" %in% effects == FALSE){ ## compute score only for variance and/or correlation parameters
            if(REML && indiv){
                stop("Not possible to compute individual score for variance and/or correlation coefficients when using REML.\n")
            }

            test.vcov <- any(n.varcoef>0)
            test.mean <- FALSE

        }else{ ## compute score all parameters
     
            test.vcov <- any(n.varcoef>0)
            test.mean <- n.mucoef>0
        }
    }

    if(test.vcov){
        ls.OmegaM1_dOmega_OmegaM1 <- attr(dOmega,"ls.OmegaM1_dOmega_OmegaM1")
        OmegaM1_dOmega_OmegaM1 <- attr(dOmega,"OmegaM1_dOmega_OmegaM1")
                
        if(REML){
            REML.num <- array(0, dim = c(n.mucoef, n.mucoef, length(name.allvarcoef)), dimnames = list(name.mucoef,name.mucoef,name.allvarcoef))
            REML.denom <- matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef, name.mucoef))
        }
    }

    ## ** compute score
    ## *** looping over individuals
    if(test.loopIndiv){

        if(test.vcov){ ## precompute
            trOmegaM1_dOmega <- stats::setNames(vector(mode = "list", length = n.pattern), U.pattern)
            for(iPattern in 1:n.pattern){ ## iPattern <- 4
                trOmegaM1_dOmega[[iPattern]]  <- stats::setNames(lapply(name.varcoef[[iPattern]], function(iVarcoef){tr(precision[[iPattern]] %*% dOmega[[iPattern]][[iVarcoef]])}), name.varcoef[[iPattern]])
            }
        }
        
        ## loop
        for(iId in 1:n.cluster){ ## iId <- 7
            iPattern <- pattern[iId]
            iIndex <- index.cluster[[iId]]
            iWeight <- weights[iId]
            iOmegaM1 <- precision[[pattern[iId]]] * scale.Omega[iId]

            iResidual <- residuals[iIndex,,drop=FALSE]
            iX <- X[iIndex,,drop=FALSE]
            tiX <- t(iX)

            if(test.mean){
                Score[iId,name.mucoef] <- iWeight * (tiX %*% iOmegaM1 %*% iResidual)
            }

            if(test.vcov){
                if(REML){
                    REML.denom <- REML.denom + iWeight * (tiX %*% iOmegaM1 %*% iX)
                }

                for(iVarcoef in name.varcoef[[iPattern]]){ ## iVarcoef <- name.varcoef[1]
                    Score[iId,iVarcoef] <- -0.5 * iWeight * trOmegaM1_dOmega[[iPattern]][[iVarcoef]] + 0.5 * iWeight * (t(iResidual) %*% ls.OmegaM1_dOmega_OmegaM1[[iPattern]][[iVarcoef]] %*% iResidual) * scale.Omega[iId]

                    if(REML){
                        REML.num[,,iVarcoef] <- REML.num[,,iVarcoef] + iWeight * (tiX %*% ls.OmegaM1_dOmega_OmegaM1[[iPattern]][[iVarcoef]] %*% iX) * scale.Omega[iId]
                    }
                }
            }
        }

        if(!indiv){
            Score <- colSums(Score)
        }
    }

    ## *** looping over covariance patterns
    if(!test.loopIndiv){

        ## loop
        for (iPattern in U.pattern) { ## iPattern <- U.pattern[1]
            iName.varcoef <- name.varcoef[[iPattern]]
            iOmegaM1 <- precision[[iPattern]]
            iTime2 <- length(iOmegaM1)

            if(test.mean){
                ## X %*% iOmega^-1 %*% residual
                Score[name.mucoef] <- Score[name.mucoef] + as.double(iOmegaM1) %*% matrix(unlist(precompute$XR[[iPattern]]), nrow = iTime2, ncol = n.mucoef, byrow = FALSE)
                ## Score[name.mucoef] <- Score[name.mucoef] + apply(precompute$XR[[iPattern]], MARGIN = 3, FUN = function(iM){sum(iM * iOmegaM1)})
            }

            if(test.vcov){
                iTrace <- sapply(dOmega[[iPattern]], FUN = function(iO){sum(iO*precision[[iPattern]])})
                ## iOmega^-1 dOmega iOmega^-1
                Score[iName.varcoef] <- Score[iName.varcoef] - 0.5 * Upattern.ncluster[iPattern] * iTrace + 0.5 * as.double(precompute$RR[[iPattern]]) %*% OmegaM1_dOmega_OmegaM1[[iPattern]]
                
                if(REML){
                    iX <- precompute$XX$pattern[[iPattern]]
                    iDouble2Mat <- as.vector(precompute$XX$key)
                    ## denominator
                    if(is.null(precompute$X.OmegaM1.X)){
                        REML.denom <- REML.denom + (as.double(iOmegaM1) %*% iX)[iDouble2Mat]
                    }else{
                        REML.denom <- REML.denom + precompute$X.OmegaM1.X[[iPattern]][iDouble2Mat]
                    }
                    ## numerator
                    iX_OmegaM1_dOmega_OmegaM1_X <- t(iX) %*% OmegaM1_dOmega_OmegaM1[[iPattern]]
                    for(iVarcoef in iName.varcoef){ ## iVarcoef <- iName.varcoef[1]
                        REML.num[,,iVarcoef] <- REML.num[,,iVarcoef] + iX_OmegaM1_dOmega_OmegaM1_X[iDouble2Mat,iVarcoef]
                    }
                }
            }
        }
    }

    ## ** export
    if(REML && test.vcov){
        REML.denomM1 <- solve(REML.denom)

        ## compute: 0.5 tr((X\OmegaM1X)^-1 (X\OmegaM1 d\Omega \OmegaM1 X)) in one go
        ## Score[name.allvarcoef] <-  Score[name.allvarcoef] + 0.5 * apply(REML.num, MARGIN = 3, function(x){tr(REML.denomM1 %*% x)})
        Score[name.allvarcoef] <-  Score[name.allvarcoef] + 0.5 * as.double(REML.denomM1) %*% matrix(REML.num, nrow = prod(dim(REML.num)[1:2]), ncol = dim(REML.num)[3], byrow = FALSE)
    }
    attr(Score,"message") <- message
    return(Score)
}


##----------------------------------------------------------------------
### score.R ends here
