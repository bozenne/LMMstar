### iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  4 2021 (10:04) 
## Version: 
## Last-Updated: aug  8 2024 (15:17) 
##           By: Brice Ozenne
##     Update #: 339
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * iid.lmm (documentation)
##' @title Extract the Influence Function From a Linear Mixed Model
##' @description Extract the influence function for an ML or REML estimator of parameters from a linear mixed model.
##' 
##' @param x a \code{lmm} object.
##' @param effects [character] Should the influence function for all coefficients be output (\code{"all"}),
##' or only for coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only for coefficients relative to the variance structure (\code{"variance"}),
##' or only for coefficients relative to the correlation structure (\code{"correlation"}).
##' Can also contain \code{"gradient"} to also output the gradient of the influence function.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the influence function. Only relevant if differs from the fitted values.
##' @param robust [logical] If \code{FALSE} the influence function is rescaled to match the model-based standard errors.
##' The correlation however will not (necessarily) match the model-based correlation.
##' @param type.information [character] Should the expected information be used  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' 
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details The influence function equals the individual score rescaled by the (inverse) information.
##' With the expected information and for a lmm fitted using ML, the information is block diagonal so the influence function for the mean and variance parameters can be computed separately.
##' Otherwise the information and individual score relative to all model parameters should be considered.
##' The later is probablematic when using REML as the REML term is the ratio of two term linear in the individual contributions which is not itself linear in the individual contributions.
##' As an add-hoc solution, the denominator is treated as fixed so the ratio is decomposed w.r.t. its numerator.
##'
##' @keywords methods
##' @return A matrix with one row per observation and one column per parameter. \itemize{
##' \item \code{df=TRUE}: with an attribute \code{"df"} containing a numeric vector with one element per parameter. \cr
##' \item \code{effects} includes \code{"gradient"}: with an attribute \code{"gradient"} containing a 3 dimensional array with dimension the number of parameters.
##' }
##' 
##' @examples
##' data(gastricbypassL)
##' 
##' #### Case WITHOUT cross terms ####
##' e.lmmREML <- lmm(weight ~ visit, method.fit = "REML", df = FALSE,
##'                  repetition = ~ visit|id, data = gastricbypassL)
##' e.lmmML <- lmm(weight ~ visit, method.fit = "ML", df = FALSE,
##'               repetition = ~ visit|id, data = gastricbypassL)
##'
##' name.mu <- names(coef(e.lmmREML, effects = "mean"))
##' name.sigmakrho <- names(coef(e.lmmREML, effects = c("variance","correlation")))
##' info.REML <- information(e.lmmREML, effects = "all", transform.names = FALSE)
##' info.ML <- information(e.lmmML, effects = "all", transform.names = FALSE)
##' info.REML2ML <- information(e.lmmML, p = coef(e.lmmREML, effects = "all"),
##'                             effects = "all", transform.names = FALSE)
##' 
##' range(info.REML[name.mu,name.sigmakrho])
##' range(info.ML[name.mu,name.sigmakrho])
##' range(info.REML[name.mu,]-info.REML2ML[name.mu,])
##' range(iid(e.lmmREML, REML2ML = TRUE) - iid(e.lmmREML, REML2ML = FALSE))
##' ## neglectable differences
##' 
##' #### Case WITH cross terms ####
##' e2.lmmREML <- lmm(glucagonAUC ~ visit + weight, method.fit = "REML", df = FALSE,
##'                  repetition = ~ visit|id, data = gastricbypassL)
##' e2.lmmML <- lmm(glucagonAUC ~ visit + weight, method.fit = "ML", df = FALSE,
##'               repetition = ~ visit|id, data = gastricbypassL)
##'
##' name2.mu <- names(coef(e2.lmmREML, effects = "mean"))
##' name2.sigmakrho <- names(coef(e2.lmmREML, effects = c("variance","correlation")))
##' info2.REML <- information(e2.lmmREML, effects = "all", transform.names = FALSE)
##' info2.ML <- information(e2.lmmML, effects = "all", transform.names = FALSE)
##' info2.REML2ML <- information(e2.lmmML, p = coef(e2.lmmREML, effects = "all"),
##'                             effects = "all", transform.names = FALSE)
##' 
##' range(info2.REML[name.mu,]-info2.REML2ML[name.mu,])
##' ## neglectable difference
##' range(info2.REML[name.mu,name.sigmakrho])
##' range(info2.ML[name.mu,name.sigmakrho])
##' range(iid(e2.lmmREML, REML2ML = TRUE) - iid(e2.lmmREML, REML2ML = FALSE))
##' ## non-neglectable differences
##' diag(crossprod(iid(e2.lmmREML, REML2ML = TRUE)))/diag(vcov(e2.lmmREML))
##' diag(crossprod(iid(e2.lmmREML, REML2ML = FALSE)))/diag(vcov(e2.lmmREML))

## * iid.lmm (code)
##' @export
iid.lmm <- function(x,
                    effects = "mean",
                    p = NULL,
                    robust = TRUE,
                    type.information = NULL,
                    transform.sigma = NULL,
                    transform.k = NULL,
                    transform.rho = NULL,
                    transform.names = TRUE,
                    REML2ML = NULL,
                    ...){

    ## ** normalize user imput

    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    dots$complete <- NULL ## for multcomp which passes an argument complete when calling vcov
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
        valid.effects <- c("mean","fixed","variance","correlation","all","gradient","dVcov")
        if(any(effects %in% valid.effects == FALSE)){
            stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
                 "Valid values: \"",paste(setdiff(valid.effects,"dVcov"), collapse ="\", \""),"\". \n")
        }
        if("dVcov" %in% effects){
            iid.dVcov <- TRUE
            effects <- setdiff(effects, "dVcov")
        }else{
            iid.dVcov <- FALSE
        }
        if("gradient" %in% effects){
            keep.grad <- TRUE
            effects <- setdiff(effects, "gradient")
            if(x$args$df==0){
                stop("Argument \'effects\' cannot contain \"gradient\" when no degrees-of-freedom have been stored. \n",
                     "Consider setting the argument \'df\' to TRUE when calling lmm. \n")
            }
            if(length(effects)==0){
                stop("Argument \'effects\' cannot only contain \"gradient\". \n",
                     "Consider setting adding \"mean\" or \"all\". \n")
            }
        }else{
            keep.grad <- FALSE
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

    ## *** type.information
    if(is.null(type.information)){
        type.information <- x$args$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    ## *** p and transform
    init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = x$reparametrize$transform.sigma, x.transform.k = x$reparametrize$transform.k, x.transform.rho = x$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform
    if(is.null(p)){
        theta <- x$param
    }else{
        theta <- init$p
    }    

    if(is.null(p) && test.notransform && x$args$type.information==type.information){
        recompute.vcov <- FALSE
    }else{
        recompute.vcov <- TRUE
    }

    ## ** get moments
    if(x$args$type.information == "expected" && x$args$method.fit == "ML" && !keep.grad && !iid.dVcov){
        effects2 <- effects
    }else{
        ## get the score and the vcov for all 
        effects2 <- c("mean","variance","correlation")
    }
    outMoments <- .moments.lmm(value = theta, design = x$design, time = x$time, method.fit = x$args$method.fit, type.information = type.information,
                               transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                               logLik = FALSE, score = TRUE, information = keep.grad || iid.dVcov, vcov = recompute.vcov, df = recompute.vcov && (keep.grad || iid.dVcov),
                               indiv = TRUE, effects = effects2, robust = FALSE,
                               trace = FALSE, precompute.moments = !is.null(x$design$precompute.XX), method.numDeriv = options$method.numDeriv, transform.names = transform.names)
    x.score <- outMoments$score
    if(recompute.vcov){
        x.vcov <- outMoments$vcov
    }else{ ## must be all effects, i.e. keep all columns
        x.vcov <- x$vcov
        dimnames(x.vcov) <- list(colnames(x.score),colnames(x.score))
    }
    if(x$args$type.information == "expected" && x$args$method.fit == "ML"){
        keep.names <- colnames(outMoments$score)
    }else{
        keep.names <- stats::model.tables(x, effects = c("param",effects),
                                          transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)$trans.name
}

    ## ** compute iid
    out <- x.score %*% x.vcov[,keep.names,drop=FALSE]
    if(robust==FALSE){
        out <- sweep(out, MARGIN = 2, FUN = "*", STATS = sqrt(diag(x.vcov[keep.names,keep.names,drop=FALSE]))/sqrt(colSums(out^2, na.rm = TRUE)))
    }

    ## ** compute gradient
    if(keep.grad || iid.dVcov){
        x.hessian <- -outMoments$information
        if(recompute.vcov){
            x.dVcov <- outMoments$dVcov
        }else{  ## must be all effects, i.e. keep all columns
            x.dVcov <- x$dVcov
            dimnames(x.dVcov) <- list(colnames(x.score),colnames(x.score),colnames(x.score))
        }

        if(keep.grad){
            attr(out,"gradient") <- array(NA, dim = c(NROW(out),length(keep.names),NCOL(x.score)), dimnames = list(NULL,keep.names,colnames(x.score)))
            for(iParam in colnames(x.score)){ ## iParam <- colnames(x.score)[1]
                attr(out,"gradient")[,,iParam] <- x.hessian[,,iParam] %*% x.vcov[,keep.names,drop=FALSE] + x.score %*% x.dVcov[,keep.names,iParam]
            }
        }
        if(iid.dVcov){
            attr(out,"dVcov") <- array(NA, dim = c(NROW(out),length(keep.names),NCOL(x.score)), dimnames = list(NULL,keep.names,colnames(x.score)))
            for(iParam in colnames(x.score)){ ## iParam <- colnames(x.score)[1]
                attr(out,"dVcov")[,,iParam] <- x.hessian[,,iParam] %*% x.dVcov[,keep.names,iParam]
            }
        }
    }

    
    ## ** name and restaure NAs
    if(!is.numeric(x$cluster$levels)){
        rownames(out) <- x$cluster$levels[match(1:NROW(out),x$cluster$index)]
    } 
    out <- restaureNA(out, index.na = x$index.na,
                      level = "cluster", cluster = x$cluster)        
    if(keep.grad){
        if(!is.numeric(x$cluster$levels)){
            dimnames(attr(out,"gradient"))[[1]] <- x$cluster$levels[match(1:dim(attr(out,"gradient"))[[1]],x$cluster$index)]

        } 
        attr(out,"gradient") <- restaureNA(attr(out,"gradient"), index.na = x$index.na,
                                           level = "cluster", cluster = x$cluster)
    }
    if(iid.dVcov){
        if(!is.numeric(x$cluster$levels)){
            dimnames(attr(out,"dVcov"))[[1]] <- x$cluster$levels[match(1:dim(attr(out,"dVcov"))[[1]],x$cluster$index)]

        } 
        attr(out,"dVcov") <- restaureNA(attr(out,"dVcov"), index.na = x$index.na,
                                        level = "cluster", cluster = x$cluster)
    }

    ## ** export
    attr(out, "message") <- attr(x.score,"message")
    return(out)
}

## * iid.lmm (documentation)
##' @title Extract the Influence Function From Multiple Linear Mixed Models
##' @description Extract the influence function for ML or REML estimators of parameters from group-specific linear mixed models.

## * iid.mlmm (code)
##' @export
iid.mlmm <- function(x, effects = "contrast", p = NULL, ordering = "by", simplify = TRUE, ...){

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
        if(any(names(p) %in% names(x$model) == FALSE)){
            stop("Incorrect names for argument \'p\': \"",paste(setdiff(names(p),names(x$model)), collapse = "\", \""),"\". \n", 
                 "Should be among \"",paste(names(x$model), collapse = "\", \""),"\". \n")
        }
    }

    ## *** ordering    
    ordering <- match.arg(ordering, c("by","parameter"))
        
    ## ** extract iid
    cluster <- x$object$cluster
    n.cluster <- length(cluster)

    if(all(effects=="contrast")){
        ls.contrast <- stats::coef(x, type = "ls.contrast") ## extract linear combination from each model
        contrast <- stats::coef(x, type = "contrast") ## contrast combinations across models
        ls.iid <- lapply(names(x$model), function(iBy){
            lava::iid(x$model[[iBy]], effects = effects2, p = p[[iBy]], ...)
        })
        ls.out <- mapply(iIID = ls.iid, iC = ls.contrast, FUN = function(iIID,iC){ ## iIID <- ls.iid[[1]] ; iC <- ls.contrast[[1]]
            iOut <- matrix(0, nrow = n.cluster, ncol = NROW(iC), dimnames = list(cluster, rownames(iC)))
            iOut[rownames(iIID),] <- iIID %*% t(iC[,colnames(iIID), drop = FALSE])
            return(iOut)
        }, SIMPLIFY = FALSE)
        out <- (do.call("cbind",ls.out) %*% t(contrast))        
    }else{
        cluster <- x$object$cluster
        n.cluster <- length(cluster)

        ls.out <- lapply(names(x$model), function(iBy){
            iO <- x$model[[iBy]]
            iIID <- lava::iid(iO, p = p[[iBy]], effects = effects, ...)

            if(simplify){
                iOut <- matrix(0, nrow = n.cluster, ncol = NCOL(iIID), dimnames = list(cluster, paste0(iBy,": ", colnames(iIID))))
                attr(iOut,"name") <- colnames(iIID)
                attr(iOut,"by") <- rep(iBy, NCOL(iIID))                
            }else{
                iOut <- matrix(0, nrow = n.cluster, ncol = NCOL(iIID), dimnames = list(cluster, colnames(iIID)))
            }
            iOut[rownames(iIID),] <- iIID
            attr(iOut,"message") <- attr(iIID,"message")
            return(iOut)
        })
        names(ls.out) <- names(x$model)
        indiv <- is.matrix(ls.out[[1]])
    }

    ## ** re-order
    if(all(effects=="contrast")){
        name.order <- names(stats::coef(x, effects = "contrast", ordering = ordering))        
        out <- out[,name.order,drop=FALSE]
        attr(out,"message") <- unique(unlist(lapply(ls.iid,attr,"message")))
    }else if(simplify){
        out <- do.call("cbind",unname(ls.out))
        attr(out,"original.name") <- unname(unlist(lapply(ls.out,attr,"name")))
        attr(out,"by") <- unname(unlist(lapply(ls.out,attr,"by")))
        attr(out,"message") <- unname(unique(unlist(lapply(ls.out,attr,"message"))))

        if(ordering == "parameter"){
            reorder <- order(factor(attr(out,"original.name"),levels = unique(attr(out,"original.name"))))
            out.save <- out
            out <- out.save[,reorder,drop=FALSE]
            attr(out, "original.name") <- attr(out.save, "original.name")[reorder]
            attr(out, "by") <- attr(out.save, "by")[reorder]
            attr(out, "message") <- attr(out.save, "message")
        }

    }else if(ordering == "by"){

        out <- ls.out

    }else if(ordering == "parameter"){

        ## unique parameters
        Uname <- unique(unlist(lapply(ls.out,colnames)))

        ## combine iid
        out <- stats::setNames(lapply(Uname, function(iName){ ## iName <- "X1"
            iOut <- do.call(cbind,lapply(ls.out, function(iOut){ iOut[,iName]}))
            colnames(iOut) <- names(ls.out)
            attr(iOut,"message") <- unname(unique(unlist(lapply(ls.out, attr, "message"))))
            return(iOut)
        }), Uname)
        
    }

    ## ** export
    return(out)
}

## * iid.rbindWald_lmm (code)
##' @export
iid.rbindWald_lmm <- function(x, ...){

    message("The influence function is not stored when combining Wald tests. \n",
            "Nothing to return. \n")

    return(invisible(NULL))
}


## * iid.Wald_lmm (documentation)
##' @title Extract the Influence Function From Wald Tests Applied to a Linear Mixed Model
##' @description Extract the influence function of linear contrasts applied to an ML or REML estimator of parameters from a linear mixed model.
##'
##' @param object a \code{Wald_lmm} object.
##' @param effects [character] should the influence function for the linear contrasts involved in the Wald test be output (\code{"Wald"}),
##' or for the linear mixed model parameters (\code{"all"})?
##' @param transform.names [logical] should the name of the coefficients be updated to reflect the transformation that has been used?
##' Only relevant when \code{effects="all"}.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @return A matrix with one row per observation and one column per parameter (\code{effects="Wald"} or \code{effects="all"}) or a logical value (\code{effects="test"}).

## * iid.Wald_lmm (code)
##' @export
iid.Wald_lmm <- function(x, effects = "Wald", transform.names = TRUE, ...){

    table.param <- stats::model.tables(x, effects = "param")
    
    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- NULL
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** object
    if(x$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling rbind.Wald_lmm. \n")
        return(invisible(NULL))
    }
    
    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character. \n")
    }
    valid.effects <- c("Wald","all","test",setdiff(names(attributes(x$glht[[1]]$iid)),c("dim","dimnames","message"))) ## add dVcov if it is there    
    if(any(effects %in% valid.effects == FALSE)){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(setdiff(valid.effects,c("test","dVcov")), collapse ="\", \""),"\". \n")
    }
    if("dVcov" %in% effects){
        iid.dVcov <- TRUE
        effects <- setdiff(effects, "dVcov")
    }else{
        iid.dVcov <- FALSE
    }
    if(length(effects)!=1){
        stop("Argument \'effects\' can only contain one of \"Wald\", \"all\",  \"test\". \n")
    }

    ## ** special case: test whether there is an iid
    if(effects == "test"){
        return(!is.null(x$glht[[1]]$iid))
    }else if(x$args$simplify!=0){
        stop("Cannot extract the influence function when it has not be stored. \n",
             "Consider setting the argument \'simplify\' to 0 or FALSE when calling anova. \n")
    }

    ## ** extract iid
    out <- x$glht[[1]]$iid

    if(effects=="Wald"){

        contrast <- stats::model.tables(x, effects = "contrast", transform.names = FALSE)
        if(iid.dVcov){
            dVcov.out <- array(NA, dim = c(NROW(out), NROW(contrast), dim(attr(out,"dVcov"))[3]),
                              dimnames = list(rownames(out), rownames(contrast), dimnames(attr(out,"dVcov"))[[3]]))
            for(iParam in 1:NCOL(contrast)){ 
                dVcov.out[,,iParam] <- attr(out,"dVcov")[,colnames(contrast),iParam] %*% t(contrast)
            }
        }
        out <- out[,colnames(contrast),drop=FALSE] %*% t(contrast)
        if(iid.dVcov){
            attr(out,"dVcov") <- dVcov.out
        } ## otherwise already removed by subsetting
        attr(out,"message") <- attr(x$glht[[1]]$iid,"message")
    }

    ## ** rename
    if(transform.names && effects == "all"){ ## when effects=="Wald" the column names reflect the linear hypothesis not the model parameters
        colnames(out) <- table.param[match(colnames(out), table.param$name),"trans.name"]
        if(iid.dVcov==FALSE){
            attr(out,"dVcov") <- NULL
        }else{
            dimnames(attr(out,"dVcov")) <- list(dimnames(attr(out,"dVcov"))[[1]],
                                                table.param[match(dimnames(attr(out,"dVcov"))[[2]], table.param$name),"trans.name"],
                                                table.param[match(dimnames(attr(out,"dVcov"))[[3]], table.param$name),"trans.name"])
        }
    }

    ## ** export
    return(out)
}



##----------------------------------------------------------------------
### iid.R ends here
