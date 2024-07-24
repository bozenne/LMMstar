### iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  4 2021 (10:04) 
## Version: 
## Last-Updated: jul 24 2024 (16:00) 
##           By: Brice Ozenne
##     Update #: 162
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
##' @description Extract the influence function from a linear mixed model.
##' @rdname iid.lmm
##' @rdname influence.lmm
##' 
##' @param x a \code{lmm} object.
##' @param effects [character] Should the influence function for all coefficients be output (\code{"all"}),
##' or only for coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only for coefficients relative to the variance structure (\code{"variance"}),
##' or only for coefficients relative to the correlation structure (\code{"correlation"}).
##' @param p [numeric vector] value of the model coefficients at which to evaluate the influence function. Only relevant if differs from the fitted values.
##' @param robust [logical] If \code{FALSE} the influence function is rescaled to match the model-based standard errors.
##' The correlation however will not (necessarily) match the model-based correlation.
##' @param type.information [character] Should the expected information be used  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param REML2ML [logical] Should the individual score contributions be calculated based on the Maximum Likelihood (ML) equations even when the model was fitted using Restricted Maximum Likelihood (REML)?
##' 
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details The influence function equals the individual score rescaled by the (inverse) information.
##' With the expected information and for a lmm fitted using ML, the information is block diagonal so the influence function for the mean and variance parameters can be computed separately.
##' Otherwise the information and individual score relative to all model parameters should be considered. The later is probablematic when using REML and two add-hoc solutions are considered:
##' \itemize{
##' \item \code{REML2ML=TRUE}: evaluate the individual score for the variance parameters using ML score equations instead of REML,
##' e.g., neglecting the product of the REML score contribution for the covariance parameters with the cross terms between mean and variance parameters in the inverse information. 
##' \item \code{REML2ML=FALSE}: only use the information and individual score relative to the mean parameters,
##' e.g., neglecting the cross-terms between mean and variance parameters in the information. 
##' }
##'
##' @keywords methods
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

    options <- LMMstar.options()

    ## ** normalize user imput

    ## *** dots
    dots <- list(...)
    dots$complete <- NULL ## for multcomp which passes an argument complete when calling vcov
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** effects
    if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","variance","correlation"), several.ok = TRUE)

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

    ## *** REML2ML
    indiv <- TRUE
    if(is.null(REML2ML)){
        REML2ML <- options$REML2ML
    }
    if(REML2ML){
        attr(indiv,"REML2ML") <- TRUE
    }

    ## ** get inverse information, score, and compute iid
    if((x$args$type.information == "expected" && x$args$method.fit == "ML") || (!REML2ML && x$args$method.fit == "REML")){
        x.score <- score.lmm(x, effects = effects, p = p, indiv = indiv, 
                             transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        x.vcov <- vcov.lmm(x, effects = effects, p = p, robust = FALSE, type.information = type.information, df = FALSE,
                           transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        out <- x.score %*% x.vcov
        if(x$args$method.fit=="ML"){
            message <- NULL
        }else{
            message <- "neglecting information cross-terms"
        }
    }else{
        keep.names <- names(coef.lmm(x, effects = effects,
                                     transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names))
        x.scoreALL <- score.lmm(x, effects = "all", p = p, indiv = indiv, 
                                transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        x.vcovALL <- vcov.lmm(x, effects = "all", p = p, robust = FALSE, type.information = type.information, df = FALSE,
                              transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        x.vcov <- x.vcovALL[keep.names,keep.names,drop=FALSE]

        out <- x.scoreALL %*% x.vcovALL[,keep.names,drop=FALSE]
        message <- "neglecting product REML score with vcov cross-terms"
    }

    if(robust==FALSE){
        out <- sweep(out, MARGIN = 2, FUN = "*", STATS = sqrt(diag(x.vcov))/sqrt(colSums(out^2, na.rm = TRUE)))
    }

    ## ** export
    attr(out,"message") <- message
    return(out)
}

## * iid.Wald_lmm (code)
##' @export
iid.Wald_lmm <- function(x, effects = "contrast", ...){

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character.")
    }
    if(length(effects)!=1){
        stop("Argument \'effects\' must have length 1.")
    }    
    valid.effects <- c("contrast","all")

    ## ** extract
    if(x$args$type=="auto"){
        message("The influence function has not been stored. \n",
                "Consider specifying the argument \'effects\' when calling anova with explicit contrast (e.g. via a matrix or equations). \n")
        return(NULL)

    }
    out <- x$glht[[1]]$iid
    if(effects=="contrast"){
        contrast <- model.tables(x, effects = "contrast")
        out <- out[,colnames(contrast),drop=FALSE] %*% t(contrast)
    }

    ## ** export
    return(out)
}

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
        contrast <- coef(x, type = "contrast") ## contrast combinations across models
        ls.iid <- lapply(names(x$model), function(iBy){
            iid(x$model[[iBy]], effects = effects2, p = p[[iBy]], ...)
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
            iIID <- iid(iO, p = p[[iBy]], effects = effects, ...)

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
        name.order <- names(coef(x, effects = "contrast", ordering = ordering))        
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


##----------------------------------------------------------------------
### iid.R ends here
