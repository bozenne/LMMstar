### iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  4 2021 (10:04) 
## Version: 
## Last-Updated: jul 26 2024 (17:57) 
##           By: Brice Ozenne
##     Update #: 211
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
##' @description Extract the influence function of linear mixed model parameters.
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
##' @return A matrix with one row per observation and one column per parameter.
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
    if(is.null(effects)){
        if((is.null(transform.sigma) || identical(transform.sigma,"none")) && (is.null(transform.k) || identical(transform.k,"none")) && (is.null(transform.rho) || identical(transform.rho,"none"))){
            effects <- options$effects
        }else{
            effects <- c("mean","variance","correlation")
        }
    }else if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }else{
        effects <- match.arg(effects, c("mean","fixed","variance","correlation","ranef"), several.ok = TRUE)
        effects[effects== "fixed"] <- "mean"
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

    ## ** get inverse information, score, and compute iid
    if(x$args$type.information == "expected" && x$args$method.fit == "ML"){
        effects2 <- effects
    }else{
        ## get the score and the vcov for all 
        effects2 <- c("mean","variance","correlation")
    }

    outMoments <- .moments.lmm(value = theta, design = x$design, time = x$time, method.fit = x$args$method.fit, type.information = type.information,
                               transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                               logLik = FALSE, score = TRUE, information = FALSE, vcov = TRUE, df = FALSE, indiv = TRUE, effects = effects2, robust = FALSE,
                               trace = FALSE, precompute.moments = !is.null(x$design$precompute.XX), transform.names = transform.names)
    x.score <- outMoments$score
    x.vcov <- outMoments$vcov
    if(x$args$type.information == "expected" && x$args$method.fit == "ML"){
        keep.names <- colnames(outMoments$score)
    }else{
        keep.names <- names(coef.lmm(x, effects = effects,
                                     transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names))
    }
    out <- x.score %*% x.vcov[,keep.names,drop=FALSE]
    if(robust==FALSE){
        out <- sweep(out, MARGIN = 2, FUN = "*", STATS = sqrt(diag(x.vcov[keep.names,keep.names,drop=FALSE]))/sqrt(colSums(out^2, na.rm = TRUE)))
    }
    attr(out, "message") <- attr(x.score,"message")


    ## ** export
    return(out)
}

## * iid.Wald_lmm
##' @title Extract the Influence Function from Wald Tests
##' @description Extract the influence function of linear mixed model parameters involved in the Wald test.
##' @rdname iid.Wald_lmm
##' @rdname influence.Wald_lmm
##'
##' @param object a \code{Wald_lmm} object.
##' @param method [character] should the influence function of the linear contrasts involved in the Wald test (\code{"none"})
##' or of the linear mixed model parameters (\code{"all"}) be output?
##' Can also be \code{"test"} to test whether the influence function has been stored.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @return A matrix with one row per observation and one column per parameter (\code{method="all"} or \code{method="contrast"}) or a logical value (\code{method="test"}).
##' 
##' @export
iid.Wald_lmm <- function(x, method = "none", ...){

    options <- LMMstar.options()
    adj.method <- options$adj.method

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** object
    if(x$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling rbind.Wald_lmm. \n")
        return(invisible(NULL))
    }

    ## *** method
    if(!is.character(method) || !is.vector(method)){
        stop("Argument \'method\' must be a character.")
    }
    if(length(method)!=1){
        stop("Argument \'method\' must have length 1.")
    }    
    valid.method <- c("none","all","test",adj.method)
    if(any(method %in% valid.method == FALSE)){
        stop("Incorrect value for argument \'method\': \"",paste(setdiff(method,valid.method), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.method, collapse ="\", \""),"\". \n")
    }    
    if(method %in% adj.method){
        method <- "none"
    }
 
    ## ** extract
    if(method == "test"){
        return(!is.null(x$glht[[1]]$iid))
    }else if(x$args$type=="auto"){
        message("The influence function has not been stored. \n",
                "Consider specifying the argument \'method\' when calling anova with explicit contrast (e.g. via a matrix or equations). \n")
        return(NULL)

    }

    out <- x$glht[[1]]$iid
    if(method=="none"){
        contrast <- model.tables(x, method = "contrast")
        out <- out[,colnames(contrast),drop=FALSE] %*% t(contrast)
        attr(out,"message") <- attr(x$glht[[1]]$iid,"message")
    }else{ ## restaure transformed names
        colnames(out) <- names(coef(x, method = "all"))
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
