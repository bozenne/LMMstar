### logLik.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (17:26) 
## Version: 
## Last-Updated: okt  3 2024 (11:05) 
##           By: Brice Ozenne
##     Update #: 439
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * logLik.lmm (documentation)
##' @title Log-Likelihood From a Linear Mixed Model
##' @description Extract or compute the log-likelihood of a linear mixed model.
##' 
##' @param object a \code{lmm} object.
##' @param newdata [data.frame] dataset relative to which the log-likelihood should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param indiv [logical] Should the contribution of each cluster to the log-likelihood be output? Otherwise output the sum of all clusters of the derivatives. 
##' @param p [numeric vector] value of the model coefficients at which to evaluate the log-likelihood. Only relevant if differs from the fitted values.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details \bold{indiv}: only relevant when using maximum likelihood. Must be \code{FALSE} when using restricted maximum likelihood.
##' 
##' @return A numeric value (total logLikelihood) or a vector of numeric values, one for each cluster (cluster specific logLikelihood).
##' 
##' @keywords methods

## * logLik.lmm (code)
##' @export
logLik.lmm <- function(object, newdata = NULL, p = NULL, indiv = FALSE, ...){


    ## ** normalize user input

    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    dots$options <- NULL
    
    ## ** extract or recompute log-likelihood
    if(is.null(newdata) && is.null(p) && indiv == FALSE){

        design <- object$design ## useful in case of NA
        out <- object$logLik

    }else{

        ## normalize newdata argument
        if(!is.null(newdata)){
            design <- stats::model.matrix(object, newdata = newdata, effects = "all", simplify = FALSE)
        }else{
            design <- object$design
        }

        ## normalize p argument
        if(!is.null(p)){
            init <- .init_transform(p = p, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, 
                                    x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho,
                                    table.param = object$design$param)
            theta <- init$p
        }else{
            theta <- object$param
        }

        ## evaluate log-lik
        out <- .moments.lmm(value = theta, design = design, time = object$time, method.fit = object$args$method.fit, type.information = object$args$type.information,
                            transform.sigma = "none", transform.k = "none", transform.rho = "none",
                            logLik = TRUE, score = FALSE, information = FALSE, vcov = FALSE, df = FALSE, indiv = indiv, 
                            trace = FALSE, precompute.moments = !is.null(object$design$precompute.XX))$logLik

    } 

    ## ** name and restaure NAs
    if(indiv){

        if(!is.numeric(object$cluster$levels)){
            names(out) <- object$cluster$levels[match(1:length(out),object$cluster$index)]
        }
        out <- restaureNA(out, index.na = object$index.na,
                          level = "cluster", cluster = object$cluster)        

    }


    ## ** export
    return(out)
}

## * logLik.mlmm (documentation)
##' @title Log-Likelihood From Multiple Linear Mixed Models
##' @description Extract or compute the log-likelihood from each group-specific Linear Mixed Model.
##' 
##' @param object a \code{mlmm} object.
##' @param newdata [data.frame] dataset relative to which the log-likelihood should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param indiv [logical] Should the contribution of each cluster to the log-likelihood be output? Otherwise output the sum of all clusters of the derivatives. 
##' @param p [list of numeric vector] values for the model parameters at which to evaluate the log-likelihood.
##' Only relevant if differs from the fitted values.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details \bold{indiv}: only relevant when using maximum likelihood. Must be \code{FALSE} when using restricted maximum likelihood.
##' 
##' @return A numeric vector when \code{indiv=FALSE} (log-likelihood per model)
##' or a matrix of numeric values when when \code{indiv=TRUE} with one line for each cluster and one column for each model.
##' 
##' @keywords methods

## * logLik.mlmm (code)
##' @export
logLik.mlmm <- function(object, newdata = NULL, p = NULL, indiv = FALSE, ...){

    ## ** normalize user input

    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    dots$options <- NULL
    
    ## *** p
    if(!is.null(p)){
        
        if(!is.list(p)){
            stop("Argument \'p\' should either be NULL or a list. \n")
        }
        if(is.null(names(p))){
            stop("Argument \'p\' should either be NULL or a named list. \n")
        }
        if(any(names(p) %in% names(object$model) == FALSE)){
            stop("Incorrect names for argument \'p\': \"",paste(setdiff(names(p),names(object$model)), collapse = "\", \""),"\". \n", 
                 "Should be among \"",paste(names(object$model), collapse = "\", \""),"\". \n")
        }

    }

    ## ** extract log-likelihood
    if(indiv){
        out <- do.call(cbind,lapply(names(object$model), function(iBy){
            stats::logLik(object$model[[iBy]], indiv = TRUE, p = p[[iBy]], newdata = newdata)
        }))
    }else{
        out <- sapply(names(object$model), function(iBy){
            stats::logLik(object$model[[iBy]], indiv = FALSE, p = p[[iBy]], newdata = newdata)
        })
    }

    ## ** export
    return(out)
}

## * .logLik
.logLik <- function(X, residuals, precision, weights, 
                    pattern, index.cluster, indiv, REML, precompute){

    ## ** extract information
    if(indiv && REML){##  https://towardsdatascience.com/maximum-likelihood-ml-vs-reml-78cf79bef2cf
        stop("Cannot compute individual likelihood contribution with REML. \n")
    }
    n.cluster <- length(pattern) 
    log2pi <- log(2*pi)
    n.mucoef <- NCOL(X)

    ## ** prepare output
    compute.indiv <- indiv || is.null(precompute$weights) || is.null(precompute$RR)
    if(compute.indiv){
        ll <- rep(NA, n.cluster)
    }else if(any(sapply(precision, inherits, "try-error"))){ ## when evaluating log-likelihood at parameter values where the residual variance-covariance matrix is singular
        return(NA)
    }else{
        ll <- 0
    }

    ## ** REML contribution
    if(REML){
        if(is.null(precompute$REML)){
            X.OmegaM1.X <- matrix(0, nrow = n.mucoef, ncol = n.mucoef)
            for(iId in 1:n.cluster){ ## iId <- 1
                iX <- X[index.cluster[[iId]],,drop=FALSE]
                iOmegaM1 <- precision[[pattern[iId]]]
                iWeights <- weights[iId]
                X.OmegaM1.X <- X.OmegaM1.X + iWeights * t(iX) %*% iOmegaM1 %*% iX               
            }
            precompute$REML <- log(det(X.OmegaM1.X))
        }
        ll <- ll - precompute$REML$logdet_X.OmegaM1.X/2 + n.mucoef * log2pi/2            
    }

    ## ** compute log-likelihood
    if(compute.indiv){

        ## *** looping over individuals
        for (iId in 1:n.cluster){ ## iId <- 1
            iResidual <- residuals[index.cluster[[iId]], , drop = FALSE]
            iOmegaM1 <- precision[[pattern[iId]]]
            ll[iId] <- - weights[iId] * (NCOL(iOmegaM1) * log2pi - attr(iOmegaM1,"logdet") + t(iResidual) %*% iOmegaM1 %*% iResidual)/2
        }
        if(!indiv){
            ll <- sum(ll)
        }

    }else{

        ## *** looping over covariance patterns
        for (iPattern in 1:length(precision)) { ## iPattern <- 1
            ## precompute$RR has already been weigthed
            iOmegaM1 <- precision[[iPattern]]
            ll <- ll - 0.5 * unname(precompute$weights[iPattern]) * (NCOL(iOmegaM1) * log2pi - attr(iOmegaM1,"logdet")) - 0.5 * sum(precompute$RR[[iPattern]] * attr(iOmegaM1,"vectorize"))
        }

    }

    ## ** export
    return(ll)
}


##----------------------------------------------------------------------
### logLik.R ends here
