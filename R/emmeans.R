### emmeans.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 10 2021 (16:08) 
## Version: 
## Last-Updated: mar 12 2024 (18:25) 
##           By: Brice Ozenne
##     Update #: 93
##----------------------------------------------------------------------
## 
### Commentary: 
##  No exported anymore to avoid dependency on emmeans
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * recover_data.lmm (code)
##' @title Link to emmeans package
##' @description Link to emmeans package. Not meant for direct use.
##' @noRd
##' 
##' @param object a \code{lmm} object.
##' @param trms see \code{emmeans::emm_basis} documentation 
##' @param xlev see \code{emmeans::emm_basis} documentation 
##' @param grid see \code{emmeans::emm_basis} documentation 
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return dataset or list used by the emmeans package.
recover_data.lmm <- function(object, ...){
    requireNamespace("emmeans")
    fcall <- object$call
    ff <- stats::formula(object, effects = "mean")
    data <- object$data
    
    oterms <- stats::terms(stats::model.frame(ff, data = data))
    
    out <- emmeans::recover_data(fcall, trms = stats::delete.response(oterms), na.action = NULL, frame = data)
    return(out)
}

## * emm_basis.lmm (code)
emm_basis.lmm <- function(object, trms, xlev, grid, ...){

    out <- list()
    m  <-  stats::model.frame(trms, grid, na.action = stats::na.pass, xlev = xlev)
    out$X <- stats::model.matrix(object, data = m, effects = "mean")
    ## out$X  <-  stats::model.matrix(trms, m, contrasts.arg = object$contrasts)
    out$bhat  <- stats::coef(object, effects = "mean")
    out$nbasis  <-  matrix(NA)  ## no rank deficiency
    out$V  <- stats::vcov(object, effects = "mean")

    if(!is.null(object$dVcov)){
        out$dffun <- function(k,dfargs){
            if(is.vector(k)){k <- rbind(k)}
            if(is.null(colnames(k))){colnames(k) <- dfargs$name.meanparam}
            do.call(dfargs$FUN, args = list(X.beta = k, vcov.param = dfargs$vcov.param, dVcov.param = dfargs$dVcov.param))
        }    
    }else{
        out$dffun <- function(k,dfargs){
            return(Inf)
        }
    }
    out$dfargs <- list(FUN = .dfX,
                       name.meanparam = colnames(out$X),
                       vcov.param = object$vcov,
                       dVcov.param = object$dVcov)
    
    return(out)
}


##----------------------------------------------------------------------
### emmeans.R ends here
