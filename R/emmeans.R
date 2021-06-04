### emmeans.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 10 2021 (16:08) 
## Version: 
## Last-Updated: Jun  4 2021 (09:08) 
##           By: Brice Ozenne
##     Update #: 45
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * recover_data.lmm (code)
##' @title Link to emmeans package
##' @description Link to emmeans package. For internal use.
##' 
##' @name LMMstar2emmeans
##' 
##' @param object a \code{lmm} object.
##' @param trms see \code{emmeans::emm_basis} documentation 
##' @param xlev see \code{emmeans::emm_basis} documentation 
##' @param grid see \code{emmeans::emm_basis} documentation 
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @method recover_data lmm
##' @export
recover_data.lmm <- function(object, ...){
    fcall <- object$call
    ff <- stats::formula(object, effects = "mean")
    data <- object$data
    
    oterms <- stats::terms(stats::model.frame(ff, data = data))
    out <- recover_data(fcall, trms = stats::delete.response(oterms), na.action = NULL, frame = data)
    return(out)
}

## * emm_basis.lmm (code)
##' @rdname LMMstar2emmeans
##' @method emm_basis lmm
##' @export
emm_basis.lmm <- function(object, trms, xlev, grid, ...){

    out <- list()
    m  <-  stats::model.frame(trms, grid, na.action = stats::na.pass, xlev = xlev)    
    out$X  <-  stats::model.matrix(trms, m, contrasts.arg = object$contrasts) 
    out$bhat  <- stats::coef(object, effects = "mean")
    out$nbasis  <-  matrix(NA)  ## no rank deficiency
    out$V  <- stats::vcov(object, effects = "mean")
    out$dffun <- function(k,dfargs){
        do.call(dfargs$FUN, args = list(X.beta = k, vcov.param = dfargs$vcov.param, dVcov.param = dfargs$dVcov.param))
    }
    out$dfargs <- list(FUN = .dfX , vcov.param = vcov(object, effects = "all"), dVcov.param = object$dVcov)

    return(out)
}


##----------------------------------------------------------------------
### emmeans.R ends here
