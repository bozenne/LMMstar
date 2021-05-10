### emmeans.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 10 2021 (16:08) 
## Version: 
## Last-Updated: May 10 2021 (19:35) 
##           By: Brice Ozenne
##     Update #: 31
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * recover_data.lmm (code)
##' @export
recover_data.lmm <- function(object, ...){
    fcall <- object$call
    ff <- formula(object, effects = "mean")
    data <- object$data
    
    oterms <- terms(model.frame(ff, data = data))
    out <- recover_data(fcall, trms = delete.response(oterms), na.action = NULL, frame = data)
    return(out)
}

## * emm_basis.lmm (code)
## The function must obtain six things and return them in a named list. They are the:
## - matrix X of linear functions for each point in the reference grid
## - the regression coefficients bhat;
## - the variance-covariance matrix V;
## - a matrix nbasis for non-estimable functions;
## - a function dffun(k,dfargs) for computing degrees of freedom for the linear function sum(k*bhat);
## - and a list dfargs of arguments to pass to dffun.
##  - Optionally, the returned list may include a model.matrix element (the model matrix for the data or a compact version thereof obtained via .cmpMM()), which, if included, enables the submodel option.
## emmeans:::emm_basis.lm
##' @export
emm_basis.lmm <- function(object, trms, xlev, grid, ...){

    out <- list()
    m  <-  model.frame(trms, grid, na.action = na.pass, xlev = xlev)    
    out$X  <-  model.matrix(trms, m, contrasts.arg = object$contrasts) 
    out$bhat  <- coef(object, effects = "mean")
    out$nbasis  <-  matrix(NA)  ## no rank deficiency
    out$V  <- vcov(object, effects = "mean")
    out$dffun <- function(k,dfargs){
        do.call(dfargs$FUN, args = list(X.beta = k, vcov.param = dfargs$vcov.param, dVcov.param = dfargs$dVcov.param))
    }
    out$dfargs <- list(FUN = .dfX , vcov.param = vcov(object, effects = "all"), dVcov.param = attr(object$df,"dVcov"))

    return(out)
}


##----------------------------------------------------------------------
### emmeans.R ends here
