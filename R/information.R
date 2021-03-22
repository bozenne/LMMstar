### information.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (22:13) 
## Version: 
## Last-Updated: mar 22 2021 (22:38) 
##           By: Brice Ozenne
##     Update #: 7
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * information.lmm (code)
##' @export
information.lmm <- function(x, data = NULL, p = NULL, ...){

    if(is.null(data) && is.null(p)){
        out <- x$information
    }else{
        if(!is.null(data)){
            ff.allvars <- c(all.vars(x$formula$mean), all.vars(x$formula$var))
            if(any(ff.allvars %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing variable(s) \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
            }

            design <- .model.matrix.lmm(formula.mean = x$formula$mean.design,
                                        formula.var = x$formula$var.design,
                                        data = data,
                                        var.outcome = x$outcome$var,
                                        var.strata = x$strata$var, U.strata = x$strata$levels,
                                        var.time = x$time$var, U.time = x$time$levels,
                                        var.cluster = x$cluster$var,
                                        structure = x$structure
                                        )
            X <- design$X.mean
            index.variance <- design$index.vargroup
            index.cluster <- design$index.cluster
            X.var <- design$X.var
        }else{
            X <- x$design$X.mean
            index.variance <- x$design$index.vargroup
            index.cluster <- x$design$index.cluster
            X.var <- x$design$X.var
        }
        if(!is.null(p)){
            if(any(names(x$param$type) %in% names(p) == FALSE)){
                stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(x$param$type)[names(x$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
            }
            Omega <- attr(X.var,"FUN.Omega")(object = X.var, sigma = p[names(x$param$sigma)], k = p[names(x$param$k)], rho = p[names(x$param$cor)])
            precision <- lapply(Omega, solve)
        }else{
            precision <- x$OmegaM1
        }
        hess <- .hessian(X = X, precision = precision,
                         index.variance = index.variance, index.cluster = index.cluster, indiv = TRUE, REML = out$method.fit=="REML")
        out <- -apply(hess, FUN = sum, MARGIN = 2:3)
        
    }
    return(out)
}


##----------------------------------------------------------------------
### information.R ends here
