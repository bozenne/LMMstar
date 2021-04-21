### residuals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:40) 
## Version: 
## Last-Updated: Apr 20 2021 (15:54) 
##           By: Brice Ozenne
##     Update #: 47
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * residuals.lmm (code)
##' @export
residuals.lmm <- function(object, data = NULL, p = NULL, type.object = "lmm", type.residual = "response", format = "wide"){
    
    ## ** normalize user imput
    type.object <- match.arg(type.object, c("lmm","gls"))
    format <- match.arg(format, c("wide","long"))
    type.residuals <- match.arg(type.residual, c("response","pearson","normalized"))

    ## ** extract
    if(type.object == "lmm"){

        if(!is.null(data)){
            ff.allvars <- c(all.vars(object$formula$mean), all.vars(object$formula$var))
            if(any(ff.allvars %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing variable(s) \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
            }

            design <- .model.matrix.lmm(formula.mean = object$formula$mean.design,
                                        formula.var = object$formula$var.design,
                                        data = data,
                                        var.outcome = object$outcome$var,
                                        var.strata = object$strata$var, U.strata = object$strata$levels,
                                        var.time = object$time$var, U.time = object$time$levels,
                                        var.cluster = object$cluster$var,
                                        structure = object$structure
                                        )
            Y <- design$Y
            X <- design$X.mean
            X.var <- design$X.var
            n.cluster <- design$cluster$n
            index.cluster <- design$index.cluster
            index.variance <- design$index.vargroup
            index.time <- design$index.time
        }else{
            Y <- object$design$Y
            X <- object$design$X.mean
            X.var <- object$design$X.var
            n.cluster <- object$design$cluster$n
            index.cluster <- object$design$index.cluster
            index.variance <- object$design$index.vargroup
            index.time <- object$design$index.time
        }

        if(!is.null(p)){
            if(any(names(object$param$mu) %in% names(p) == FALSE)){
                stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(object$param$mu)[names(object$param$mu) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
            }
            beta <- p[names(object$param$mu)]
            if(type.residuals %in% c("pearson","normalized")){
                Omega <- attr(X.var,"FUN.Omega")(object = X.var, sigma = p[names(object$param$sigma)], k = p[names(object$param$k)], rho = p[names(object$param$cor)])
                precision <- lapply(Omega, solve)
            }
        }else{
            beta <- object$param$mu
            if(type.residuals %in% c("pearson","normalized")){
                Omega <- object$Omega
                precision <- object$OmegaM1
            }
        }
        if(!is.null(data) || !is.null(p)){
            res <-  Y - X.mean %*% mu
        }else{
            res <- object$residuals
        }

        ## normalization
        if(type.residuals %in% c("pearson","normalized")){
            for(iId in 1:n.cluster){ ## iId <- 7
                iIndex <- which(index.cluster==iId)
                iOrder <- order(object$design$index.time[iIndex])
                iResidual <- res[iIndex,,drop=FALSE][iOrder]
                if(type.residuals == "pearson"){
                    resnorm <- res[index.cluster==iId]/sqrt(diag(Omega[[index.variance[iId]]]))
                }else if(type.residuals == "normalized"){
                    resnorm <- as.double(res[index.cluster==iId] %*% precision[[index.variance[iId]]])
                }
                res[index.cluster==iId] <- resnorm[sort(iOrder)]
            }
        }

        if(format=="long"){
            res <- reshape2::dcast(data = data.frame(residuals = res, cluster = index.cluster, time = index.time),
                                   formula = cluster~time, value.var = "residuals")
            names(res) <- c(object$cluster$var, object$time$levels)
        }
        return(res)

    }else if(type.object == "gls"){
        if(object$strata$n == 1){
            return(residuals(object$gls, type = type.residual))
        }else {
            return(lapply(object$gls, residuals, type = type.residual))
        }
    }

}

##----------------------------------------------------------------------
### residuals.R ends here
