### df.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 18 2021 (10:34) 
## Version: 
## Last-Updated: sep  7 2021 (17:11) 
##           By: Brice Ozenne
##     Update #: 35
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .df
.df <- function(value, reparametrize,
                design, time, method.fit, type.information,
                transform.sigma, transform.k, transform.rho, effects, robust, diag,
                precompute.moments = precompute.moments, method.numDeriv = method.numDeriv){

    ## ** prepare vector of parameters
    param.value <- value
    param.type <- design$param$type
    param.strata <- design$param$strata
    name.allcoef <- names(param.type)
    n.allcoef <- length(param.type)

    param.nameVar <- name.allcoef[param.type %in% c("sigma","k","rho")]
    param.nameMean <- name.allcoef[param.type %in% c("mu")]

    test.transform <- (transform.sigma != "none") || (transform.k != "none") || (transform.rho != "none")

    param.trans.value <- c(param.value[param.nameMean],reparametrize$p)[name.allcoef]

    ## ** warper for computing information
    FUN_information <- function(p, effects, as.double){

        if(test.transform){ ## back-transform
            backp <- .reparametrize(p = p[param.nameVar], type = param.type[param.nameVar], strata = param.strata[param.nameVar], 
                                    Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                    transform.sigma = transform.sigma,
                                    transform.k = transform.k,
                                    transform.rho = transform.rho,
                                    transform.names = FALSE)
            p[param.nameVar] <- backp$p
        }

        iMoment <- .moments.lmm(value = p, design = design, time = time, method.fit = method.fit, type.information = type.information,
                                transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                logLik = FALSE, score = FALSE, information = TRUE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = effects, robust = robust,
                                trace = FALSE, precompute.moments = precompute.moments, method.numDeriv = method.numDeriv, transform.names = FALSE)

        if(as.double){
            return(as.double(iMoment$information))
        }else{
            return(iMoment$information)
        }
    }

    ## ** variance-covariance matrix
    ## full variance covariance matrix
    info <- FUN_information(param.trans.value, effects = c("mean","variance","correlation"), as.double = FALSE)
    vcov <- solve(info)

    ## only for the coefficient of interest
    info.effects <- FUN_information(param.trans.value, effects = effects, as.double = FALSE)
    name.effects <- colnames(info.effects)
    n.effects <- length(name.effects)
    vcov.effects <- vcov[name.effects,name.effects,drop=FALSE]
    ## matrix(FUN_information(param.trans.value, as.double = TRUE), nrow = n.allcoef, ncol = n.allcoef)

    ## ** derivative of the information using numerical derivative
    if(type.information == "observed"){
        M.dInfo <- numDeriv::jacobian(func = function(p){FUN_information(p, effects = effects,as.double = TRUE)}, x = param.trans.value, method = method.numDeriv)
        colnames(M.dInfo) <- name.allcoef
    }else{
        M.dInfo <- numDeriv::jacobian(func = function(p){FUN_information(c(param.value[param.nameMean],p)[name.allcoef], effects = effects, as.double = TRUE)}, x = param.trans.value[param.nameVar], method = method.numDeriv)
        colnames(M.dInfo) <- param.nameVar
    }

    A.dVcov <- array(0, dim = c(n.effects,n.effects,n.allcoef), dimnames = list(name.effects,name.effects,name.allcoef))
    for(iParam in 1:NCOL(M.dInfo)){ ## iParam <- 1
        iName <- colnames(M.dInfo)[iParam]
        A.dVcov[,,iName] <- - vcov.effects %*% matrix(M.dInfo[,iName], nrow = n.effects, ncol = n.effects) %*% vcov.effects
    }

    ## solve(crossprod(model.matrix(e.lmm, effects = "mean")))
    ## 4*coef(e.lmm)["sigma"]^2/stats::nobs(e.lmm)[1]
    ## ** degrees of freedom
    if(diag){
        df <- stats::setNames(sapply(1:n.effects, function(iP){
            2 * vcov.effects[iP,iP]^2 / (A.dVcov[iP,iP,] %*% vcov %*% A.dVcov[iP,iP,])
        }), name.effects)
    }else{
        df <- matrix(NA, nrow = n.effects, ncol = n.effects, dimnames = list(name.effects, name.effects))
        for(iParam in 1:n.effects){
            for(iiParam in 1:iParam){
                df[iParam,iiParam] <- 2 * vcov.effects[iParam,iiParam]^2 / (A.dVcov[iParam,iiParam,] %*% vcov %*% A.dVcov[iiParam,iParam,])
                if(iParam != iiParam){
                    df[iiParam,iParam] <- df[iParam,iiParam]
                }
            }
        }
    }
    
    ## ** export
    attr(df,"dVcov") <- A.dVcov
    return(df)
}

##----------------------------------------------------------------------
### df.R ends here
