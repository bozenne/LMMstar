### getVarCov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:57) 
## Version: 
## Last-Updated: mar  5 2021 (21:48) 
##           By: Brice Ozenne
##     Update #: 22
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * getVarCov.lmm
getVarCov.lmm <- function(object, individual = NULL, simplifies = TRUE, type = c("lmm","gls")){
    type <- match.arg(type, c("lmm","gls"))

    if(type == "lmm"){
    
        if(is.null(individual)){
            out <- object$resvcov$full
        }else{
            valid.group <- names(object$lmm$index.vargroup)
            if(any(individual %in% valid.group == FALSE)){
                stop("Invalid individual. \n")
            }else{
                out <- object$resvcov$all[object$lmm$index.vargroup[individual]]
            }
        }
        if(is.list(out) && length(out)==1 && simplifies){
            return(out[[1]])
        }else{
            return(out)
        }
        
    }else{
        if(is.null(object$variable$strata)){
            if(is.null(individual)){
                return(getVarCov(object$gls[[1]]))
            }else{
                return(getVarCov(object$gls[[1]], individual = individual))
            }
        }else{
            stop("Not implemented yet!") ## not to see which model contain the individual
        }

    }
        

    return(variance)
}

## * .getVarCov
.getVarCov <- function(object, sigma, k, rho,
                       name.time, Sigma.pattern, Sigma.pattern.full){
    
    ## ** full residual variance-covariance matrix
    n.time <- length(name.time)
    if(is.null(object$modelStruct$varStruct)){
        ksigma <- sigma
    }else{
        ksigma <- sigma * k
    }

    if(is.null(object$modelStruct$corStruct)){
        variance <- lapply(ksigma, function(iKsigma){matrix( iKsigma^2, nrow = 1, ncol = 1, dimnames = list(name.time, name.time) )})
    }else if(is.null(object$modelStruct$varStruct)){
        Mcor <- matrix(1, nrow = n.time, ncol = n.time, dimnames = list(name.time, name.time))
        Mcor[lower.tri(Mcor)] <- rho
        Mcor[upper.tri(Mcor)] <- t(Mcor)[upper.tri(Mcor)]
        variance <- list(sigma^2*Mcor)
    }else{
        Mcor <- matrix(1, nrow = n.time, ncol = n.time, dimnames = list(name.time, name.time))
        Mcor[lower.tri(Mcor)] <- rho
        Mcor[upper.tri(Mcor)] <- t(Mcor)[upper.tri(Mcor)]
        variance <- setNames(lapply(Sigma.pattern, function(x){
            tcrossprod(ksigma[x]) * Mcor
        }), names(Sigma.pattern))
    }
    
    ## ** (missing data) residual variance-covariance matrix    
    if(is.null(object$modelStruct$varStruct)){
        variance.vargroup <- variance
    }else if(is.null(object$modelStruct$corStruct)){ 
        variance.vargroup <- variance
    }else{ ## variance and correlation structure
        name.vargroup <- names(Sigma.pattern)
        n.vargroup <- length(name.vargroup)
        variance.vargroup <- setNames(vector(mode = "list", length = n.vargroup), name.vargroup)
        MSigma.pattern.full <- do.call(rbind,Sigma.pattern.full)

        for(i.vargroup in 1:n.vargroup){ ## i.vargroup <- 2
            iname.vargroup <- name.vargroup[i.vargroup]
            if(iname.vargroup %in% names(Sigma.pattern.full)){
                variance.vargroup[[i.vargroup]] <- variance[[iname.vargroup]]
            }else{
                test <- t(apply(MSigma.pattern.full, 1, `%in%`, Sigma.pattern[[iname.vargroup]])) ## t() because apply(,1,) returns the transposed version
                index <- which.max(rowSums(test))
                variance.vargroup[[i.vargroup]] <- variance[[names(index)]][test[index,],test[index,],drop=FALSE]
                colnames(variance.vargroup[[i.vargroup]]) <- name.time[index]
                rownames(variance.vargroup[[i.vargroup]]) <- name.time[index]
            }
        }
    }

    ## ** export
    return(list(variance.full = variance,
                variance.vargroup = variance.vargroup))

}

##----------------------------------------------------------------------
### getVarCov.R ends here
