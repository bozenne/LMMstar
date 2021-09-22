### precompute.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 22 2021 (13:47) 
## Version: 
## Last-Updated: sep 22 2021 (13:48) 
##           By: Brice Ozenne
##     Update #: 2
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .precomputeXX
## Precompute square of the design matrix
.precomputeXX <- function(X, pattern, pattern.time, pattern.cluster, index.cluster){
    p <- NCOL(X)
    n.pattern <- length(pattern)
    n.time <- lapply(pattern.time,length)

    out <- list(pattern = stats::setNames(lapply(pattern, function(iPattern){array(0, dim = c(n.time[[iPattern]],n.time[[iPattern]],p*(p+1)/2))}), pattern),
                key = matrix(as.numeric(NA),nrow=p,ncol=p,dimnames=list(colnames(X),colnames(X))),
                Xpattern = stats::setNames(vector(mode = "list", length = n.pattern),pattern))

    ## key
    out$key[lower.tri(out$key,diag = TRUE)] <- 1:sum(lower.tri(out$key,diag = TRUE))
    out$key[upper.tri(out$key)] <- t(out$key)[upper.tri(out$key)]

    ## fill matrix
    for(iPattern in pattern){ ## iPattern <- pattern[1]
        iTime <- length(pattern.time[[iPattern]])

        if(iTime==1){
            out$Xpattern[[iPattern]] <- do.call(rbind,lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){X[iIndex,,drop=FALSE]}))
            iX.summary <- crossprod(out$Xpattern[[iPattern]])
            ## out$key[lower.tri(out$key,diag = TRUE)]
            out$pattern[[iPattern]][1,1,] <- iX.summary[lower.tri(iX.summary, diag = TRUE)]
        }else{
            out$Xpattern[[iPattern]] <- array(unlist(lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){X[iIndex,,drop=FALSE]})),
                                              dim = c(iTime,NCOL(X),length(index.cluster[pattern.cluster[[iPattern]]])))

            for(iCol1 in 1:p){ ## iCol1 <- 1
                for(iCol2 in 1:iCol1){ ## iCol2 <- 2
                    ## for(iId in pattern.cluster[[iPattern]]){
                    ##     out$pattern[[iPattern]][,,out$key[iCol1,iCol2]] <- out$pattern[[iPattern]][,,out$key[iCol1,iCol2]] + tcrossprod(X[index.cluster[[iId]],iCol1,drop=FALSE],X[index.cluster[[iId]],iCol2,drop=FALSE])
                    ## }
                    out$pattern[[iPattern]][,,out$key[iCol1,iCol2]] <- tcrossprod(out$Xpattern[[iPattern]][,iCol1,],out$Xpattern[[iPattern]][,iCol2,])
                }
            }
        }
    }
    return(out)
}

## * .precomputeXR
## Precompute design matrix times residuals
.precomputeXR <- function(X, residuals, pattern, pattern.time, pattern.cluster, index.cluster){
    p <- NCOL(X[[1]])
    name.mucoef <- colnames(X)
    n.pattern <- length(pattern)
    n.time <- lapply(pattern.time,length)

    out <- stats::setNames(lapply(pattern, function(iPattern){array(0, dim = c(n.time[[iPattern]], n.time[[iPattern]], ncol = p), dimnames = list(NULL,NULL,name.mucoef))}), pattern)

    for(iPattern in pattern){ ## iPattern <- pattern[1]

        iTime <- length(pattern.time[[iPattern]])
        iResiduals <- do.call(cbind, lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){residuals[iIndex,,drop=FALSE]}))
        iX <- X[[iPattern]]

        if(iTime == 1){
            out[[iPattern]][1,1,] <- iResiduals %*% iX
        }else{
            for(iCol in 1:p){ ## iCol1 <- 1
                ## for(iId in 1:length(pattern.cluster[[iPattern]])){ ## iId <- 1
                ##     out[[iPattern]][,,iCol] <- out[[iPattern]][,,iCol] + tcrossprod(iX[[iId]][,iCol,drop=FALSE],iResiduals[[iId]])
                ## }
                out[[iPattern]][,,iCol] <- tcrossprod(iX[,iCol,], iResiduals)
            }
        }
    }

    return(out)
}

## * .precomputeRR
## Precompute square of the residuals
.precomputeRR <- function(residuals, pattern.time, pattern, pattern.cluster, index.cluster){

    n.pattern <- length(pattern)
    n.time <- stats::setNames(lapply(pattern.time,length), pattern)
    out <- stats::setNames(lapply(pattern, function(iPattern){matrix(0, nrow = n.time[[iPattern]], ncol = n.time[[iPattern]])}), pattern)

    for(iPattern in pattern){ ## iPattern <- pattern[1]
        
        ## for(iId in pattern.cluster[[iPattern]]){
        ##     out[[iPattern]] <- out[[iPattern]] + tcrossprod(residuals[index.cluster[[iId]],,drop=FALSE])
        ## }
        out[[iPattern]] <- tcrossprod(do.call(cbind,lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){residuals[iIndex,,drop=FALSE]})))
    }
    return(out)
}


##----------------------------------------------------------------------
### precompute.R ends here
