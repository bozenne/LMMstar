### hessian.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (20:56) 
## Version: 
## Last-Updated: mar 22 2021 (22:36) 
##           By: Brice Ozenne
##     Update #: 24
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .hessian
.hessian <- function(X, precision,
                     index.variance, index.cluster, indiv, REML){

    ## ** prepare
    n.obs <- length(index.cluster)
    n.cluster <- length(index.variance)
    n.allcoef <- NCOL(X)
    name.allcoef <- colnames(X)
    Hessian <- array(NA, dim = c(n.cluster, n.allcoef, n.allcoef),
                      dimnames = list(NULL, name.allcoef, name.allcoef))
    

    ## ** compute score
    for(iId in 1:n.cluster){ ## iId <- 7
        iX <- X[index.cluster==iId,,drop=FALSE]
        iOmega <- precision[[index.variance[iId]]]
        Hessian[iId,,] <- - t(iX) %*% iOmega %*% iX
    }

    ## ** export
    if(indiv){
        return(Hessian)
    }else{
        return(apply(Hessian, MARGIN = 2:3,FUN = sum))
    }

}


##----------------------------------------------------------------------
### hessian.R ends here
