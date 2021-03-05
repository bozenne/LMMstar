### score.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:59) 
## Version: 
## Last-Updated: mar  5 2021 (21:25) 
##           By: Brice Ozenne
##     Update #: 11
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .score
.score <- function(Y, X, beta, sigma, k, rho, precision,
                   index.variance, index.cluster, indiv, REML){

    ## ** prepare
    n.obs <- length(index.cluster)
    n.cluster <- length(index.variance)
    n.allcoef <- NCOL(X)
    name.allcoef <- colnames(X)
    Score <- matrix(NA, nrow = n.cluster, ncol = n.allcoef,
                    dimnames = list(NULL, name.allcoef))
    
    residuals <- Y - X %*% beta

    ## ** compute score
    for(iId in 1:n.cluster){ ## iId <- 7
        iResidual <- residuals[index.cluster==iId,,drop=FALSE]
        iX <- X[index.cluster==iId,,drop=FALSE]
        iOmega <- precision[[index.variance[iId]]]
        Score[iId,] <- t(iX) %*% iOmega %*% iResidual
    }

    ## ** export
    if(indiv){
        return(Score)
    }else{
        return(rowSums(Score))
    }

}


##----------------------------------------------------------------------
### score.R ends here
