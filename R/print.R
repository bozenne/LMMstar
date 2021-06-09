### print.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: Jun  8 2021 (10:04) 
##           By: Brice Ozenne
##     Update #: 42
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * print.lmm (code)
##' @export
print.lmm <- function(x, ...){

    param.mu <- x$param$value[x$param$type=="mu"]
    param.sigma <- x$param$value[x$param$type=="sigma"]
    param.k <- x$param$value[x$param$type=="k"]
    param.rho <- x$param$value[x$param$type=="rho"]
    structure <- x$structure
    logLik <- stats::logLik(x)
    nobs <- stats::nobs(x)
    
    ## type of model
    if(length(param.rho) == 0){
        if(length(c(param.sigma,param.k))==1){
            cat("  Univariate Gaussian Model \n \n")
        }else{
            cat("  Univariate Gaussian Model with heterogeneous residual variance \n \n")
        }
    }else{
        if(structure=="UN"){
            cat("  Multivariate Gaussian Model with an unstructured covariance matrix \n \n")
        }else if(structure=="CS"){
            cat("  Multivariate Gaussian Model with a compound symmetry covariance matrix \n \n")
        }
    }

    ## dataset
    cat("data           : ",nobs["obs"], " observations and distributed in ", nobs["cluster"], " clusters \n",  sep = "")

    ## log-likelihood
    cat("log-likelihood : ", as.double(logLik),"\n",sep="")

    ## parameters
    cat( "parameters     : ",length(param.mu)," mean (",paste0(names(param.mu),collapse=" "),") \n",
         "                 ",length(c(param.sigma,param.k))," variance (",paste0(names(c(param.sigma,param.k)),collapse=" "),") \n",sep="")
    if(length(param.rho)>0){
        cat("                 ",length(param.rho)," correlation (",paste0(names(param.rho),collapse=" "),") \n",sep="")
    }
}


##----------------------------------------------------------------------
### print.R ends here
