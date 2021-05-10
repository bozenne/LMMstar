### print.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: May 10 2021 (19:29) 
##           By: Brice Ozenne
##     Update #: 36
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

    ## type of model
    if(length(x$param$rho) > 0){
        if(length(c(x$param$sigma,x$param$k))==1){
            cat("  Linear model \n \n")
        }else{
            cat("  Linear model with heterogeneous residual variance \n \n")
        }
    }else{
        if(x$structure=="UN"){
            cat("  Linear mixed effect model with an unstructured covariance matrix \n \n")
        }else if(x$structure=="CS"){
            cat("  Linear mixed effect model with a compound symmetry covariance matrix \n \n")
        }
    }

    ## dataset
    cat("data           : ",sum(x$design$cluster$nobs), " observations and distributed in ", x$design$cluster$n, " clusters \n",  sep = "")

    ## log-likelihood
    cat("log-likelihood : ", as.double(x$logLik),"\n",sep="")

    ## parameters
    cat( "parameters     : ",length(x$param$mu)," mean (",paste0(names(x$param$mu),collapse=" "),") \n",
         "                 ",length(x$param[c("sigma","k")])," variance (",paste0(names(c(x$param$sigma,x$param$k)),collapse=" "),") \n",sep="")
    if(length(x$param$rho)>0){
        cat("                 ",length(x$param$rho)," correlation (",paste0(names(x$param$rho),collapse=" "),") \n",sep="")
    }
}


##----------------------------------------------------------------------
### print.R ends here
