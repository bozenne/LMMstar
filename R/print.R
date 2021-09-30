### print.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: sep 30 2021 (15:54) 
##           By: Brice Ozenne
##     Update #: 85
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
    structure <- x$design$vcov
    logLik <- stats::logLik(x)
    nobs <- stats::nobs(x)

    ## ** prepare
    M.print <- NULL
    
    ## ** type of model
    if(length(param.rho) == 0){
        if(length(c(param.sigma,param.k))==1){
            cat("     Linear regression \n")
        }else{
            cat("     Linear regression with heterogeneous residual variance \n")
        }
    }else{
        if(structure$type=="UN"){
            cat("     Linear Mixed Model with an unstructured covariance matrix \n")
        }else if(structure$type=="CS"){
            cat("     Linear Mixed Model with a compound symmetry covariance matrix \n")
        }
    }

    ## ** outcome/cluster/time
    txt.var <- "outcome"
    value.var <- x$outcome$var
    if(!is.null(x$cluster$var)){
        txt.var <- c(txt.var,"cluster")
        value.var <- c(value.var,x$cluster$var)
    }
    if(!is.null(x$time$var)){
        txt.var <- c(txt.var,"time")
        value.var <- c(value.var,x$time$var)
    }
    Ctxt.var <- paste(txt.var,collapse="/")

    M.print <- rbind(M.print,
                     cbind(Ctxt.var,": ",paste(value.var, collapse="/")))

    ## ** dataset
    M.print <- rbind(M.print,
                     cbind("data",": ",paste(nobs["obs"], " observations and distributed in ", nobs["cluster"], " clusters",sep="")))

    ## ** parameters
    M.print <- rbind(M.print,
                     cbind("parameters",": ",paste(length(param.mu)," mean (",paste0(names(param.mu),collapse=" "),")", sep="")))
    M.print <- rbind(M.print,
                     cbind("","  ",paste(length(c(param.sigma,param.k))," variance (",paste0(names(c(param.sigma,param.k)),collapse=" "),")", sep="")))
    if(length(param.rho)>0){
    M.print <- rbind(M.print,
                     cbind("","  ",paste(length(param.rho)," correlation (",paste0(names(param.rho),collapse=" "),")", sep="")))
    }

    ## ** log-likelihood
    M.print <- rbind(M.print,
                     cbind("log-likelihood",": ",as.double(logLik)))

    ## ** optimisation
    if(!is.null(x$opt)){
        M.print <- rbind(M.print,
                         cbind("convergence",": ",paste0(x$opt$cv," (",x$opt$n.iter," iterations)")))
    }

    ## ** print
    need.blank <- max(nchar(M.print[,1]))-nchar(M.print[,1])
    add.blank <- sapply(need.blank, function(iN){paste(rep(" ",iN),collapse="")})
    M.print[,1] <- paste0(M.print[,1],add.blank)
    txt.print <- unname(sapply(apply(M.print,1,paste,collapse=""),paste,"\n"))
    cat("\n ")
    cat(txt.print)
    return(invisible(NULL))
}


##----------------------------------------------------------------------
### print.R ends here
