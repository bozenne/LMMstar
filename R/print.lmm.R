### print.lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: jun 24 2022 (17:43) 
##           By: Brice Ozenne
##     Update #: 110
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

    ## ** extract from object
    param.value <- x$param
    param.type <- stats::setNames(x$design$param$type, x$design$param$name)
    
    param.mu <- param.value[names(which(param.type=="mu"))]
    param.sigma <- param.value[names(which(param.type=="sigma"))]
    param.k <- param.value[names(which(param.type=="k"))]
    param.rho <- param.value[names(which(param.type=="rho"))]
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
            if(is.na(structure$name$strata)){
                txt.strata <- "an"
            }else{
                txt.strata <- "a stratified"
            }
            cat("     Linear Mixed Model with ",txt.strata," unstructured covariance matrix \n", sep = "")
        }else if(structure$type=="CS"){
            if(is.na(structure$name$strata)){
                txt.strata <- "a"
            }else{
                txt.strata <- "a stratified"
            }
            if(all(is.na(structure$name$cor[[1]]))){
                cat("     Linear Mixed Model with ",txt.strata," compound symmetry covariance matrix \n", sep = "")
            }else if(structure$heterogeneous){
                cat("     Linear Mixed Model with ",txt.strata," block unstructured covariance matrix \n", sep = "")
            }else{
                cat("     Linear Mixed Model with ",txt.strata," block compound symmetry covariance matrix \n", sep = "")
            }
        }
    }

    ## ** outcome/cluster/time
    txt.var <- "outcome"
    value.var <- x$outcome$var
    if(!is.null(x$cluster$var)){
        txt.var <- c(txt.var,"cluster")
        value.var <- c(value.var,x$cluster$var)
    }
    if(!all(is.na(x$time$var))){
        txt.var <- c(txt.var,"time")
        value.var <- c(value.var,paste(x$time$var, collapse = ", "))
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
    if(x$method.fit=="ML"){
        M.print <- rbind(M.print,
                         cbind("log-likelihood",": ",as.double(logLik)))
    }else if(x$method.fit=="REML"){
        M.print <- rbind(M.print,
                         cbind("log-restr.likelihood",": ",as.double(logLik)))
    }

    ## ** optimisation
    if(x$opt$name!="gls"){
        if(!is.na(x$opt$n.iter)){
            M.print <- rbind(M.print,
                             cbind("convergence",": ",paste0(x$opt$cv>0," (",x$opt$n.iter," iterations)")))
        }else if(!is.null(attr(x$opt$n.iter,"eval"))){
            M.print <- rbind(M.print,
                             cbind("convergence",": ",paste0(x$opt$cv>0," (evaluations: ",attr(x$opt$n.iter,"eval")["logLik"]," likelihood, ",attr(x$opt$n.iter,"eval")["score"]," score)")))
        }else{
            M.print <- rbind(M.print,
                             cbind("convergence",": ",x$opt$cv>0))
        }
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
### print.lmm.R ends here
