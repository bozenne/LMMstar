### print.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: Nov 12 2022 (16:15) 
##           By: Brice Ozenne
##     Update #: 199
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
            cat("\t\tLinear regression \n")
        }else{
            cat("\t\tLinear regression with heterogeneous residual variance \n")
        }
    }else{
        if(structure$type=="UN"){
            if(is.na(structure$name$strata)){
                txt.strata <- "an"
            }else{
                txt.strata <- "a stratified"
            }
            cat("\t\tLinear Mixed Model with ",txt.strata," unstructured covariance matrix \n", sep = "")
        }else if(structure$type=="CS"){
            if(is.na(structure$name$strata)){
                txt.strata <- "a"
            }else{
                txt.strata <- "a stratified"
            }
            if(!structure$block){
                cat("\t\tLinear Mixed Model with ",txt.strata," compound symmetry covariance matrix \n", sep = "")
            }else if(structure$heterogeneous){
                cat("\t\tLinear Mixed Model with ",txt.strata," block unstructured covariance matrix \n", sep = "")
            }else{
                cat("\t\tLinear Mixed Model with ",txt.strata," block compound symmetry covariance matrix \n", sep = "")
            }
        }else if(structure$type=="TOEPLITZ"){
            if(is.na(structure$name$strata)){
                txt.strata <- "a"
            }else{
                txt.strata <- "a stratified"
            }
            if(!structure$block){
                cat("\t\tLinear Mixed Model with ",txt.strata," Toeplitz covariance matrix \n", sep = "")
            }else if(structure$heterogeneous == "UN"){
                cat("\t\tLinear Mixed Model with ",txt.strata," unstructured covariance matrix with constant subdiagonal \n", sep = "")
            }else if(structure$heterogeneous == "LAG"){
                cat("\t\tLinear Mixed Model with ",txt.strata," block Toeplitz covariance matrix \n", sep = "")
            }else if(structure$heterogeneous == "CS"){
                cat("\t\tLinear Mixed Model with ",txt.strata," block compound symmetry covariance matrix with specific subdiagonal \n", sep = "")
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
    ls.printparam <- list(c("parameter",": "),
                          c("","  "),
                          c("","  "))
    iParam <- 1
    if(length(param.mu)>0){
        M.print <- rbind(M.print,
                         cbind(ls.printparam[[iParam]][1],ls.printparam[[iParam]][2],
                               paste(length(param.mu)," mean (",paste0(names(param.mu),collapse=" "),")", sep="")))
        iParam <- iParam + 1
    }
    if(length(c(param.sigma,param.k))>0){
        M.print <- rbind(M.print,
                         cbind(ls.printparam[[iParam]][1],ls.printparam[[iParam]][2],
                               paste(length(c(param.sigma,param.k))," variance (",paste0(names(c(param.sigma,param.k)),collapse=" "),")", sep="")))
        iParam <- iParam + 1
    }
    if(length(param.rho)>0){
        M.print <- rbind(M.print,
                         cbind(ls.printparam[[iParam]][1],ls.printparam[[iParam]][2],
                               paste(length(param.rho)," correlation (",paste0(names(param.rho),collapse=" "),")", sep="")))
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


## * print.confint_lmm
##' @export
print.confint_lmm <- function(x, digits = 3, detail = FALSE, ...){

    ## ** display table
    print(as.data.frame(x), digits = digits, ...)

    ## ** caption
    message.backtransform <- attr(x,"backtransform")
    test.backtransform <- !is.null(message.backtransform) && any(!is.na(message.backtransform$FUN))

    adj.method <- c(stats::p.adjust.methods,"single-step", "Westfall", "Shaffer", "free", "single-step2")
    method.multiplicity <- attr(x,"method")
    test.multiplicity <- !is.null(method.multiplicity) && any(method.multiplicity %in% setdiff(adj.method,"none"))

    if(detail>0 && (test.multiplicity || test.backtransform)){

        cat("\n   Note: ")
        add.space <- ""

        if(test.multiplicity){
            if(!is.null(message.backtransform)){
                names_x <- intersect(names(x),names(message.backtransform))
            }else{
                names_x <- names(x)
            }
            if(NROW(x)==1){
                short2text <- stats::setNames(c("confidence interval","confidence interval","p-value"),c("lower","upper","p.value"))
                txt <- unique(short2text[intersect(names(short2text),names_x)])
            }else{
                short2text <- stats::setNames(c("confidence intervals","confidence intervals","p-values"),c("lower","upper","p.value"))
                txt <- unique(short2text[intersect(names(short2text),names_x)])
            }
            if("bonferroni" %in% method.multiplicity){
                txt.method <- "Bonferroni"
            }else if("single-step" %in% method.multiplicity || "single-step2" %in% method.multiplicity){
                txt.method <- "max-test adjustment"
            }else{
                txt.method <- method.multiplicity[method.multiplicity %in% adj.method]
            }
            add.space <- "         "
            cat(paste(txt,collapse = ", ")," have been adjusted for multiplicity using ",txt.method,". \n",sep="")
        }
        
        if(test.backtransform){
            message.backtransform <- message.backtransform[!is.na(message.backtransform$FUN),,drop=FALSE]

            if(any(message.backtransform[,setdiff(names(message.backtransform), "FUN")] == FALSE)){
                warning("Could not back-transform everything.\n")
            }

            if(NROW(x)==1){
                short2text <- stats::setNames(c("estimate","standard error","confidence interval","confidence interval"),c("estimate","se","lower","upper"))
                txt <- unique(short2text[intersect(names(short2text),intersect(names(x),names(message.backtransform)))])
            }else{
                short2text <- stats::setNames(c("estimates","standard errors","confidence intervals","confidence intervals"),c("estimate","se","lower","upper"))
                txt <- unique(short2text[intersect(names(short2text),intersect(names(x),names(message.backtransform)))])
            }
            cat(add.space,paste(txt,collapse = ", ")," have been back-transformed",sep="")
            if(detail>=0.5){
                cat(" (",paste0(paste(rownames(message.backtransform),collapse = "/")," parameters with ",paste(message.backtransform$FUN,collapse="/")),"). \n", sep ="")
            }
            cat("\n")
        }
    }
    return(invisible(NULL))
}

## * print.Wald_lmm
##' @export
print.Wald_lmm <- function(x, ...){
    dots <- list(...)
    dots$print <- c(1,0)
    return(do.call(summary, c(list(object = x, legend = FALSE), dots)))
}

## * print.LRT_lmm
##' @export
print.LRT_lmm <- function(x, ...){

    ## ** display
    cat("\t\tLikelihood ratio test \n\n")
    out <- as.data.frame(x)
    out$null <- NULL
    rownames(out) <- ""
    print(as.data.frame(out))
    cat("\n")

    ## ** export
    return(invisible(NULL))
}



## * print.mlmm
##' @export
print.mlmm <- function(x, ...){

    print(x$model)

    return(invisible(NULL))
}

## * print.partialCor
##' @export
print.partialCor <- function(x, digits = 3, ...){
    out <- do.call("print.confint_lmm", c(list(x, detail = FALSE, digits = digits), ...))
    return(invisible(NULL))
}

## * print.resample
##' @export
print.resample <- function(x, digits = 3, ...){
    args <- attr(x,"args")
    n.sample <- attr(x,"n.sample")
    
    base::print.data.frame(x, digits = digits)
    if(n.sample!=args$n.sample){
        cat(paste0("(based on ",n.sample," samples - ",round((1-n.sample/args$n.sample)*100, digits = digits),"% failed) \n"))
    }
    return(invisible(NULL))
}

##----------------------------------------------------------------------
### print.R ends here
