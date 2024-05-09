### print.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: May  9 2024 (13:18) 
##           By: Brice Ozenne
##     Update #: 252
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
    }else if(inherits(structure,"RE")){
        structure.ranef <- structure$ranef
        if(structure.ranef$crossed==FALSE && structure.ranef$nested==FALSE){
            cat("\t\tLinear Mixed Model with a random intercept \n", sep = "")
        }else if(structure.ranef$crossed==FALSE && structure.ranef$nested==TRUE){
            cat("\t\tLinear Mixed Model with nested random intercepts \n", sep = "")
        }else if(structure.ranef$crossed==TRUE && structure.ranef$nested==FALSE){
            cat("\t\tLinear Mixed Model with cross random intercepts \n", sep = "")
        }else{
            cat("\t\tLinear Mixed Model with random effects \n", sep = "")
        }        
    }else{
        if(inherits(structure,"UN")){
            if(is.na(structure$name$strata)){
                txt.strata <- "an"
            }else{
                txt.strata <- "a stratified"
            }
            cat("\t\tLinear Mixed Model with ",txt.strata," unstructured covariance matrix \n", sep = "")
        }else if(inherits(structure,"CS")){
            if(is.na(structure$name$strata)){
                txt.strata <- "a"
            }else{
                txt.strata <- "a stratified"
            }
            if(all(is.na(structure$name$cor))){
                cat("\t\tLinear Mixed Model with ",txt.strata," compound symmetry covariance matrix \n", sep = "")
            }else if(structure$type == "heterogeneous"){
                cat("\t\tLinear Mixed Model with ",txt.strata," block unstructured covariance matrix \n", sep = "")
            }else if(structure$type == "homogeneous"){
                cat("\t\tLinear Mixed Model with ",txt.strata," block compound symmetry covariance matrix \n", sep = "")
            }else if(structure$type == "heterogeneous0"){
                cat("\t\tLinear Mixed Model with ",txt.strata," crossed unstructured covariance matrix \n", sep = "")
            }else if(structure$type == "homogeneous0"){
                cat("\t\tLinear Mixed Model with ",txt.strata," crossed compound symmetry covariance matrix \n", sep = "")
            }
        }else if(inherits(structure,"TOEPLITZ")){            
            if(is.na(structure$name$strata)){
                txt.strata <- "a"
            }else{
                txt.strata <- "a stratified"
            }
            if(all(is.na(structure$name$cor))){
                cat("\t\tLinear Mixed Model with ",txt.strata," Toeplitz covariance matrix \n", sep = "")
            }else if(structure$type == "heterogeneous"){
                cat("\t\tLinear Mixed Model with ",txt.strata," unstructured covariance matrix with constant subdiagonal \n", sep = "")                
            }else if(tolower(structure$type) == "lag"){
                cat("\t\tLinear Mixed Model with ",txt.strata," block Toeplitz covariance matrix \n", sep = "")
            }else if(structure$type == "homogeneous"){
                cat("\t\tLinear Mixed Model with ",txt.strata," block compound symmetry covariance matrix with specific subdiagonal \n", sep = "")
            }
        }else if(inherits(structure,"CUSTOM")){
            cat("\t\tLinear Mixed Model with user-defined covariance matrix \n", sep = "")
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
                     cbind("data",": ",paste(nobs["obs"], " observations from ", nobs["cluster"], " clusters",sep="")))

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
    if(x$args$method.fit=="ML"){
        M.print <- rbind(M.print,
                         cbind("log-likelihood",": ",as.double(logLik)))
    }else if(x$args$method.fit=="REML"){
        M.print <- rbind(M.print,
                         cbind("log-restr.likelihood",": ",as.double(logLik)))
    }

    ## ** optimisation
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


## * print.Wald_lmm
##' @export
print.effect_lmm <- function(x, ...){
    dots <- list(...)
    dots$print <- c(0,0.5)
    return(do.call(summary, c(list(object = x, legend = FALSE), dots)))
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

## * print.residuals_lmm
#' @export
print.residuals_lmm <- function(x, ...){
    x.print <- x
    attr(x.print, "args") <- NULL
    class(x.print) <- setdiff(class(x.print),"residuals_lmm")
    print(x.print)
    
}
## * print.summarize
#' @export
print.summarize <- function(x,...){
    ## remove duplicated values
    if(length(unique(x$outcome))==1){
        x$outcome <- NULL
    }else{
        x$outcome[duplicated(x$outcome)] <- ""
    }
    if("pc.missing" %in% eval(attr(x,"call")$columns) == FALSE){
        x$pc.missing <- NULL
    }
    name.X <-  attr(x,"name.X")
    if(length(name.X)>0){
        for(iX in name.X){ ## iX <- name.X[2]
            iX.value <- as.numeric(as.factor(x[[iX]]))
            iLevels <- cumsum(iX.value!=c(0,iX.value[-length(iX.value)]))
            if(any(duplicated(iLevels))){
                x[[iX]] <- as.character(x[[iX]])
                x[[iX]][duplicated(iLevels)] <- ""
            }
        }
    }

    if(!is.null(attr(x,"digits")) && ("digits" %in% names(list(...)) == FALSE)){
        print(as.data.frame(x), digits = attr(x,"digits"), ...)
    }else{
        print(as.data.frame(x), ...)
    }
    if(!is.null(attr(x,"correlation"))){
        cat("\n Pearson's correlation: \n")
        ls.cor <- attr(x,"correlation")
        if(length(ls.cor)==1){ ## outcome
            ls.cor <- ls.cor[[1]]
            if(length(ls.cor)==1){ ## group
                ls.cor <- ls.cor[[1]]
            }
        }else{
            ls.cor <- unlist(ls.cor, recursive = FALSE)
        }
        print(ls.cor, ...)
    }
    return(invisible(NULL))
}

## * print.summarizeNA
#' @export
print.summarizeNA <- function(x,...){

    ## total <- c(list(frequency = NA), list(missing.pattern = "any"), as.list(colSums(sweep(x[,-(1:2)], FUN = "*", MARGIN = 1, STATS = x$frequency))))
    newx <- x
    newnames <- attr(x, "args")$newnames

    if(newnames[1] %in% names(x)){
        newx[[newnames[1]]][duplicated(x[[newnames[1]]])] <- ""
    }

    ## newx <- rbind(x, as.data.frame(total))
    toprint <- format.data.frame(newx, digits = NULL, na.encode = FALSE)
    toprint[is.na(newx)] <- NA
    print(toprint, na.print="", quote=FALSE, row.names = FALSE)

    
    return(NULL)
    
}
##----------------------------------------------------------------------
### print.R ends here
