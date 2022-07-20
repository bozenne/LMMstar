### rbind.anova_lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:51) 
## Version: 
## Last-Updated: jul 20 2022 (15:11) 
##           By: Brice Ozenne
##     Update #: 250
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * rbind.anova_lmm (documentation)
##' @title Linear Hypothesis Testing Across Linear Mixed Models
##' @description Linear hypothesis testing accross linear mixed model.
##'
##' @param model a \code{anova_lmm} object (output of \code{anova} applied to a \code{lmm} object)
##' @param ...  possibly other \code{anova_lmm} objects
##' @param effects [numeric matrix] matrix indicating how to combine the left-hand side of the hypotheses. By default identity matrix.
##' @param rhs [numeric vector] the right hand side of the hypothesis. Should have the same length as the number of row of argument \code{effects}.
##' @param name [character vector or NULL] character used to identify each model in the output.
##' By default, use the name of the outcome of the model.
##' @param sep [character] character used to separate the outcome and the covariate when naming the tests.
##'
##' @examples
##' ## simulate data
##' set.seed(10)
##' dL <- sampleRem(1e2, n.times = 3, format = "long")
##'
##' ## estimate mixed models
##' e.lmm1 <- lmm(Y ~ X1+X2+X3, repetition = ~visit|id, data = dL,
##'               structure = "CS", df = FALSE)
##' e.lmm2 <- lmm(Y ~ X1+X8+X9, repetition = ~visit|id, data = dL,
##'               structure = "CS", df = FALSE)
##'
##' ## select null hypotheses
##' AAA <- anova(e.lmm1, ci = TRUE, effect = c("X1|X2,X3"="X1=0","X2|X1,X3"="X2=0"))
##' BBB <- anova(e.lmm2, ci = TRUE, effect = c("X1|X8,X9"="X1=0"))
##'
##' ## combine
##' ZZZ <- rbind(AAA,BBB)
##' summary(ZZZ)

## * rbind.anova_lmm (code)
##' @export
rbind.anova_lmm <- function(model, ..., effects = NULL, rhs = NULL, name = NULL, sep = ": "){
    default <- LMMstar.options()
    call <- match.call()
    
    ## ** Check user input
    dots <- list(...)
    if(length(dots)==0){
        return(model) ## nothing to combine
    }else if(any(sapply(dots,inherits,"anova_lmm")==FALSE)){
        stop("Extra arguments should inherit from anova_lmm. \n")
    }
    ls.object <- c(list(model),dots)
    n.object <- length(ls.object)
    if(!is.null(name) && length(name)!=n.object){
        stop("Argument \'name\' should have length ",n.object,", i.e. as many elements as anova objects. \n")
    }
    
    table.args <- cbind(do.call(rbind,lapply(ls.object,"[[","args")),
                        do.call(rbind,lapply(ls.object, function(iO){data.frame(n.test = NROW(iO$univariate),
                                                                                iO$object[c("outcome","method.fit","type.information")])
                        })))

    Utype <- unique(unlist(table.args$type))
    newtable.args <- data.frame(type = ifelse(length(Utype)>1,"all",Utype),
                                table.args[1,c("robust","df","ci","transform.sigma","transform.k","transform.rho","transform.names")]
                                )
    newtable.args$n.test <- sum(table.args$n.test)
    if(length(unique(table.args$outcome))==1){
        newtable.args$outcome <- table.args$outcome[1]
    }
    if(length(unique(table.args$method.fit))==1){
        newtable.args$method.fit <- table.args$method.fit[1]
    }
    if(length(unique(table.args$type.information))==1){
        newtable.args$type.information <- table.args$type.information[1]
    }
    
    if(any(table.args$ci==FALSE)){
        stop("All argument should contain a \"glht\" object, i.e. call anova with argument ci=TRUE. \n")
    }
    if(any(newtable.args$robust!=table.args$robust[-1])){
        stop("Element \'robust\' should take the same value for all objects. \n")
    }
    if(any(newtable.args$df!=table.args$df[-1])){
        stop("Element \'df\' should take the same value for all objects. \n")
    }
    if(!all(is.na(table.args$transform.sigma)) && any(newtable.args$transform.sigma!=table.args$transform.sigma[-1])){
        stop("Element \'transform.sigma\' should take the same value for all objects. \n")
    }
    if(!all(is.na(table.args$transform.k)) && any(newtable.args$transform.k!=table.args$transform.k[-1])){
        stop("Element \'transform.k\' should take the same value for all objects. \n")
    }
    if(!all(is.na(table.args$transform.rho)) && any(newtable.args$transform.rho!=table.args$transform.rho[-1])){
        stop("Element \'transform.rho\' should take the same value for all objects. \n")
    }
    if(any(newtable.args$transform.names!=table.args$transform.names[-1])){
        stop("Element \'transform.names\' should take the same value for all objects. \n")
    }
    ls.glht <- unname(unlist(unlist(lapply(ls.object,"[[","glht"),recursive = FALSE),recursive = FALSE))
    alternative <- ls.glht[[1]]$alternative
    if(any(alternative != sapply(ls.glht,"[[","alternative")[-1])){
        stop("Element \'alternative\' should take the same value for all glht objects. \n")
    }

    if(is.null(name)){
        name.unitest <- unlist(lapply(ls.object, function(iO){rownames(iO$univariate)}))
        if(any(duplicated(name.unitest))){
            name.unitest <- unlist(lapply(ls.object, function(iO){paste0(iO$object$outcome, sep, rownames(iO$univariate))}))

            if(is.na(sep) || any(duplicated(name.unitest))){
                stop("Univariate test should have distinct names. \n",
                     "Consider naming them when calling anova, e.g. anova(object, effect = c(\"myname1\"=\"X1=0\",\"myname2\"=\"X2=0\"))")
            }else{
                for(iO in 1:n.object){
                    rownames(ls.object[[iO]]$univariate) <- paste0(ls.object[[iO]]$object$outcome, sep, rownames(ls.object[[iO]]$univariate))
                }
            }
        }
    }else{
        name.unitest <- unlist(lapply(1:n.object, function(iO){paste0(name[iO],sep,rownames(ls.object[[iO]]$univariate))}))
        for(iO in 1:n.object){
            rownames(ls.object[[iO]]$univariate) <- paste0(name[iO],sep,rownames(ls.object[[iO]]$univariate))
        }
        if(any(duplicated(name.unitest))){
            stop("Univariate test should have distinct names. \n",
                 "Consider naming them when calling anova, e.g. anova(object, effect = c(\"myname1\"=\"X1=0\",\"myname2\"=\"X2=0\"))")
        }
    }

    contrast <- effects
    if(!is.null(contrast)){

        browser()
        contrMat(n=c(1,1,1), type='Tukey')
        if(NCOL(contrast)!=newtable.args$n.test){
            stop("Incorrect contrast matrix: should have ",newtable.args$n.test," columns.\n",
                 "(one for each univariate test) \n")
        }
        if(is.null(colnames(contrast))){
            colnames(contrast) <- name.unitest
        }else if(!identical(sort(colnames(contrast)), name.unitest)){
            stop("Incorrect column names for argument contrast.\n",
                 "Should match \"",paste(name.unitest, collapse="\" \""),"\".\n")
        }else{
            contrast <- contrast[,name.unitest,drop=FALSE]
        }
    }else{
        contrast <- diag(1, nrow = newtable.args$n.test, ncol = newtable.args$n.test)
        dimnames(contrast) <- list(name.unitest, name.unitest)
    }
    if(!is.null(rhs)){
        if(length(rhs)!=newtable.args$n.test){
            stop("Incorrect rhs: should have ",newtable.args$n.test," values.\n",
                 "(one for each univariate test) \n")
        }
    }else{
        rhs <- stats::setNames(rep(0, newtable.args$n.test), name.unitest)
    }


    cluster.var <- ls.object[[1]]$object$cluster.var
    if(any(cluster.var!=unlist(lapply(ls.object[-1], function(iO){iO$object$cluster.var})))){
        stop("Cluster variable differs among objects. \n")
    }
    
    ## ** Try to find unique names

    ## *** for the model parameters
    ## find all names
    ls.testname <- lapply(ls.object, function(iO){iO$multivariate$test})
    vec.testname <- unlist(ls.testname)

    if(!is.null(name)){
        outcome <- name
    }else if(all(sapply(ls.testname, function(iO){all(is.character(iO))})) && all(duplicated(vec.testname)==FALSE)){
        outcome <- vec.testname
    }else{
        outcome <- sapply(ls.object, function(iO){iO$object$outcome})
        if(any(duplicated(outcome))){
            outcome <- c(deparse(call$model),
                         sapply(setdiff(which(names(call)==""),1), function(iIndex){deparse(call[[iIndex]])}))
        }
    }

    ## *** for the hypotheses
    gridTest <- do.call(rbind,lapply(1:n.object, function(iO){
        cbind(outcome = outcome[iO], ls.object[[iO]]$multivariate[,c("type","test"),drop=FALSE])
    }))

    if(is.na(sep) && all(duplicated(gridTest$test)==FALSE)){
        gridTest <- gridTest[,"test",drop=FALSE]
    }else if(!is.na(sep) && all(duplicated(outcome)==FALSE)){
        gridTest <- gridTest[,"outcome",drop=FALSE]
    }else if(all(duplicated(gridTest$test)==FALSE)){
        gridTest <- gridTest[,"test",drop=FALSE]
    }else if(all(duplicated(gridTest$type)==FALSE)){
        gridTest <- gridTest[,"type",drop=FALSE]
    }else{
        if(all(duplicated(outcome))){
            gridTest$outcome <- NULL
        }
        if(all(duplicated(outcome))){
            gridTest$type <- NULL
        }
        if(all(duplicated(test))){
            gridTest$test <- NULL
        }
    }
    col.nametest <- colnames(gridTest)
    name.test <- unique(as.character(interaction(gridTest[,col.nametest,drop=FALSE])))

    ## ** Extract elements from anova object
    ## *** univariate Wald test
    newtable.univariate <- do.call(rbind,lapply(1:n.object,function(iO){ ## iO <- 1
        iTable <- cbind(outcome = outcome[iO], ls.object[[iO]]$univariate)
        iTable$name.test <- factor(as.character(interaction(iTable[,col.nametest,drop=FALSE])),levels = name.test)
        return(iTable)
    }))
    newtable.univariate$lower <- NA
    newtable.univariate$upper <- NA
    newtable.univariate$p.value <- NA
    
    ## *** cluster
    ls.cluster <- lapply(ls.object, function(iO){iO$object$cluster})
    seq.cluster <- unique(unlist(ls.cluster))
    n.cluster <- length(seq.cluster)

    independence <- (cluster.var=="XXclusterXX" || all(duplicated(unlist(ls.cluster))==FALSE)) ## each individual only appear once
    if(all(sapply(ls.object, function(iO){attr(iO$args$robust,"call")})==FALSE)){
        newtable.args$robust <- independence==FALSE
    }

    ## *** estimate
    beta.estimate <- stats::setNames(newtable.univariate$estimate, name.unitest)
    if(independence){
        
        beta.vcov <- as.matrix(do.call(Matrix::bdiag,lapply(ls.object, "[[", "vcov")))
        
        if(all(diag(contrast)==1) && all(contrast[upper.tri(contrast, diag = FALSE)] == 0) && all(contrast[lower.tri(contrast, diag = FALSE)] == 0)){
            beta.df <- newtable.univariate$df
        }else{
            beta.df <- rep(Inf, length(beta.estimate))
            newtable.args$df <- FALSE
        }
        
    }else{
        if(any(sapply(ls.object, function(iO){is.null(iO$iid)}))){
            stop("Cannot evaluate the joint influence function of the variance-covariance parameters. \n",
                 "Consider using method.fit=\"ML\" when calling lmm. \n")
        }
        M.iid <- matrix(0, nrow = n.cluster, ncol = length(beta.estimate),
                        dimnames = list(seq.cluster, names(beta.estimate)))
        for(iO in 1:n.object){            
            M.iid[ls.cluster[[iO]],rownames(ls.object[[iO]]$univariate)] <- ls.object[[iO]]$iid
        }
        
        if(any(is.na(M.iid))){
            beta.vcov <- tcrossprod(sqrt(apply(M.iid^2, 2, sum, na.rm = TRUE))) * stats::cor(M.iid, use = "pairwise")
            ## usually better compared to formula 11.43 from chapter 11.4 of the book High-dimensional statistics by WAINWRIGHT
            ## iIDD0 <- M.iid/(1-mean(is.na(M.iid)))
            ## iIDD0[is.na(M.iid)] <- 0
            ## out$vcov <- crossprod(iIDD0) - mean(is.na(M.iid))*diag(diag(crossprod(iIDD0)))
            ## out$vcov - crossprod(M.iid)
        }else{
            beta.vcov <- crossprod(M.iid)
        }
        beta.df <- rep(Inf, length(beta.estimate))
        newtable.args$df <- FALSE
    }
    dimnames(beta.vcov) <- list(name.unitest,name.unitest)
    names(beta.df) <- name.unitest

    ## ** Combine elements

    ## *** multivariate tests
    outSimp <- .simplifyContrast(contrast, rhs) ## remove extra lines
    C.vcov.C <- outSimp$C %*% beta.vcov %*% t(outSimp$C)
    C.vcov.C_M1 <- try(solve(C.vcov.C), silent = TRUE)
    if(inherits(C.vcov.C_M1,"try-error")){
        multistat <- NA
        attr(multistat,"error") <- "\n  Could not invert the covariance matrix for the proposed contrast."
    }else{
        multistat <- as.double(t(outSimp$C %*% beta.estimate - outSimp$rhs) %*% C.vcov.C_M1 %*% (outSimp$C %*% beta.estimate - outSimp$rhs))/outSimp$dim 
    }

    Utype2 <- unique(unlist(lapply(ls.object, function(iO){iO$multivariate$type})))
    newtable.multivariate <- data.frame(type = ifelse(length(Utype2)>1,"all",Utype2),
                                        test = 1,
                                        null = paste(unlist(lapply(ls.object,function(iO){iO$multivariate$null})), collapse = ", "),
                                        statistic = multistat,
                                        df.num = newtable.args$n.test,
                                        df.denom = Inf,
                                        partial.r2 = NA,
                                        p.value = 1 - stats::pf(multistat, df1 = newtable.args$n.test, df2 = Inf))


    ## *** univariate
    newtable.univariate$se <- sqrt(diag(C.vcov.C))
    newtable.univariate$statistic <- newtable.univariate$estimate/newtable.univariate$se
    if(newtable.args$df){
        e.glht <- list(linfct = contrast, rhs = rhs, coef = beta.estimate, vcov = beta.vcov, df = ceiling(stats::median(beta.df)), alternative = alternative)
    }else{
        e.glht <- list(linfct = contrast, rhs = rhs, coef = beta.estimate, vcov = beta.vcov, df = FALSE, alternative = alternative)
    }
    class(e.glht) <- "glht"

    ## ** Export
    out <- list(multivariate = newtable.multivariate[,names(model$multivariate),drop=FALSE],
                univariate = newtable.univariate[,names(model$univariate),drop=FALSE],
                glht = list(all = list("1" = e.glht)),
                args = newtable.args,
                vcov = C.vcov.C)

    attr(out, "test") <- "Wald"
    class(out) <- append(c("manova_lmm","anova_lmm"),class(out))
    return(out)
}


## * rbind.manova_lmm (code)
##' @export
rbind.manova_lmm <- function(...){
    stop("Cannot use rbind on the output of rbind.anova_lmm. \n")
}

##----------------------------------------------------------------------
### rbind.anova_lmm.R ends here
