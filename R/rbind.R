### rbind.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:51) 
## Version: 
## Last-Updated: May 12 2024 (17:06) 
##           By: Brice Ozenne
##     Update #: 505
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * rbind.Wald_lmm (documentation)
##' @title Linear Hypothesis Testing Across Linear Mixed Models
##' @description Linear hypothesis testing accross linear mixed model.
##'
##' @param model a \code{Wald_lmm} object (output of \code{anova} applied to a \code{lmm} object)
##' @param ...  possibly other \code{Wald_lmm} objects
##' @param effects [character or numeric matrix] how to combine the left-hand side of the hypotheses.
##' By default identity matrix but can also be \code{"Dunnett"},  \code{"Tukey"}, or  \code{"Sequen"} (see function \code{multcomp::contrMat} from the multcomp package).
##' @param rhs [numeric vector] the right hand side of the hypothesis. Should have the same length as the number of row of argument \code{effects}.
##' @param name [character vector or NULL] character used to identify each model in the output.
##' By default, use the name of the outcome of the model.
##' @param sep [character] character used to separate the outcome and the covariate when naming the tests.
##'
##' @details WARNING: in presence of measurements from the same cluster across several models,
##' the influence function is used to estimate the covariance between the model parameters.
##' This is why the (robust) standard errors may not match the (model-based) standard error from the linear mixed
##' Setting the argument \code{robust} to \code{FALSE} when calling \code{anova.lmm} will rescale the (robust) standard errors to mimic the original model-based standard errors.
##' @keywords methods
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
##' model.tables(e.lmm1) ## model-based standard errors
##' model.tables(e.lmm1, robust = TRUE) ## robust standard errors
##' 
##' ## select null hypotheses & combine (robust standard errors)
##' AAA <- anova(e.lmm1, ci = TRUE, effect = c("X1|X2,X3"="X1=0","X2|X1,X3"="X2=0"))
##' BBB <- anova(e.lmm2, ci = TRUE, effect = c("X1|X8,X9"="X1=0"))
##'
##' ZZZ <- rbind(AAA,BBB)
##'
##' ## select null hypotheses & combine (model-based like standard errors)
##' AA <- anova(e.lmm1, ci = TRUE, effect = c("X1|X2,X3"="X1=0","X2|X1,X3"="X2=0"),
##'              robust = FALSE)
##' BB <- anova(e.lmm2, ci = TRUE, effect = c("X1|X8,X9"="X1=0"),
##'              robust = FALSE)
##' ZZ <- rbind(AA,BB)

## * rbind.Wald_lmm (code)
##' @export
rbind.Wald_lmm <- function(model, ..., effects = NULL, rhs = NULL, name = NULL, sep = ": "){

    default <- LMMstar.options()
    call <- match.call()
    
    ## ** Check user input
    dots <- list(...)
    if(length(dots)==0){
        return(model) ## nothing to combine
    }else if(any(sapply(dots,inherits,"Wald_lmm")==FALSE)){
        stop("Extra arguments should inherit from Wald_lmm. \n")
    }
    ls.object <- c(list(model),dots)
    n.object <- length(ls.object)
    if(!is.null(name) && length(name)!=n.object){
        stop("Argument \'name\' should have length ",n.object,", i.e. as many elements as anova objects. \n")
    }
    if(any(sapply(ls.object, function(iO){NROW(iO$multivariate)})!=1)){
        stop("Cannot handle multiple Multivariate Wald test in a single object. \n")
    }

    table.args <- cbind(do.call(rbind,lapply(ls.object,"[[","args")),
                        do.call(rbind,lapply(ls.object, function(iO){data.frame(n.test = NROW(iO$univariate),
                                                                                iO$object[c("outcome","method.fit","type.information")])
                        })))

    Utype <- unique(unlist(table.args$type))

    newtable.args <- data.frame(type = ifelse(length(Utype)>1,"all",Utype),
                                sep = sep, 
                                table.args[1,c("robust","df","ci","transform.sigma","transform.k","transform.rho","backtransform")])

    newobject <- list()
    if(length(unique(table.args$outcome))==1){
        newobject$outcome <- table.args$outcome[1]
    }
    if(length(unique(table.args$method.fit))==1){
        newobject$method.fit <- table.args$method.fit[1]
    }
    if(length(unique(table.args$type.information))==1){
        newobject$type.information <- table.args$type.information[1]
    }
    
    if(any(table.args$ci==FALSE)){
        stop("All argument should contain a \"glht\" object, i.e. call anova with argument ci=TRUE. \n")
    }
    if(any(newtable.args$ACO!=table.args$ACO[-1])){
        stop("Element \'ACO\' should take the same value for all objects. \n")
    }
    if(any(newtable.args$ATE!=table.args$ATE[-1])){
        stop("Element \'ATE\' should take the same value for all objects. \n")
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
    if(any(newtable.args$backtransform!=table.args$backtransform[-1])){
        stop("Element \'backtransform\' should take the same value for all objects. \n")
    }
    ls.glht <- unname(unlist(unlist(lapply(ls.object,"[[","glht"),recursive = FALSE),recursive = FALSE))
    alternative <- ls.glht[[1]]$alternative
    if(any(alternative != sapply(ls.glht,"[[","alternative")[-1])){
        stop("Element \'alternative\' should take the same value for all glht objects. \n")
    }

    if(is.null(name)){
        name.modelparam <- unlist(lapply(ls.object, function(iO){rownames(iO$univariate)}))
        if(any(duplicated(name.modelparam))){
            name.modelparam <- unlist(lapply(ls.object, function(iO){paste0(iO$object$outcome, sep, rownames(iO$univariate))}))

            if(is.na(sep) || any(duplicated(name.modelparam))){
                stop("Univariate test should have distinct names. \n",
                     "Consider naming them when calling anova, e.g. anova(object, effect = c(\"myname1\"=\"X1=0\",\"myname2\"=\"X2=0\"))")
            }else{
                for(iO in 1:n.object){
                    rownames(ls.object[[iO]]$univariate) <- paste0(ls.object[[iO]]$object$outcome, sep, rownames(ls.object[[iO]]$univariate))
                }
            }
        }
    }else{
        name.modelparam <- unlist(lapply(1:n.object, function(iO){paste0(name[iO],sep,rownames(ls.object[[iO]]$univariate))}))
        for(iO in 1:n.object){ ## iO <- 1
            rownames(ls.object[[iO]]$univariate) <- paste0(name[iO],sep,rownames(ls.object[[iO]]$univariate))
            if(!is.null(attr(ls.object[[iO]]$univariate, "backtransform"))){
                names(attr(ls.object[[iO]]$univariate, "backtransform")) <- paste0(name[iO],sep,attr(ls.object[[iO]]$univariate, "backtransform"))
            }
        }
        if(any(duplicated(name.modelparam))){
            stop("Univariate test should have distinct names. \n",
                 "Consider naming them when calling anova, e.g. anova(object, effect = c(\"myname1\"=\"X1=0\",\"myname2\"=\"X2=0\"))")
        }
    }

    n.modelparam <- length(name.modelparam)
    if(!is.null(effects)){

        if(is.matrix(effects)){
            contrast <- effects
            if(NCOL(contrast)!=n.modelparam){
                stop("Incorrect contrast matrix: should have ",n.modelparam," columns.\n",
                     "(one for each univariate test) \n")
            }
            if(is.null(colnames(contrast))){
                colnames(contrast) <- name.modelparam
            }else if(!identical(sort(colnames(contrast)), sort(name.modelparam))){
                stop("Incorrect column names for argument contrast.\n",
                     "Should match \"",paste(sort(name.modelparam), collapse="\" \""),"\".\n")
            }else{
                contrast <- contrast[,name.modelparam,drop=FALSE]
            }
        }else if(all(is.character(effects))){
            if(length(effects)>1){
                stop("When a character, argument \'effects\' should have length 1. \n")
            }
            valid.contrast <- c("Dunnett","Tukey","Sequen")
            if(effects %in% valid.contrast == FALSE){
                stop("When a character, argument \'effects\' should be one of \"",paste(valid.contrast, collapse = "\" \""),"\". \n")
            }

            contrast <- multcomp::contrMat(rep(1,length(name.modelparam)), type = effects)
            colnames(contrast) <- name.modelparam
            try(rownames(contrast) <- unlist(lapply(strsplit(split = "-",rownames(contrast),fixed=TRUE), function(iVec){
                paste(name.modelparam[as.numeric(trimws(iVec))], collapse = " - ")
            })), silent = TRUE)
        }
    }else{
        contrast <- diag(1, nrow = n.modelparam, ncol = n.modelparam)
        dimnames(contrast) <- list(name.modelparam, name.modelparam)
    }
    n.test <- NROW(contrast)
    name.test <- rownames(contrast)

    if(!is.null(rhs)){
        if(length(rhs)!=n.test){
            stop("Incorrect rhs: should have ",n.test," values.\n",
                 "(one for each univariate test) \n")
        }
    }else if(is.null(effects)){
        rhs <- stats::setNames(unlist(lapply(ls.object,function(iO){iO$univariate$null})), name.test)
    }else{
        rhs <- rep(0, NROW(effects))
    }

    cluster.var <- ls.object[[1]]$object$cluster.var
    if(any(cluster.var!=unlist(lapply(ls.object[-1], function(iO){iO$object$cluster.var})))){
        stop("Cluster variable differs among objects. \n")
    }
    newobject$cluster.var <- cluster.var
    
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
        if(all(duplicated(gridTest$outcome))){
            gridTest$outcome <- NULL
        }
        if(all(duplicated(gridTest$type))){
            gridTest$type <- NULL
        }
        if(all(duplicated(gridTest$test))){
            gridTest$test <- NULL
        }
    }
    col.nametest <- colnames(gridTest)
    name.test <- unique(nlme::collapse(gridTest[,col.nametest,drop=FALSE], as.factor = FALSE))

    ## ** Extract elements from anova object
    ## *** univariate Wald test
    ls.univariate <- lapply(1:n.object,function(iO){ ## iO <- 1
        iTable <- cbind(outcome = unname(outcome[iO]), ls.object[[iO]]$univariate)
        iTable$name.test <- factor(nlme::collapse(iTable[,col.nametest,drop=FALSE], as.factor = FALSE),levels = name.test)
        rownames(iTable) <- rownames(ls.object[[iO]]$univariate)
        attr(iTable,"backtransform") <- attr(ls.object[[iO]]$univariate, "backtransform")
        return(iTable)
    })
    table.univariate <- do.call(rbind,ls.univariate)
    attr(table.univariate, "backtransform") <- unlist(lapply(ls.univariate,attr,"backtransform"))

    ## *** cluster
    ls.cluster <- lapply(ls.object, function(iO){iO$object$cluster})
    seq.cluster <- unique(unlist(ls.cluster))
    n.cluster <- length(seq.cluster)
    newobject$cluster <- seq.cluster

    independence <- (cluster.var=="XXclusterXX" || all(duplicated(unlist(ls.cluster))==FALSE)) ## each individual only appear once
    if(all(sapply(ls.object, function(iO){attr(iO$args$robust,"call")})==FALSE)){
        newtable.args$robust <- independence==FALSE
    }

    ## *** estimate
    beta.estimate <- stats::setNames(table.univariate$estimate, name.modelparam)
    if(independence){
        
        beta.vcov <- as.matrix(do.call(Matrix::bdiag,lapply(ls.object, "[[", "vcov")))
        
        if(all(diag(contrast)==1) && all(contrast[upper.tri(contrast, diag = FALSE)] == 0) && all(contrast[lower.tri(contrast, diag = FALSE)] == 0)){
            beta.df <- table.univariate$df
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
        if(any(!is.infinite(table.univariate$df)) && any(!is.na(table.univariate$df)) && abs(diff(range(table.univariate$df, na.rm = TRUE)))<0.1){
            beta.df <- rep(mean(table.univariate$df), length(beta.estimate))
            newtable.args$df <- TRUE
        }else{
            beta.df <- rep(Inf, length(beta.estimate))
            newtable.args$df <- FALSE
        }
    }
    dimnames(beta.vcov) <- list(name.modelparam,name.modelparam)
    names(beta.df) <- name.modelparam

    ## ** Combine elements

    ## *** multivariate tests
    outSimp <- .simplifyContrast(contrast, rhs = rhs) ## remove extra lines
    C.vcov.C <- outSimp$C %*% beta.vcov %*% t(outSimp$C)
    C.vcov.C_M1 <- try(solve(C.vcov.C), silent = TRUE)

    if(inherits(C.vcov.C_M1,"try-error")){
        multistat <- NA
        attr(multistat,"error") <- "\n  Could not invert the covariance matrix for the proposed contrast."
        warning(attr(multistat,"error"))
    }else{
        multistat <- as.double(t(outSimp$C %*% beta.estimate - outSimp$rhs) %*% C.vcov.C_M1 %*% (outSimp$C %*% beta.estimate - outSimp$rhs))/outSimp$dim 
    }
    Utype2 <- unique(unlist(lapply(ls.object, function(iO){iO$multivariate$type})))
    vec.null <- paste0(rownames(outSimp$C),"=",outSimp$rhs)

    newtable.multivariate <- data.frame(type = ifelse(length(Utype2)>1,"all",Utype2),
                                        test = 1,
                                        null = paste(vec.null, collapse = ", "),
                                        statistic = multistat,
                                        df.num = outSimp$dim,
                                        df.denom = Inf,
                                        p.value = 1 - stats::pf(multistat, df1 = n.test, df2 = Inf))

    ## *** univariate
    if(!is.null(effects)){
        C.vcov.C <- contrast %*% beta.vcov %*% t(contrast)
        newtable.univariate <- data.frame(outcome = NA,
                                          type = ifelse(length(Utype2)>1,"all",Utype2),
                                          test = 1,
                                          estimate = (contrast %*% beta.estimate)[,1],
                                          se = sqrt(diag(C.vcov.C)),
                                          df = Inf,
                                          statistic = NA,
                                          lower = NA,
                                          upper = NA,
                                          null = rhs,
                                          p.value = NA)
        rownames(newtable.univariate) <- rownames(contrast)
        newtable.univariate$outcome <- lapply(1:NROW(contrast), function(iRow){table.univariate[colnames(contrast)[contrast[iRow,]!=0],"outcome"]})
        newtable.univariate$statistic <- newtable.univariate$estimate/newtable.univariate$se
        newtable.args$df <- FALSE

        ## take care of parameters subject to transformation for statistical inference
        if(!is.null(attr(table.univariate,"backtransform"))){ ## previously already a linear combination 
            param.untransformed <- stats::setNames(attr(table.univariate,"backtransform"),name.modelparam)
            attr(newtable.univariate,"backtransform") <- (contrast %*% param.untransformed[colnames(contrast)])[,1]

        }else if(newtable.args$backtransform && (any(contrast %in% 0:1 == FALSE) || any(rowSums(contrast != 0)>1))){ ## only now a linear combinations (backtransform no more ok)
            index.nobacktransform <- apply(contrast, MARGIN = 1, function(iRow){sum(iRow %in% 0:1 == FALSE) + sum(iRow != 0)>1})
            param.untransformed <- stats::setNames(unlist(lapply(ls.object, coef.Wald_lmm, backtransform = TRUE)), name.modelparam)
            attr(newtable.univariate,"backtransform") <- (contrast %*% param.untransformed[colnames(contrast)])[index.nobacktransform,1]
        }
    }else{
        C.vcov.C <- beta.vcov
        newtable.univariate <- table.univariate
        ## update se and statistic in case robust standard errors have been used
        newtable.univariate$se <- sqrt(diag(beta.vcov))
        newtable.univariate$statistic <- newtable.univariate$estimate/newtable.univariate$se
    }
    if(newtable.args$df){
        e.glht <- list(linfct = contrast, rhs = rhs,
                       coef = table.univariate$estimate, vcov = beta.vcov, df = floor(stats::median(newtable.univariate$df)), alternative = alternative)
        if(!is.na(default$df)){
            e.glht$df <- pmax(e.glht$df, default$df)
        }
    }else{
        e.glht <- list(linfct = contrast, rhs = rhs,
                       coef = table.univariate$estimate, vcov = beta.vcov, df = FALSE, alternative = alternative)
    }
    class(e.glht) <- "glht"

    ## ** retrieve parameter associated to each test
    nParam.hypo <- unlist(lapply(ls.glht,function(iO){rowSums(iO$linfct!=0)}))
    if(all(contrast %in% 0:1) && all(rowSums(contrast)==1) && all(nParam.hypo==1)){ ## single parameter per test
        index.param <- which(contrast==1, arr.ind = TRUE)[,"col"]
        newtable.univariate[,"parameter"] <- names(nParam.hypo[index.param])
    }else{
        ls.names <- do.call(c, lapply(ls.glht, FUN = function(iO){
            lapply(1:NROW(iO$linfct), FUN = function(iRow){colnames(iO$linfct)[iO$linfct[iRow,]!=0]})
        }))
        if(all(contrast %in% 0:1) && all(rowSums(contrast)==1)){  ## several parameters per test but all from the same model
            newtable.univariate$parameter <- ls.names
        }else if(length(unique(ls.names))==1){  
            if(length(ls.names[[1]])==1){ ## contrasting one parameter accross models
                newtable.univariate$parameter <- ls.names[[1]]
            }else{ ## contrasting one set of parameters accross models
                newtable.univariate$parameter <- ls.names[1]
            }
        }
    }

    ## ** export
    e.glht$linfct.original <- lapply(ls.glht, "[[", "linfct")
    out <- list(multivariate = newtable.multivariate[,names(model$multivariate),drop=FALSE],
                univariate = newtable.univariate[,c("outcome","parameter",names(model$univariate)),drop=FALSE],
                glht = list(all = list("1" = e.glht)),
                object = newobject,
                args = newtable.args,
                vcov = C.vcov.C)
    attr(out$univariate,"backtransform") <- attr(newtable.univariate,"backtransform")
    attr(out$object,"independence") <- independence
    attr(out,"call") <- list(anova = lapply(ls.object,attr,"call"),
                             rbind = match.call())
    class(out) <- append(c("rbindWald_lmm","Wald_lmm"),class(out))
    return(out)
}


## * rbind.rbindWald_lmm (code)
##' @export
rbind.rbindWald_lmm <- function(...){
    stop("Cannot use rbind on the output of rbind.Wald_lmm. \n")
}

##----------------------------------------------------------------------
### rbind.R ends here
