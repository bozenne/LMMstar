### rbind.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:51) 
## Version: 
## Last-Updated: jul 17 2024 (09:41) 
##           By: Brice Ozenne
##     Update #: 770
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
##' @param name.short [logical] use short names for the output coefficients, e.g., omit the regression variable name when the same regression variable is used in all models.
##' @param sep [character] character used to separate the name/outcome and the covariate when identifying the linear hypotheses.
##'
##' @details In presence of measurements from the same cluster across several models,
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
rbind.Wald_lmm <- function(model, ..., effects = NULL, rhs = NULL,
                           name = NULL, name.short = TRUE, sep = ": "){

    default <- LMMstar.options()
    call <- match.call()

    ## ** Check user input
    ## *** dots
    dots <- list(...)

    ## special case of single lmm
    if(length(dots)==0){
        if(NROW(model$multivariate)!=1){
            stop("Cannot handle multiple Multivariate Wald test in a single object. \n")
        }
        if(is.null(model$iid) || is.null(model$vcov)){ ## stop if effects is mean, variance, covariance, all
            stop("No iid or variance covariance matrix was stored in the Wald_lmm object. \n")
        }
        ls.contrast <- list(model$glht[[1]][[1]]$linfct)
        M.contrast <- diag(1, nrow = NROW(model$glht[[1]][[1]]$linfct))
        dimnames(M.contrast) <- list(rownames(M.contrast),rownames(M.contrast))

        if(is.null(name)){
            model$univariate <- cbind(model = 1, parameter = rownames(model$univariate), model$univariate)
        }else{
            model$univariate <- cbind(model = name, parameter = rownames(model$univariate), model$univariate)
        }
        model$glht[[1]][[1]]$linfct.original <- ls.contrast
        model$glht[[1]][[1]]$linfct <- M.contrast

        return(model) ## nothing to combine
    }else if(any(sapply(dots,inherits,"Wald_lmm")==FALSE)){
        stop("Extra arguments should inherit from Wald_lmm. \n")
    }else{
        ## combine input
        ls.object <- c(list(model),dots)
        n.object <- length(ls.object)
        ls.glht <- lapply(ls.object, function(iO){iO$glht[[1]][[1]]})
        all.coefnames <- names(unlist(lapply(ls.object,coef)))
        all.coefmodel <- unlist(lapply(1:n.object, function(iO){rep(iO, length(coef(ls.object[[iO]])))}))
        all.coeflabels <- unlist(lapply(ls.object, function(iO){rownames(iO$univariate)}))
        table.args <- cbind(do.call(rbind,lapply(ls.object,"[[","args")),
                            do.call(rbind,lapply(ls.object, function(iO){data.frame(n.univariate = NROW(iO$univariate),
                                                                                    iO$object[c("outcome","method.fit","type.information")])
                            })))
        alternative <- unique(sapply(ls.glht,"[[","alternative"))
    }

    ## *** content and available information in model and dots
    if(any(sapply(ls.object, function(iO){NROW(iO$multivariate)})!=1)){
        stop("Cannot handle multiple Multivariate Wald test in a single object. \n")
    }
    if(any(sapply(ls.object, function(iO){is.null(iO$iid) || is.null(iO$vcov)}))){ ## stop if effects is mean, variance, covariance, all
        stop("No iid or variance covariance matrix was stored in the Wald_lmm objects. \n")
    }
    if(any(sapply(ls.object, function(iO){iO$args$ci})==FALSE)){
        stop("All argument should contain a \"glht\" object, i.e. call anova with argument ci=TRUE. \n")
    }

    ## *** compatibility between model and dots
    test.compatibility <- c("ci","robust","df","transform.sigma","transform.k","transform.rho","backtransform","method.fit","type.information")
    if(NROW(unique(table.args[test.compatibility]))>1){
        pb <- names(which(lengths(apply(table.args[test.compatibility], MARGIN = 2, unique, simplify = FALSE))>1))
        stop("Element(s) \"",paste(pb, collapse = "\", \""),"\" should take the same value for all objects. \n")
    }

    if(length(alternative)>1){
        stop("Element \'alternative\' should take the same value for all glht objects. \n")
    }

    if(length(unique(lapply(ls.object, function(iO){iO$object$cluster.var})))>1){
        stop("Cluster variable differs between objects. \n")
    }
    cluster.var <- ls.object[[1]]$object$cluster.var

    ## *** name: find unique names for the models parameters
    if(!is.null(name) && length(name) == n.object){
        name.modelparam <- paste(name[all.coefmodel], all.coefnames, sep = sep[1])
        names(all.coefmodel) <- name[all.coefmodel]

        if(any(duplicated(name.modelparam))){
            stop("Duplicated names for the linear hypotheses between Wald_lmm objects. \n",
                 "Consider naming them when calling anova, e.g. anova(object, effect = c(\"myname1\"=\"X1=0\",\"myname2\"=\"X2=0\")). \n")
        }

        if(name.short && length(unique(all.coefnames))==1){
            attr(name.modelparam,"short") <- name
        }else{
            attr(name.modelparam,"short") <- name.modelparam
        }

    }else if(all(sapply(ls.object, function(iO){is.character(iO$multivariate$test)}))){ ## all hypotheses have been named in anova
        name.modelparam <- all.coeflabels

        if(any(duplicated(name.modelparam))){
            stop("Duplicated names for the linear hypotheses between Wald_lmm objects. \n",
                 "Consider using distinct names when calling anova, e.g. anova(object, effect = c(\"myname1\"=\"X1=0\",\"myname2\"=\"X2=0\")). \n")
        }

        attr(name.modelparam,"short") <- name.modelparam

    }else if(is.null(name) || identical(name,TRUE)){
        name.modelparam <- paste(table.args$outcome[all.coefmodel], all.coefmodel ,sep = sep)

        if(any(duplicated(name.modelparam))){
            stop("Duplicated names for the linear hypotheses between Wald_lmm objects. \n",
                 "Consider providing a distinct name for each object via the argument \'name\' \n",
                 "or naming the hypotheses when calling anova, e.g. anova(object, effect = c(\"myname1\"=\"X1=0\",\"myname2\"=\"X2=0\")). \n")
        }

        if(name.short && length(unique(all.coefnames))==1){  ## should be of correct size since when several coef per model the if condition would not be met
            attr(name.modelparam,"short") <- table.args$outcome
        }else if(name.short && length(unique(table.args$outcome))==1){
            attr(name.modelparam,"short") <- all.coefnames
        }else{
            attr(name.modelparam,"short") <- name.modelparam
        }

    }else if(all(is.na(name)) || identical(name,FALSE)){
        name.modelparam <- all.coefnames

        if(any(duplicated(name.modelparam))){
            stop("Duplicated names for the linear hypotheses between Wald_lmm objects. \n",
                 "Consider providing a distinct name for each object via the argument \'name\' \n",
                 "or naming the hypotheses when calling anova, e.g. anova(object, effect = c(\"myname1\"=\"X1=0\",\"myname2\"=\"X2=0\")). \n")
        }

        attr(name.modelparam,"short") <- name.modelparam

    }else if(length(name)!=n.object){
        stop("Argument \'name\' should have length ",n.object,", i.e. as many elements as Wald_lmm objects. \n")
    }
    n.modelparam <- length(name.modelparam)

    ## *** effects
    if(!is.null(effects)){

        if(is.matrix(effects)){
            contrast <- effects
            if(NCOL(contrast)!=n.modelparam){
                stop("Incorrect contrast matrix: should have ",n.modelparam," columns.\n",
                     "(one for each univariate test) \n")
            }
            if(is.null(rownames(contrast))){
                rownames(contrast) <- unname(apply(contrast, MARGIN = 1, function(iC){paste(names(iC)[iC!=0], collapse = ", ")}))
                if(any(duplicated(rownames(contrast)))){
                    stop("Missing rownames in argument \'effects\'. \n")
                }
            }else if(any(duplicated(rownames(contrast)))){
                stop("Duplicated rownames in argument \'effects\'.\n",
                     "They should be unique as they will be used to name the corresponding estimates. \n")
            }
            if(is.null(colnames(contrast))){
                colnames(contrast) <- name.modelparam
            }else if(!identical(sort(colnames(contrast)), sort(name.modelparam))){
                stop("Incorrect column names for argument \'effects\'.\n",
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
            contrast <- multcomp::contrMat(stats::setNames(rep(1,length(name.modelparam)),attr(name.modelparam,"short")), type = effects)
            colnames(contrast) <- name.modelparam
        }
    }else{
        contrast <- diag(1, nrow = n.modelparam, ncol = n.modelparam)
        dimnames(contrast) <- list(attr(name.modelparam,"short"), name.modelparam)
    }
    n.test <- NROW(contrast)

    ## *** rhs
    if(!is.null(rhs)){
        if(length(rhs)!=n.test){
            stop("Incorrect rhs: should have ",n.test," values.\n",
                 "(one for each univariate test) \n")
        }
    }else{
        if(is.null(effects)){
            rhs <- unlist(lapply(ls.object,function(iO){iO$univariate$null}))
        }else{
            rhs <- rep(0, NROW(contrast))
        }
        names(rhs) <- rownames(contrast)
    }

    
    ## ** merge information from all objects
    newtable.args <- cbind(type = ifelse(length(unique(unlist(table.args$type))>1), "all", table.args$type[[1]]),
                           sep = sep,
                           table.args[1,setdiff(test.compatibility, c("method.fit","type.information"))])

    newobject <- list(outcome = table.args$outcome,
                      cluster.var = cluster.var,
                      method.fit = table.args$method.fit[1],
                      type.information = table.args$type.information[1])
    if(!is.null(name) && length(name) == n.object){
        names(newobject$outcome) <- name
    }

    ## all parameter names (newname, oldname)
    if(!is.null(name) && length(name)==n.object){
        ls.theta.name <- lapply(1:n.object, function(iO){
            stats::setNames(paste(name[iO], sep[1], colnames(ls.object[[iO]]$glht[[1]][[1]]$linfct),sep=""), colnames(ls.object[[iO]]$glht[[1]][[1]]$linfct))
        })
    }else{
        ls.theta.name <- lapply(1:n.object, function(iO){
            stats::setNames(paste(iO, sep[1], colnames(ls.object[[iO]]$glht[[1]][[1]]$linfct),sep=""),colnames(ls.object[[iO]]$glht[[1]][[1]]$linfct))
        })
    }
    theta.name <- unlist(ls.theta.name)

    ## all contrast
    theta.contrast <- contrast %*% as.matrix(Matrix::bdiag(lapply(ls.object, function(iO){iO$glht[[1]][[1]]$linfct})))
    colnames(theta.contrast)  <- theta.name

    ## ** Extract elements from anova object
    ## *** univariate Wald test
    table.univariate <- do.call(rbind,lapply(ls.object, "[[", "univariate"))
    ## combine back-transformed estimates when a combination e.g. k\sigma^2\rho as each element should be transformed separately
    attr(table.univariate, "backtransform") <- unlist(lapply(ls.object, function(iO){attr(iO$univariate,"backtransform")}))

    ## *** cluster
    ls.cluster <- lapply(ls.object, function(iO){iO$object[["cluster"]]}) ## prefer [[ to $ to avoid partial matching (i.e. not output cluster.var if cluster is missing)
    seq.cluster <- unique(unlist(ls.cluster))
    n.cluster <- length(seq.cluster)
    newobject$cluster <- seq.cluster

    if(any(duplicated(unlist(ls.cluster)))){
        if(cluster.var=="XXclusterXX"){
            stop("Unable to decide whether observations from different models are matched or independent. \n",
                 "Consider specifying the \'repetition\' argument fitting the linear mixed model. \n")
        }
        independence <- FALSE
    }else{
        independence <- TRUE
    }
    if(all(sapply(ls.object, function(iO){attr(iO$args$robust,"call")})==FALSE)){
        newtable.args$robust <- independence==FALSE
    }
    
    ## *** estimate
    beta.estimate <- stats::setNames(table.univariate$estimate, name.modelparam)

    if(independence){
        
        beta.vcov <- as.matrix(do.call(Matrix::bdiag,lapply(ls.object, "[[", "vcov")))
        
    }else{
        if(any(sapply(ls.object, function(iO){is.null(iO$iid)}))){
            stop("Cannot evaluate the joint influence function of the variance-covariance parameters. \n",
                 "Consider using method.fit=\"ML\" when calling lmm. \n")
        }
        M.iid <- matrix(0, nrow = n.cluster, ncol = length(beta.estimate),
                        dimnames = list(seq.cluster, names(beta.estimate)))
        for(iO in 1:n.object){ ## iO <- 1            
            M.iid[ls.cluster[[iO]],which(all.coefmodel == iO)] <- ls.object[[iO]]$iid[,all.coefnames[which(all.coefmodel == iO)],drop=FALSE]
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
    }
    dimnames(beta.vcov) <- list(name.modelparam,name.modelparam)

    ## ** Combine elements

    ## *** multivariate tests
    outSimp <- simplifyContrast(contrast, rhs = rhs) ## remove extra lines
    C.vcov.C <- outSimp$C %*% beta.vcov %*% t(outSimp$C)
    C.vcov.C_M1 <- try(solve(C.vcov.C), silent = TRUE)

    if(inherits(C.vcov.C_M1,"try-error")){
        multistat <- NA
        attr(multistat,"error") <- "\n  Could not invert the covariance matrix for the proposed contrast."
        warning(attr(multistat,"error"))
    }else{
        multistat <- as.double(t(outSimp$C %*% beta.estimate - outSimp$rhs) %*% C.vcov.C_M1 %*% (outSimp$C %*% beta.estimate - outSimp$rhs))/outSimp$dim 
    }
    vec.null <- paste0(rownames(outSimp$C),"=",outSimp$rhs)

    Utype <- unique(unlist(lapply(ls.object, function(iO){iO$multivariate$type})))
    if(length(Utype)>1){
        Utype <- "all"
    }
    newtable.multivariate <- data.frame(type = Utype,
                                        test = 1,
                                        null = paste(vec.null, collapse = ", "),
                                        statistic = multistat,
                                        df.num = outSimp$dim,
                                        df.denom = Inf,
                                        p.value = 1 - stats::pf(multistat, df1 = n.test, df2 = Inf))

    ## *** univariate
    C.vcov.C <- contrast %*% beta.vcov %*% t(contrast)
    newtable.univariate <- data.frame(model = NA,
                                      type = Utype,
                                      estimate = (contrast %*% beta.estimate)[,1],
                                      se = sqrt(diag(C.vcov.C)),
                                      df = Inf,
                                      statistic = NA,
                                      lower = NA,
                                      upper = NA,
                                      null = rhs,
                                      p.value = NA)
    rownames(newtable.univariate) <- rownames(contrast)
    newtable.univariate$statistic <- newtable.univariate$estimate/newtable.univariate$se

    ## **** degree of freedom calculation
    if(all(contrast[lower.tri(contrast)+upper.tri(contrast)>0]==0)){
        newtable.univariate$df <- table.univariate$df
    }else if(diff(range(table.univariate$df))<0.1){
        newtable.univariate$df <- round(mean(table.univariate$df))
    }else if(all(lengths(lapply(ls.object,"[[","dVcov"))!=0)){
        ## when dVcov is available and higher level contrasts only select parameters
        ## ADD-HOC APPROXIMATION (ignores correlation in dVcov across models)
        theta.iid <- matrix(0, nrow = n.cluster, ncol = length(theta.name), dimnames = list(seq.cluster, theta.name))
        theta.dVcov <- array(0, dim = rep(length(theta.name),3), dimnames = list(theta.name,theta.name,theta.name))
        for(iO in 1:n.object){ ## iO <- 1
            theta.dVcov[ls.theta.name[[iO]],ls.theta.name[[iO]],ls.theta.name[[iO]]] <- ls.object[[iO]]$dVcov
            theta.iid[rownames(attr(ls.object[[iO]]$dVcov,"iid")),ls.theta.name[[iO]]] <- attr(ls.object[[iO]]$dVcov,"iid")
        }
        newtable.univariate$df <- .dfX(X.beta = theta.contrast, vcov.param = crossprod(theta.iid), dVcov.param = theta.dVcov)
        if(!independence){
            attr(newtable.univariate, "message.df") <- "neglecting correlation between parameters from different models in dVcov"
        }
    }
    newtable.args$df <- any(is.infinite(newtable.univariate$df)==FALSE)

    ## **** add columns
    if(!is.null(name) && length(name) == n.object){
        newtable.univariate$model <- lapply(1:NROW(contrast), function(iRow){names(all.coefmodel)[contrast[iRow,]!=0]})
    }else{
        newtable.univariate$model <- lapply(1:NROW(contrast), function(iRow){all.coefmodel[contrast[iRow,]!=0]})
    }
    newtable.univariate$test <- "1"
    rownames(newtable.univariate) <- rownames(contrast)
    newtable.univariate$parameter <- lapply(contrast2name(theta.contrast, ignore.value = TRUE), function(iE){unique(names(theta.name)[theta.name %in% iE$name])})

    ## **** take care of parameters subject to transformation for statistical inference
    if(!is.null(attr(table.univariate,"backtransform"))){ ## previously already a linear combination 
        param.untransformed <- stats::setNames(attr(table.univariate,"backtransform"),name.modelparam)
        attr(newtable.univariate,"backtransform") <- (contrast %*% param.untransformed[colnames(contrast)])[,1]
    }else if(newtable.args$backtransform && (any(contrast %in% 0:1 == FALSE) || any(rowSums(contrast != 0)>1))){ ## only now a linear combinations (backtransform no more ok)
        index.nobacktransform <- apply(contrast, MARGIN = 1, function(iRow){sum(iRow %in% 0:1 == FALSE) + sum(iRow != 0)>1})
        param.untransformed <- stats::setNames(unlist(lapply(ls.object, coef.Wald_lmm, backtransform = TRUE)), name.modelparam)
        attr(newtable.univariate,"backtransform") <- (contrast %*% param.untransformed[colnames(contrast)])[index.nobacktransform,1]
    }


    ## *** glht
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

    ## ** export
    e.glht$linfct.original <- lapply(ls.glht, "[[", "linfct")
    out <- list(multivariate = newtable.multivariate[,names(model$multivariate),drop=FALSE],
                univariate = newtable.univariate[,c("model","parameter",names(model$univariate)),drop=FALSE],
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
