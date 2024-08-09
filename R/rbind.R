### rbind.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:51) 
## Version: 
## Last-Updated: aug  8 2024 (17:17) 
##           By: Brice Ozenne
##     Update #: 1150
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * rbind.Wald_lmm (documentation)
##' @title Combine Wald Tests From Linear Mixed Models
##' @description Combine linear hypothesis tests from possibly different linear mixed models.
##'
##' @param model a \code{Wald_lmm} object (output of \code{anova} applied to a \code{lmm} object)
##' @param ...  possibly other \code{Wald_lmm} objects
##' @param effects [character or numeric matrix] how to combine the left-hand side of the hypotheses.
##' By default identity matrix but can also be \code{"Dunnett"},  \code{"Tukey"}, or  \code{"Sequen"} (see function \code{multcomp::contrMat} from the multcomp package).
##' @param rhs [numeric vector] the right hand side of the hypothesis. Should have the same length as the number of row of argument \code{effects}.
##' @param univariate [logical] Should an estimate, standard error, confidence interval, and p-value be output for each hypothesis?
##' @param multivariate [logical] Should all hypotheses be simultaneously tested using a multivariate Wald test?
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
##' 
##' ## combine null hypotheses
##' ## - model-based standard errors
##' AAA <- anova(e.lmm1, effect = c("X1|X2,X3"="X1=0","X2|X1,X3"="X2=0"))
##' BBB <- anova(e.lmm2, effect = c("X1|X8,X9"="X1=0"))
##' ZZZ <- rbind(AAA,BBB)
##' summary(ZZZ) ## adjusted for multiple testing
##' rbind(model.tables(e.lmm1)[2:3,], model.tables(e.lmm2)[2,,drop=FALSE])
##'
##' ## select null hypotheses & combine (model-based like standard errors)
##' AA <- anova(e.lmm1, effect = c("X1|X2,X3"="X1=0","X2|X1,X3"="X2=0"),
##'              robust = TRUE)
##' BB <- anova(e.lmm2, effect = c("X1|X8,X9"="X1=0"),
##'              robust = TRUE)
##' ZZ <- rbind(AA,BB)
##' summary(ZZ)  ## adjusted for multiple testing
##' rbind(model.tables(e.lmm1, robust = TRUE)[2:3,],
##'       model.tables(e.lmm2, robust = TRUE)[2,,drop=FALSE])

## * rbind.Wald_lmm (code)
##' @export
rbind.Wald_lmm <- function(model, ..., effects = NULL, rhs = NULL,
                           univariate = TRUE, multivariate = TRUE, 
                           name = NULL, name.short = TRUE, sep = ": "){

    mycall <- match.call()

    ## ** Check user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    
    ## special case of single lmm
    if(length(dots)==0){
        if(model$args$univariate==FALSE){
            stop("Cannot combine univariate Wald tests when they have not been stored. \n",
                 "Consider setting the argument \'univariate\' to TRUE when calling anova. \n")
        }
        if(length(unique(model$univariate$name))>1){
            stop("Cannot handle multiple multivariate Wald tests in a single object. \n")
        }

        ## add model to object
        if(!is.null(name)){
            if(length(name) != 1){
                stop("Argument \'name\' should have length 1, i.e., as many elements as Wald_lmm objects. \n")
            }
            if(any(name=="all")){
                stop("Argument \'name\' should not contain the value \"all\". \n")
            }
            model$object$model <- name        
            rownames(model$glht[[1]]$linfct) <- paste(name, rownames(model$glht[[1]]$linfct), sep = sep)
        }else{
            model$object$model <- model$object$outcome
        }
        model$object$sep <- sep
        model$object$param$model <- model$object$model

        class(model) <- append(c("rbindWald_lmm","Wald_lmm"),class(model))
        return(model) 
    }else if(any(sapply(dots,inherits,"Wald_lmm")==FALSE)){
        stop("Extra arguments should inherit from Wald_lmm. \n")
    }else{
        ## combine input
        ls.object <- c(list(model),dots)
        n.object <- length(ls.object)

        table.args <- cbind(do.call(rbind,lapply(ls.object,"[[","args")),
                            do.call(rbind,lapply(ls.object, function(iO){data.frame(n.univariate = NROW(iO$univariate),
                                                                                    iO$object[c("outcome","method.fit","type.information")])
                            })))
        alternative <- unique(sapply(ls.object, function(iO){iO$glht[[1]]$alternative}))

        ls.contrast <- lapply(ls.object, function(iO){stats::model.tables(iO, effects = "contrast", transform.names = FALSE, simplify = FALSE)[[1]]})

        all.table.param <- do.call(rbind,lapply(1:n.object, function(iO){cbind(model = iO, stats::model.tables(ls.object[[iO]], effects = "param"))}))
        rownames(all.table.param) <- NULL
        Wald.table.param <- do.call(rbind,lapply(1:n.object, function(iO){
            iCoef <- stats::coef(ls.object[[iO]], effects = "Wald", options = options)
            return(data.frame(model = iO, name = names(iCoef), value = iCoef))
        }))
    }

    ## *** content and available information in model and dots
    if(any(sapply(ls.object, function(iO){iO$args$univariate})==FALSE)){
        stop("Cannot combine univariate Wald tests when they have not been stored. \n",
             "Consider setting the argument \'univariate\' to TRUE when calling anova. \n")
    }
    if(any(sapply(ls.object, function(iO){length(unique(iO$univariate$name))})!=1)){
        stop("Cannot handle multiple multivariate Wald tests in a single object. \n")
    }

    ## *** compatibility between model and dots
    test.compatibility <- c("robust","df","transform.sigma","transform.k","transform.rho","transform.all","method.fit","type.information", "univariate", "multivariate","simplify")
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

    ## *** name (object)
    if(is.null(name)){
        duplicated.coefnames <- duplicated(paste(table.args$outcome[all.table.param$model], all.table.param$trans.name, sep = sep))
        duplicated.hypo <- duplicated(paste(table.args$outcome, Wald.table.param$name, sep = sep))
        if(any(duplicated.hypo) || any(duplicated.coefnames)){
            name <- 1:n.object
        }else{
            name <- table.args$outcome     
        }
    }else{
        if(length(name) != n.object){
            stop("Argument \'name\' should have length ",n.object,", i.e., as many elements as Wald_lmm objects. \n")
        }else if(any(duplicated(paste(name[all.table.param$model], all.table.param$trans.name, sep = sep))) || any(duplicated(paste(name[Wald.table.param$model], Wald.table.param$name, sep = sep)))){
            stop("Argument \'name\' should not contain duplicated values. \n")
        }else if(any(name=="all")){
            stop("Argument \'name\' should not contain the value \"all\". \n")
        }
    }
    all.coefUnames <- paste(name[all.table.param$model], all.table.param$trans.name, sep = sep) ## no more duplicates if same param in different models (e.g. (Intercept), age, or sigma)
    all.coefUnamesO <- paste(name[all.table.param$model], all.table.param$name, sep = sep) ## same but without transformation in the name
    
    if(name.short && all(duplicated(Wald.table.param$name)==FALSE)){  
        hypo.name <- Wald.table.param$name
    }else if(name.short && all(duplicated(name[Wald.table.param$model])==FALSE)){
        hypo.name <- name[Wald.table.param$model]
    }else{
        hypo.name <- paste(name[Wald.table.param$model], Wald.table.param$name, sep = sep)
    }

    ## *** effects
    if(!is.null(effects)){

        if(is.matrix(effects)){
            contrast <- effects

            ## number of columns
            if(NCOL(contrast)!=NROW(Wald.table.param)){
                stop("Incorrect contrast matrix: should have ",NROW(Wald.table.param)," columns.\n",
                     "(one for each univariate test) \n")
            }

            ## column ordering
            if(!is.null(colnames(contrast))){
                if(all(colnames(contrast) %in% Wald.table.param$name)){
                    if(all(duplicated(Wald.table.param$name)==FALSE)){
                        contrast <- contrast[,Wald.table.param$name,drop=FALSE]
                    }else{
                        stop("Ambiguous names for argument \'effect\' due to duplicated model parameter names. \n",
                             "Consider providing a distinct name for each object via the argument \'name\' \n",
                             "or naming the hypotheses when calling anova, e.g. anova(object, effect = c(\"myname1\"=\"X1=0\",\"myname2\"=\"X2=0\")). \n")
                    }
                }else if(all(colnames(contrast) %in% hypo.name)){
                    contrast <- contrast[,hypo.name,drop=FALSE]
                }else if(all(colnames(contrast) %in% paste(name, Wald.table.param$name, sep = sep))){
                    contrast <- contrast[,paste(name[Wald.table.param$model], Wald.table.param$name, sep = sep),drop=FALSE]
                }else{
                    if(is.null(call$name)){
                        stop("Incorrect column names for argument \'effects\'.\n",
                             "Should match \"",paste(paste(name[Wald.table.param$model], Wald.table.param$name, sep = sep), collapse="\" \""),"\".\n")
                    }else{
                        stop("Incorrect column names for argument \'effects\'.\n",
                             "Should match \"",paste(Wald.table.param$name, collapse="\" \""),"\".\n")
                    }
                }
            }
            colnames(contrast) <- NULL

            ##  rows
            if(is.null(rownames(contrast))){
                rownames(contrast) <- unname(apply(contrast, MARGIN = 1, function(iC){paste(hypo.name[iC!=0], collapse = ", ")}))
                if(any(duplicated(rownames(contrast)))){
                    stop("Missing rownames in argument \'effects\'. \n")
                }
            }else if(any(duplicated(rownames(contrast)))){
                stop("Duplicated rownames in argument \'effects\'.\n",
                     "They should be unique as they will be used to name the corresponding estimates. \n")
            }
            
        }else if(all(is.character(effects))){
            if(length(effects)>1){
                stop("When a character, argument \'effects\' should have length 1. \n")
            }
            valid.contrast <- c("Dunnett","Tukey","Sequen")
            if(effects %in% valid.contrast == FALSE){
                stop("When a character, argument \'effects\' should be one of \"",paste(valid.contrast, collapse = "\" \""),"\". \n")
            }
            contrast <- multcomp::contrMat(stats::setNames(rep(1,NROW(Wald.table.param)),hypo.name), type = effects)
            colnames(contrast) <- NULL
        }
    }else{
        contrast <- diag(1, nrow = NROW(Wald.table.param), ncol = NROW(Wald.table.param))
        dimnames(contrast) <- list(hypo.name, NULL)
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

    ## ** Extract elements from anova object

    ## *** cluster
    ls.cluster <- lapply(ls.object, function(iO){iO$object[["cluster"]]}) ## prefer [[ to $ to avoid partial matching (i.e. not output cluster.var if cluster is missing)
    seq.cluster <- unique(unlist(ls.cluster))
    n.cluster <- length(seq.cluster)
    
    if(any(duplicated(unlist(ls.cluster)))){
        if(cluster.var=="XXclusterXX"){
            stop("Unable to decide whether observations from different models are matched or independent. \n",
                 "Consider specifying the \'repetition\' argument fitting the linear mixed model. \n")
        }
        independence <- FALSE
        if(any(sapply(ls.object, lava::iid, effects = "test") == FALSE)){ ## stop if effects is mean, variance, covariance, all
            stop("Cannot combine tests from overlapping population without the influence function. \n",
                 "Consider setting argument \'simplify\' to FALSE when calling anova. \n")
        }
    }else{
        independence <- TRUE
    }

    ## *** estimate
    all.coefvalues <- stats::setNames(all.table.param$trans.value, all.coefUnames)
    all.coefvaluesO <- stats::setNames(all.table.param$value, all.coefUnamesO)

    ## *** contrast
    all.contrast <- contrast %*% as.matrix(Matrix::bdiag(ls.contrast))
    colnames(all.contrast)  <- all.coefUnamesO

    ## *** vcov
    if(independence){

        ls.vcov <- lapply(ls.object, stats::vcov, effects = "all", options = options)
        all.vcov <- as.matrix(do.call(Matrix::bdiag,ls.vcov))
        dimnames(all.vcov) <- list(all.coefUnamesO, all.coefUnamesO)
        all.iid <- NULL

    }else{

        all.iid <- matrix(0, nrow = n.cluster, ncol = NROW(all.table.param),
                          dimnames = list(seq.cluster, all.coefUnamesO))

        for(iO in 1:n.object){ ## iO <- 1
            ## iid
            iAll.iid <- lava::iid(ls.object[[iO]], effects = "all", options = options)
            iAll.coefindex <- match(paste(name[iO], colnames(iAll.iid), sep = sep), all.coefUnames) ## handle the transformation of the names
            all.iid[rownames(iAll.iid),iAll.coefindex] <- iAll.iid
            attr(all.iid,"message") <- union(attr(all.iid,"message"), attr(iAll.iid,"message"))
        }

        if(any(is.na(all.iid))){
            all.vcov <- tcrossprod(sqrt(apply(all.iid^2, 2, sum, na.rm = TRUE))) * stats::cor(all.iid, use = "pairwise")
            ## usually better compared to formula 11.43 from chapter 11.4 of the book High-dimensional statistics by WAINWRIGHT
            ## iIDD0 <- M.iid/(1-mean(is.na(M.iid)))
            ## iIDD0[is.na(M.iid)] <- 0
            ## out$vcov <- crossprod(iIDD0) - mean(is.na(M.iid))*diag(diag(crossprod(iIDD0)))
            ## out$vcov - crossprod(M.iid)
        }else{
            all.vcov <- crossprod(all.iid)
        }
        
    }

    ## *** dVcov
    if(table.args$df[1]){

        red.coefUnamesO <- unlist(lapply(1:n.object, function(iO){paste(name[iO],names(which(colSums(ls.contrast[[iO]]!=0)>0)), sep = sep)}))
        all.dVcov <- array(0, dim = c(rep(length(red.coefUnamesO),2), NROW(all.table.param)),
                           dimnames = list(red.coefUnamesO,red.coefUnamesO,all.coefUnamesO))
        
        for(iO in 1:n.object){ ## iO <- 1
            idVcov <- attr(stats::vcov(ls.object[[iO]], effects = c("all","gradient"), options = options),"gradient")
            iNewnames <- all.coefUnamesO[match(paste(name[iO], dimnames(idVcov)[[1]], sep = sep), all.coefUnames)]
            dimnames(idVcov) <- list(iNewnames, iNewnames, all.coefUnamesO[match(paste(name[iO], dimnames(idVcov)[[3]], sep = sep), all.coefUnames)])
            
            iDimnames <- lapply(1:3, function(iD){intersect(dimnames(all.dVcov)[[iD]],dimnames(idVcov)[[iD]])})           
            all.dVcov[iDimnames[[1]],iDimnames[[2]],iDimnames[[3]]] <- idVcov[iDimnames[[1]],iDimnames[[2]],iDimnames[[3]]]
        }

        ## add model-based vcov when robust=1, i.e., compute df from model-based vcov even though one uses robust vcov
        if(table.args$robust[1]==1){
            original.se <- sqrt(unlist(lapply(ls.object, function(iO){
                iGrad <- attr(stats::vcov(iO, effects = c("all","gradient"), options = options),"gradient")
                return(diag(attr(iGrad,"vcov")))
            })))
            attr(all.dVcov,"vcov") <- tcrossprod(original.se) * all.vcov / tcrossprod(sqrt(diag(all.vcov)))
        }

    }else{
        all.dVcov <- NULL
    }

    ## ** Combine elements
    out <- .anova_Wald(param = all.coefvalues,
                       param.notrans = all.coefvaluesO,
                       vcov.param = all.vcov,
                       dVcov.param = all.dVcov,
                       iid.param = NULL,
                       type.param = stats::setNames(all.table.param$type, all.coefUnamesO),
                       contrast = list(user = list(user = all.contrast)),
                       null = list(user = list(user = rhs)),
                       robust = table.args$robust[1],
                       df = table.args$df[1],
                       multivariate = multivariate,
                       univariate = univariate,
                       simplify = table.args$simplify[1],
                       transform.sigma = table.args$transform.sigma[1],
                       transform.k = table.args$transform.k[1],
                       transform.rho = table.args$transform.rho[1],
                       backtransform = any(do.call(rbind,lapply(ls.object,"[[","univariate"))$tobacktransform))

    ## update term in univariate
    out$univariate$term <- apply(contrast, MARGIN = 1, function(iRow){
        iOut <- which(iRow!=0)
        if(length(iOut)>1){
            return("user")
        }else{
            return(Wald.table.param$name[iOut])
        }
    })

    ## add model to object
    univariate.model <- apply(contrast, MARGIN = 1, function(iRow){
        iOut <- unique(Wald.table.param$model[iRow!=0])
        if(length(iOut)>1){
            return("all")
        }else{
            return(as.character(iOut))
        }
    })

    out$univariate <- cbind(model = univariate.model, out$univariate)

    ## add extra information about the original lmm
    out$object <- list(outcome =  table.args$outcome,
                       model = as.character(name),
                       method.fit = table.args$method.fit[1],
                       type.information = table.args$type.information[1],
                       cluster.var = cluster.var,
                       cluster = seq.cluster,
                       independence = independence,
                       linfct = ls.contrast,
                       param = cbind(model = as.character(all.table.param$model),
                                     all.table.param[setdiff(names(all.table.param),"model")],
                                     Uname = all.coefUnamesO, trans.Uname = all.coefUnames)
                       )

    ## check wheter parameters from different hypotheses are combined
    test.hypoCross <- apply(contrast, MARGIN = 1, function(iRow){length(unique(Wald.table.param$model[iRow!=0]))})>1
    if(any(test.hypoCross)){
        attr(out$univariate,"message.se") <- "and the influence function"
        if(!is.null(attr(all.iid,"message"))){
            attr(out$univariate,"message.se") <- paste0(attr(out$univariate,"message.se")," \n  (",attr(all.iid,"message"),")")
        }
    }else if(table.args$df[1]>0){
        if(!independence){
            if(!is.null(attr(all.iid,"message")) && all(test.hypoCross==FALSE)){
                attr(out$univariate,"df") <- paste0("Satterthwaite approximation of the degrees-of-freedom \n  (",attr(all.iid,"message")," and no correlation between models in dVcov)")
            }else{
                attr(out$univariate,"df") <- "Satterthwaite approximation of the degrees-of-freedom \n  (neglecting parameter correlation between models in dVcov)"
            }
        }else if(!is.null(attr(all.iid,"message")) && all(test.hypoCross==FALSE)){
            attr(out$univariate,"df") <- paste0("Satterthwaite approximation of the degrees-of-freedom \n  (",attr(all.iid,"message"),")")
        }
    }

    ## ** export
    attr(out,"call") <- list(anova = lapply(ls.object,attr,"call"),
                             rbind = mycall)
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
