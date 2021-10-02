### anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:38) 
## Version: 
## Last-Updated: okt  1 2021 (17:22) 
##           By: Brice Ozenne
##     Update #: 581
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * anova.lmm (documentation)
##' @title Multivariate Wald Tests For Linear Mixed Model
##' @description Perform a Wald test testing simultaneously several null hypotheses corresponding to linear combinations of the model paramaters. 
##' @name anova
##' 
##' @param object a \code{lmm} object. Only relevant for the anova function.
##' @param x an \code{anova_lmm} object. Only relevant for print and confint functions.
##' @param effects [character] Should the Wald test be computed for all variables (\code{"all"}),
##' or only variables relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only variables relative to the variance structure (\code{"variance"}),
##' or only variables relative to the correlation structure (\code{"correlation"}).
##' Can also be use to specify linear combinations of coefficients, similarly to the \code{linfct} argument of the \code{multcomp::glht} function.
##' @param rhs [numeric vector] the right hand side of the hypothesis. Only used when the argument effects is a matrix.
##' @param ci [logical] Should a confidence interval be output for each hypothesis?
##' @param level [numeric, 0-1] nominal coverage of the confidence intervals.
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.
##' @param print.null [logical] should the null hypotheses be printed in the console?
##' @param df [logical] Should a F-distribution be used to model the distribution of the Wald statistic. Otherwise a chi-squared distribution is used.
##' @param method [character] type of adjustment for multiple comparisons: one of \code{"none"}, \code{"bonferroni"}, \code{"single-step"}.
##' Not relevant for the global test (F-test or Chi-square test) - only relevant when testing each hypothesis and adjusting for multiplicity.
##' @param transform.sigma,transform.k,transform.rho,transform.names are passed to the \code{vcov} method. See details section in \code{\link{coef.lmm}}.
##' @param parm Not used. For compatibility with the generic method.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A list of matrices containing the following columns:\itemize{
##' \item \code{null}: null hypothesis
##' \item \code{statistic}: value of the test statistic
##' \item \code{df.num}: degrees of freedom for the numerator (i.e. number of hypotheses)
##' \item \code{df.denom}: degrees of freedom for the denominator (i.e. Satterthwaite approximation)
##' \item \code{p.value}: p-value.
##' }
##' as well as an attribute contrast containing the contrast matrix encoding the linear combinations of coefficients (in columns) for each hypothesis (in rows).
##' 
##' @details By default confidence intervals and p-values are adjusted based on the distribution of the maximum-statistic.
##' This is refered to as a single-step Dunnett multiple testing procedures in table II of Dmitrienko et al. (2013) and is performed using the multcomp package with the option \code{test = adjusted("single-step")}.
##'
##' @references Dmitrienko, A. and D'Agostino, R., Sr (2013), Traditional multiplicity adjustment methods in clinical trials. Statist. Med., 32: 5172-5218. https://doi.org/10.1002/sim.5990.
##'  
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit Linear Mixed Model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL)
##' 
##' ## chi-2 test
##' anova(eUN.lmm, df = FALSE)
##' 
##' ## F-test
##' anova(eUN.lmm)
##' anova(eUN.lmm, effects = "all")
##' anova(eUN.lmm, effects = c("X1=0","X2+X5=10"), ci = TRUE)
##' 

## * anova.lmm (code)
##' @rdname anova
##' @export
anova.lmm <- function(object, effects = NULL, rhs = NULL, df = !is.null(object$df), ci = FALSE, 
                      type.object = "lmm",
                      transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){
    
    
    ## ** normalized user input    
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(is.null(effects)){
        effects <- options$effects
    }

    if(inherits(effects,"lmm")){ ## likelihood ratio test
        out <- .anova_LRT(object1 = object, object2 = effects)
    }else{ ## Wald test
        out <- .anova_Wald(object, effects = effects, rhs = rhs, df = df, ci = ci, 
                           type.object = type.object,
                           transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
    }
    
    ## ** export
    class(out) <- append("anova_lmm",class(out))
    return(out)
}

## * .anova_Wald
.anova_Wald <- function(object, effects, rhs, df, ci, 
                        type.object,
                        transform.sigma, transform.k, transform.rho, transform.names){
    
    ## ** normalized user input
    terms.mean <- attr(stats::terms(object$formula$mean.design),"term.labels")
    subeffect <- NULL
    if(length(effects)==1){
        if(identical(effects,"all")){
            effects <- c("mean","variance","correlation")
        }else if(grepl("^mean_",effects)){
            iLabels <- attr(stats::terms(object$formula$mean.design),"term.labels")
            if(any(effects == paste0("mean_",iLabels))){
                subeffect <- iLabels[effects == paste0("mean_",iLabels)]
                effects <- "mean"
            }
        }else if(grepl("^variance_",effects) && !is.null(object$design$vcov$X$var)){
            iLabels <- attr(stats::terms(object$formula$var.design),"term.labels")
            if(any(effects == paste0("variance_",iLabels))){
                subeffect <- iLabels[effects == paste0("variance_",iLabels)]
                effects <- "variance"
            }
        }else if(grepl("^cor_",effects) && !is.null(object$design$vcov$X$cor)){
            iLabels <- attr(stats::terms(object$formula$cor.design),"term.labels")
            if(any(effects == paste0("correlation_",iLabels))){
                subeffect <- iLabels[effects == paste0("correlation_",iLabels)]
                effects <- "correlation"
            }
        }
    }

    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho)
    
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho

    name.coef <- names(stats::coef(object, effects = "all"))
    if(is.matrix(effects)){
        ## try to re-size the matrix if necessary
        if(NCOL(effects)!=length(name.coef)){
            if(is.null(colnames(effects))){
                stop("Argument \'effect\' should have column names when a matrix. \n")
            }
            if(any(duplicated(colnames(effects)))){
                stop("Argument \'effect\' should not have duplicated column names when a matrix. \n")
            }
            if(any(colnames(effects) %in% name.coef == FALSE)){
                stop("Argument \'effect\' should not have column names matching the coefficient names when a matrix. \n")
            }
            effects.save <- effects
            effects <- matrix(0, nrow = NROW(effects.save), ncol = length(name.coef), dimnames = list(rownames(effects.save),name.coef))
            effects[,colnames(effects.save)] <- effects.save
        }
        if(is.null(rhs)){
            rhs <- rep(0, NROW(effects))
        }
        ## run glht
        out.glht <- try(multcomp::glht(object, linfct = effects, rhs = rhs,
                                       coef. = function(iX){coef.lmm(iX, effects = "all")},
                                       vcov. = function(iX){vcov.lmm(iX, effects = "all")}),
                        silent = TRUE)
        if(inherits(out.glht,"try-error")){
            stop("Possible mispecification of the argument \'effects\' as running mulcomp::glht lead to the following error: \n",
                 out.glht)
        }
        out <- list(all = NULL)
        ls.nameTerms <- list(all = NULL)
        ls.nameTerms.num <- list(all = 1)
        ls.contrast <- list(all = out.glht$linfct)
        ls.null  <- list(all = out.glht$rhs)        

    }else if(all(tolower(effects) %in% c("mean","fixed","variance","correlation"))){
        
        if(transform.k %in% c("sd","var","logsd","logvar")){
            stop("Cannot use \'transform.rho\' equal \"sd\", \"var\", \"logsd\", or \"logvar\". \n",
                 "anova does not handle tests where the null hypothesis is at a boundary of the support of a random variable. \n")
        }
        if(transform.rho %in% c("cov")){
            stop("Cannot use \'transform.rho\' equal \"cov\". \n",
                 "anova does not handle tests where the null hypothesis is at a boundary of the support of a random variable. \n")
        }
        
        effects <- match.arg(effects, c("mean","fixed","variance","correlation"), several.ok = TRUE)
        effects[effects=="fixed"] <- "mean"
        
        out <- list()
        ls.assign <- list()
        ls.nameTerms <- list()
        ls.contrast <- list()
        ls.null <- list()
        if("mean" %in% effects){
            out <- c(out,list(mean = NULL))
            ls.assign$mean <- attr(object$design$mean,"assign")
            ls.nameTerms$mean <- attr(stats::terms(object$formula$mean.design),"term.labels")
            ls.contrast <- c(ls.contrast,list(mean = NULL))
            null.mean <- 0
            ls.null$mean <- rep(null.mean,length(ls.nameTerms$mean))            
        }
        if("variance" %in% effects){
            out <- c(out,list(variance = NULL))
            ls.assign$variance <- attr(object$design$vcov$X$var,"assign")            
            ls.nameTerms$variance <- attr(stats::terms(object$formula$var.design),"term.labels")
            if(!is.na(object$design$vcov$name$strata)){
                ls.assign$variance[ls.nameTerms$variance[ls.assign$variance]==object$design$vcov$name$strata] <- 0
            }
            ls.contrast <- c(ls.contrast,list(variance = NULL))
            null.variance <- switch(transform.k,
                                    "none" = 1,
                                    "square" = 1,
                                    "log" = 0,
                                    "logsquare" = 0)
            ls.null$variance <- rep(null.variance,length(ls.nameTerms$variance))
        }
        if("correlation" %in% effects){
            out <- c(out,list(correlation = NULL))
            ls.assign$correlation <- attr(object$design$vcov$X$cor,"assign")
            ls.nameTerms$correlation <- if(!is.null(ls.assign$correlation)){object$time$var}else{NULL}
            ls.contrast <- c(ls.contrast,list(correlation = NULL))
            null.correlation <- switch(transform.rho,
                                       "none" = 0,
                                       "atanh" = 0,
                                       "cov" = 0)
            ls.null$correlation <- rep(null.correlation,length(ls.nameTerms$correlation))
        }
        ls.nameTerms.num <- lapply(ls.nameTerms, function(iName){as.numeric(factor(iName, levels = iName))})
        
    }else if(all(grepl("=",effects)==FALSE)){
        stop("Incorrect argument \'effects\': can be \"mean\", \"variance\", \"correlation\", \"all\", \n",
             "or something compatible with the argument \'linfct\' of multcomp::glht. \n ")
    }else{
        
        ## run glht
        out.glht <- try(multcomp::glht(object, linfct = effects,
                                       coef. = function(iX){coef.lmm(iX, effects = "all")},
                                       vcov. = function(iX){vcov.lmm(iX, effects = "all")}),
                        silent = TRUE)
        newname.coef <- names(stats::coef(object, effects = "all"))
        
        if(inherits(out.glht,"try-error")){
            test.reparametrize <- grepl("log", c(object$reparametrize$transform.sigma,object$reparametrize$transform.k)) || grep("atanh", object$reparametrize$transform.rho)
            
            ## restaure untransformed parametrization (glht does not handle log(k). or atanh(cor))
            if(test.reparametrize){
                object2 <- object
                index.var <- which(object$param$type %in% c("sigma","k","rho"))
                object2$reparametrize <- .reparametrize(p = object$param$value[index.var], type = object$param$type[index.var], strata = object$param$strata[index.var], time.levels = object$time$levels,
                                                        time.k = object$design$param$time.k, time.rho = object$design$param$time.rho,
                                                        Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                                        transform.sigma = "none",
                                                        transform.k = "none",
                                                        transform.rho = "none",
                                                        transform.names = TRUE)
                out.glht <- try(multcomp::glht(object2, linfct = effects,
                                               coef. = function(iX){coef.lmm(iX, effects = "all")},
                                               vcov. = function(iX){vcov.lmm(iX, effects = "all")}), silent = TRUE)
                if(inherits(out.glht,"try-error")){
                    stop("Possible mispecification of the argument \'effects\' as running mulcomp::glht lead to the following error: \n",
                         out.glht)
                }
                oldname.coef <- colnames(out.glht$linfct)
                newname.hypo <- rownames(out.glht$linfct)
                for(iSub in 1:length(oldname.coef)){
                    newname.hypo <- gsub(pattern = oldname.coef[iSub], replacement = newname.coef[iSub], x = newname.hypo, fixed = TRUE)
                }
                dimnames(out.glht$linfct) <- list(newname.hypo,newname.coef)
                message("The coefficient names have been transformed but not the null hypotheses. \n")
            }else{
                stop("Possible mispecification of the argument \'effects\' as running mulcomp::glht lead to the following error: \n",
                     out.glht)
            }
        }
        out <- list(all = NULL)
        ls.nameTerms <- list(all = NULL)
        ls.nameTerms.num <- list(all = 1)
        ls.contrast <- list(all = out.glht$linfct)
        ls.null  <- list(all = out.glht$rhs)        
    }
    type.information <- attr(object$information,"type.information")    
    type.object <- match.arg(type.object, c("lmm","gls"))

    ## ** prepare
    if(type.object=="gls"){
        if(length(object$gls)==1){
            return(anova(object$gls))
        }else{
            return(lapply(object$gls, anova))
        }
    }

    param <- coef(object, effects = "all",
                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
    newname <- names(coef(object, effects = "all",
                          transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names))
    name.param <- names(param)
    n.param <- length(param)
    vcov.param <- vcov(object, df = df*2, effects = "all",
                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
    dVcov.param <- attr(vcov.param,"dVcov")
    if(type.information != "observed"){
        warning("when using REML with expected information, the degree of freedom of the F-statistic may depend on the parametrisation of the variance parameters. \n")
    }

    ## ** F-tests
    type <- names(out)
    for(iType in type){ ## iType <- "correlation"
        ## skip empty type
        if(length(ls.nameTerms.num[[iType]])==0 || (is.null(ls.contrast[[iType]]) && all(ls.assign[[iType]]==0))){ next }

        iParam <- coef(object, effects = iType,
                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        name.iParam <- names(iParam)

        iLs <- lapply(ls.nameTerms.num[[iType]], function(iTerm){ ## iTerm <- 1
            ## *** contrast matrix
            if(is.null(ls.contrast[[iType]])){
                if(all(ls.assign[[iType]]!=iTerm)){return(NULL)}
                if(!is.null(subeffect) && ls.nameTerms$mean[iTerm]!=subeffect){return(NULL)}
                iIndex.param <- which(ls.assign[[iType]]==iTerm)
                iN.hypo <- length(iIndex.param)
                iNull <- rep(ls.null[[iType]][iTerm],iN.hypo)
                iName.hypo <- paste(paste0(name.iParam[iIndex.param],"==",iNull), collapse = ", ")
                iC <- matrix(0, nrow = iN.hypo, ncol = n.param, dimnames = list(name.iParam[iIndex.param], newname))
                if(length(iIndex.param)==1){
                    iC[name.iParam[iIndex.param],name.iParam[iIndex.param]] <- 1
                }else{
                    diag(iC[name.iParam[iIndex.param],name.iParam[iIndex.param]]) <- 1
                }
                iC2 <- iC
                colnames(iC2) <- name.param[match(colnames(iC),newname)]
            }else{
                iC <- ls.contrast[[iType]]
                iC2 <- iC
                iN.hypo <- NROW(iC)
                iNull <- ls.null[[iType]]
                iName.hypo  <- paste(paste0(rownames(iC),"==",iNull),collapse=", ")
            }

            ## *** statistic
            iC.vcov.C_M1 <- try(solve(iC %*% vcov.param %*% t(iC)), silent = TRUE)
            if(inherits(iC.vcov.C_M1,"try-error")){
                stop("Could not invert the covariance matrix for the proposed contrast. \n")
            }
            iStat <- as.double(t(iC %*% param - iNull) %*% iC.vcov.C_M1 %*% (iC %*% param - iNull))

            ## *** degree of freedom
            if(df){
                svd.tempo <- eigen(iC.vcov.C_M1)
                D.svd <- diag(svd.tempo$values, nrow = iN.hypo, ncol = iN.hypo)
                P.svd <- svd.tempo$vectors
                contrast.svd <- sqrt(D.svd) %*% t(P.svd) %*% iC
                colnames(contrast.svd) <- name.param

                iNu_m <- dfSigma(contrast = contrast.svd,
                                 vcov = vcov.param,
                                 dVcov = dVcov.param,
                                 keep.param = name.param)
                
                iEQ <- sum(iNu_m/(iNu_m - 2))
                iDf <- 2 * iEQ/(iEQ - iN.hypo)
            }else{
                iDf <- Inf
            }

            ## *** confidence interval
            if(ci){
                if(df){
                    ci.df <-  .dfX(X.beta = iC2, vcov.param = vcov.param, dVcov.param = dVcov.param)
                }else{
                    ci.df <- Inf
                }
                CI <- data.frame(estimate = as.double(iC %*% param),
                                 se = sqrt(diag(iC %*% vcov.param %*% t(iC))),
                                 df = ci.df,
                                 statistic = NA,
                                 lower = NA,
                                 upper = NA,
                                 null = iNull,
                                 p.value = NA)
                CI$statistic <- CI$estimate/CI$se
                rownames(CI) <- rownames(iC)
                if(!is.null(names(effects))){
                    indexName <- intersect(which(names(effects)!=""),which(!is.na(names(effects))))
                    rownames(CI)[indexName] <- names(effects)[indexName]
                }
                CI.glht <- multcomp::glht(object, linfct = iC, rhs = iNull, df = ceiling(max(ci.df)),
                                          coef. = function(iX){coef.lmm(iX, effects = "all")},
                                          vcov. = function(iX){vcov.lmm(iX, effects = "all")})
            }else{
                CI <- NULL
                CI.glht <- NULL
            }

            ## *** test
            iRes <- data.frame("null" = iName.hypo,
                               "statistic" = iStat/iN.hypo,
                               "df.num" = iN.hypo,
                               "df.denom" = iDf,
                               "p.value" = 1 - stats::pf(iStat/iN.hypo, df1 = iN.hypo, df2 = iDf))
            attr(iRes, "contrast") <- iC
            attr(iRes, "CI") <- CI
            attr(iRes, "glht") <- CI.glht
            return(iRes)
            
        })
        if(!is.null(subeffect)){
            names(iLs) <- ls.nameTerms[[iType]]
            iLs <- iLs[subeffect]
        }else if(!is.null(ls.nameTerms[[iType]])){
            names(iLs) <- ls.nameTerms[[iType]]
        }

        out[[iType]] <- do.call(rbind, iLs)
        attr(out[[iType]], "contrast") <- lapply(iLs,attr,"contrast")
        attr(out[[iType]], "CI") <- lapply(iLs,attr,"CI")
        attr(out[[iType]], "glht") <- lapply(iLs,attr,"glht")
        
    }

    ## ** export
    attr(out, "test") <- "Wald"
    return(out)
}

## * .anova_LRT
.anova_LRT <- function(object1,object2){

    ## ** normalize user input
    if(all(names(coef(object1)) %in% names(coef(object2)))){
        objectH1 <- object2
        objectH0 <- object1
    }else if(all(names(coef(object2)) %in% names(coef(object1)))){
        objectH1 <- object1
        objectH0 <- object2
    }else{
        stop("One model must be nest in the other model to perform a likelihood ratio test. \n")
    }

    paramH1 <-  names(objectH1$param$type)
    typeH1 <-  objectH1$param$type
    paramH0 <-  names(objectH0$param$type)
    typeH0 <-  objectH0$param$type
    test.X <- identical(objectH0$design$mean[,paramH0[typeH0=="mu"],drop=FALSE], objectH1$design$mean[,paramH0[typeH0=="mu"],drop=FALSE])
    if(test.X==FALSE){
        stop("Mismatch between the design matrices for the mean of the two models - one should be nested in the other. \n")
    }

    paramTest <- setdiff(paramH1,paramH0)
    if(objectH1$method.fit!=objectH1$method.fit){
        stop("The two models should use the same type of objective function for the likelihood ratio test to be valid. \n")
    }
    if(objectH1$method.fit=="REML" && any(typeH1[paramTest]=="mu")){
        stop("Cannot test mean parameters when the objective function is REML. \n")
    }

    ## ** LRT
    out <- data.frame(null = paste(paste0(paramTest,"==0"), collapse = ", "),
                      logLikH1 = stats::logLik(objectH1),
                      logLikH0 = stats::logLik(objectH0),
                      statistic = NA,
                      df = length(paramTest),
                      p.value = NA)
    out$statistic <- 2*(out$logLikH1 - out$logLikH0)
    out$p.value <- 1 - stats::pchisq(out$statistic, df = out$df)
               
    ## ** export
    attr(out, "test") <- "LRT"
    return(out)
}

## * confint.anova_lmm
##' @rdname anova
##' @export
confint.anova_lmm <- function(object, parm, level = 0.95, method = "single-step", ...){

    ## ** normalize user input
    if(attr(object,"test") == "LRT"){
        message("No confidence interval available for likelihood ratio tests.")
        return(NULL)
    }
    if(!missing(parm)){
        stop("Argument \'parm\' is not used - only there for compatibility with the generic method. \n")
    }
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    alpha <- 1-level
    method <- match.arg(method, c("none","bonferroni","single-step"))

    ## ** extract info and compute CI
    out <- lapply(object, function(iO){ ## iO <- object[[1]]
        iTable <- attr(iO,"CI")
        if(is.null(iTable) || all(sapply(iTable,is.null))){return(NULL)}
        iOut <- stats::setNames(vector(mode = "list", length = length(iTable)),names(iTable))

        for(iTest in 1:length(iTable)){ ## iTest <- 1
            iOut[[iTest]] <- iTable[[iTest]]
            iOut[[iTest]]$df <- pmax(iOut[[iTest]]$df, options$min.df)
            if(method == "none" || NROW(iOut[[iTest]])==1){
                iOut[[iTest]]$lower <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * stats::qt(alpha/2, df = iOut[[iTest]]$df)
                iOut[[iTest]]$upper <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * stats::qt(1-alpha/2, df = iOut[[iTest]]$df)
                iOut[[iTest]]$p.value <- 2*(1-stats::pt( abs((iOut[[iTest]]$estimate-iOut[[iTest]]$null) / iOut[[iTest]]$se), df = iOut[[iTest]]$df))
            }else if(method == "single-step"){
                iGlht <- attr(iO,"glht")[[iTest]]
                iCi <- confint(iGlht)
                iOut[[iTest]]$lower <- iCi$confint[,"lwr"]
                iOut[[iTest]]$upper <- iCi$confint[,"upr"]
                iOut[[iTest]]$p.value <- summary(iGlht, test = multcomp::adjusted("single-step"))$test$pvalues
                iOut[[iTest]]$df <- iGlht$df
            }else if(method == "bonferroni"){
                p <- NROW(iOut[[iTest]])
                iOut[[iTest]]$lower <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * stats::qt(alpha/(2*p), df = iOut[[iTest]]$df)
                iOut[[iTest]]$upper <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * stats::qt(1-alpha/(2*p), df = iOut[[iTest]]$df)
                iOut[[iTest]]$p.value <- pmin(1,2*p*(1-stats::pt( abs((iOut[[iTest]]$estimate-iOut[[iTest]]$null) / iOut[[iTest]]$se), df = iOut[[iTest]]$df)))
            }
        }
        return(iOut)
    })
    return(out)
}

## * print.anova_lmm
##' @rdname anova
##' @export
print.anova_lmm <- function(x, level = 0.95, method = "single-step", print.null = FALSE, ...){

    
    if(attr(x,"test")=="Wald"){    
        type <- names(x)
        ci <- stats::confint(x, level = level, method = method)
        for(iType in type){

            if(is.null(x[[iType]])){next}

            if(!print.null){
                x[[iType]][["null"]] <- NULL
            }
            iNoDf <- is.infinite(x[[iType]]$df.denom)
            txt.test <- ifelse(all(iNoDf),"Chi-square test","F-test")
            if(iType == "all"){
                cat("                     ** User-specified hypotheses ** \n", sep="")
                cat(" - ",txt.test,"\n", sep="")
                print(x[[iType]], row.names = FALSE)
            }else{
                cat("                     ** ",iType," coefficients ** \n", sep="")
                cat(" - ",txt.test,"\n",sep="")
                print(x[[iType]])
            }
            if(!is.null(ci[[iType]])){
                options <- LMMstar.options()
                if(all(sapply(ci[[iType]],NROW)==1)){ ## always only one hypothesis in each global test
                    cat("\n - P-values and confidence interval \n", sep="")
                }else if(length(ci[[iType]])==1){ ## only one global test
                    cat("\n - P-values and confidence interval (adjusted for multiplicity) \n", sep="")
                }else{
                    cat("\n - P-values and confidence interval (adjusted for multiplicity within each global test) \n", sep="")
                }
                print(do.call(rbind,unname(ci[[iType]]))[,options$columns.anova])
            }
            cat("\n")
        }
    }else if(attr(x,"test")=="LRT"){
        cat(" - Likelihood ratio test \n")
        x.print <- as.data.frame(x)
        if(print.null==FALSE){x.print[["null"]] <- NULL}
        print(x.print)
    }

    return(invisible(NULL))
}
## * dfSigma
##' @title Degree of Freedom for the Chi-Square Test
##' @description Computation of the degrees of freedom of the chi-squared distribution
##' relative to the model-based variance. Copied of lavaSearch2:::dfSigmaRobust.
##' @noRd
##' 
##' @param contrast [numeric vector] the linear combination of parameters to test
##' @param vcov [numeric matrix] the variance-covariance matrix of the parameters.
##' @param dVcov [numeric array] the first derivative of the variance-covariance matrix of the parameters.
##' @param keep.param [character vector] the name of the parameters with non-zero first derivative of their variance parameter.
##' 
dfSigma <- function(contrast, vcov, dVcov, keep.param){
    ## iLink <- "LogCau~eta"
    C.vcov.C <- rowSums(contrast %*% vcov * contrast) ## variance matrix of the linear combination
    ## C.vcov.C - vcov[iLink,iLink]

    C.dVcov.C <- sapply(keep.param, function(x){
        rowSums(contrast %*% dVcov[,,x] * contrast)
    })
    ## C.dVcov.C - dVcov[iLink,iLink,]
    numerator <- 2 *(C.vcov.C)^2
    ## numerator - 2*vcov[iLink,iLink]^2
    denom <- rowSums(C.dVcov.C %*% vcov[keep.param,keep.param,drop=FALSE] * C.dVcov.C)
    ## denom - t(dVcov[iLink,iLink,]) %*% vcov[keep.param,keep.param,drop=FALSE] %*% dVcov[iLink,iLink,]
    df <- numerator/denom
    return(df)
}

##----------------------------------------------------------------------
### anova.R ends here
