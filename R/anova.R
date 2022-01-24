### anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:38) 
## Version: 
## Last-Updated: jan 24 2022 (16:37) 
##           By: Brice Ozenne
##     Update #: 655
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
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors. 
##' @param rhs [numeric vector] the right hand side of the hypothesis. Only used when the argument effects is a matrix.
##' @param ci [logical] Should a confidence interval be output for each hypothesis?
##' @param level [numeric, 0-1] nominal coverage of the confidence intervals.
##' @param print.null [logical] should the null hypotheses be printed in the console?
##' @param df [logical] Should a F-distribution be used to model the distribution of the Wald statistic. Otherwise a chi-squared distribution is used.
##' @param transform.sigma,transform.k,transform.rho,transform.names are passed to the \code{vcov} method. See details section in \code{\link{coef.lmm}}.
##' @param columns [character vector] Columns to be output. Can be any of \code{"estimate"}, \code{"se"}, \code{"df"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param transform [function] function to backtransform the estimates, standard errors, null hypothesis, and the associated confidence intervals
##' (e.g. \code{exp} if the outcomes have been log-transformed).
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
##' if(require(multcomp)){
##' amod <- lmm(breaks ~ tension, data = warpbreaks)
##' e.glht <- glht(amod, linfct = mcp(tension = "Tukey"))
##' summary(e.glht, test = Chisqtest()) ## 0.000742
##'
##' print(anova(amod, effect = mcp(tension = "Tukey"), df = FALSE), print.null = TRUE)
##' 
##' anova(amod, effect = mcp(tension = "Tukey"), ci = TRUE)
##' }

## * anova.lmm (code)
##' @rdname anova
##' @export
anova.lmm <- function(object, effects = NULL, robust = FALSE, rhs = NULL, df = !is.null(object$df), ci = FALSE, 
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
        out <- .anova_Wald(object, effects = effects, robust = robust, rhs = rhs, df = df, ci = ci, 
                           transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
    }
    
    ## ** export
    out$call <- match.call()
    class(out) <- append("anova_lmm",class(out))
    return(out)
}

## * .anova_Wald
.anova_Wald <- function(object, effects, robust, rhs, df, ci, 
                        transform.sigma, transform.k, transform.rho, transform.names){
    
    ## ** normalized user input
    terms.mean <- attr(stats::terms(object$formula$mean.design),"term.labels")
    subeffect <- NULL
    if(!inherits(effects,"mcp") && length(effects)==1){
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
    if(inherits(effects,"mcp")){        
        out.glht <- try(multcomp::glht(object, linfct = effects), ## only used for generating contrast matrix
                        silent = TRUE)
        if(inherits(out.glht,"try-error")){
            stop("Possible mispecification of the argument \'effects\' as running mulcomp::glht lead to the following error: \n",
                 out.glht)
        }
        out <- list(all = NULL)
        ls.nameTerms <- list(all = NULL)
        ls.nameTerms.num <- list(all = 1)
        ls.contrast <- list(all = matrix(0, nrow = NROW(out.glht$linfct), ncol = length(name.coef), dimnames = list(rownames(out.glht$linfct),name.coef)))
        ls.contrast$all[,colnames(out.glht$linfct)] <- out.glht$linfct
        ls.null  <- list(all = out.glht$rhs)        
    }else if(is.matrix(effects)){
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
        out.glht <- try(multcomp::glht(object, linfct = effects, rhs = rhs,  ## only used for generating contrast matrix
                                       coef. = function(iX){coef.lmm(iX, effects = "all")},
                                       vcov. = function(iX){vcov.lmm(iX, robust = robust, effects = "all")}),
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
    }else{ ## symbolic definition of effects using equations (characters)
        
        ## run glht
        out.glht <- try(multcomp::glht(object, linfct = effects,  ## only used for generating contrast matrix
                                       coef. = function(iX){coef.lmm(iX, effects = "all")},
                                       vcov. = function(iX){vcov.lmm(iX, robust = robust, effects = "all")}),
                        silent = TRUE)
        newname.coef <- names(stats::coef(object, effects = "all"))
        
        if(inherits(out.glht,"try-error")){
            test.reparametrize <- grepl("log", c(object$reparametrize$transform.sigma,object$reparametrize$transform.k)) || grep("atanh", object$reparametrize$transform.rho)
            
            ## restaure untransformed parametrization (glht does not handle log(k). or atanh(cor))
            if(test.reparametrize){
                object2 <- object
                index.var <- which(object$param$type %in% c("sigma","k","rho"))
                object2$reparametrize <- .reparametrize(p = object$param$value[index.var],
                                                        type = object$param$type[index.var], strata = object$param$strata[index.var], 
                                                        time.k = object$design$param$time.k, time.rho = object$design$param$time.rho,
                                                        name2sd = stats::setNames(object$design$vcov$param$name2,object$design$vcov$param$name),
                                                        Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                                        transform.sigma = "none",
                                                        transform.k = "none",
                                                        transform.rho = "none",
                                                        transform.names = TRUE)
                out.glht <- try(multcomp::glht(object2, linfct = effects,  ## only used for generating contrast matrix
                                               coef. = function(iX){coef.lmm(iX, effects = "all")},
                                               vcov. = function(iX){vcov.lmm(iX, robust = robust, effects = "all")}), silent = TRUE)
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

    ## ** prepare
    param <- coef(object, effects = "all",
                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
    newname <- names(coef(object, effects = "all",
                          transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names))
    name.param <- names(param)
    n.param <- length(param)
    vcov.param <- vcov(object, df = df*2, effects = "all", robust = robust,
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
                iC.uni <- iC
                colnames(iC.uni) <- name.param[match(colnames(iC),newname)]
            }else{
                iC <- ls.contrast[[iType]]
                iC.uni <- iC
                iN.hypo <- NROW(iC)
                iNull <- ls.null[[iType]]
                iName.hypo  <- paste(paste0(rownames(iC),"==",iNull),collapse=", ")
            }

            ## *** statistic
            outSimp <- .simplifyContrast(iC, iNull) ## remove extra lines
            iC.vcov.C_M1 <- try(solve(outSimp$C %*% vcov.param %*% t(outSimp$C)), silent = TRUE)
            if(inherits(iC.vcov.C_M1,"try-error")){
                iStat <- NA
                iDf <- c(iN.hypo,Inf)
                attr(iStat,"error") <- "\n  Could not invert the covariance matrix for the proposed contrast."
            }else{
                iStat <- as.double(t(outSimp$C %*% param - outSimp$rhs) %*% iC.vcov.C_M1 %*% (outSimp$C %*% param - outSimp$rhs))/outSimp$dim 
                iDf <- c(outSimp$dim,Inf)
                if(outSimp$rm>0){
                    iName.hypo <- paste(paste0(rownames(outSimp$C),"==",outSimp$rhs),collapse=", ")
                }
            }

            ## *** degree of freedom
            if(df && !inherits(iC.vcov.C_M1,"try-error")){

                svd.tempo <- eigen(iC.vcov.C_M1)
                D.svd <- diag(svd.tempo$values, nrow = outSimp$dim, ncol = outSimp$dim)
                P.svd <- svd.tempo$vectors
                contrast.svd <- sqrt(D.svd) %*% t(P.svd) %*% outSimp$C
                colnames(contrast.svd) <- name.param

                iNu_m <- dfSigma(contrast = contrast.svd,
                                 vcov = vcov.param,
                                 dVcov = dVcov.param,
                                 keep.param = name.param)
                
                iEQ <- sum(iNu_m/(iNu_m - 2))
                iDf[2] <- 2 * iEQ/(iEQ - outSimp$dim)
            }

            ## *** confidence interval
            if(ci){
                if(df){
                    ci.df <-  .dfX(X.beta = iC.uni, vcov.param = vcov.param, dVcov.param = dVcov.param)
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
                                 p.value = NA,
                                 stringsAsFactors = FALSE)
                CI$statistic <- CI$estimate/CI$se
                rownames(CI) <- rownames(iC)
                if(!is.null(names(effects)) && !inherits(effects,"mcp")){
                    indexName <- intersect(which(names(effects)!=""),which(!is.na(names(effects))))
                    rownames(CI)[indexName] <- names(effects)[indexName]
                }
                CI.glht <- multcomp::glht(object, linfct = iC, rhs = iNull, df = ceiling(max(ci.df)),
                                          coef. = function(iX){coef.lmm(iX, effects = "all")},
                                          vcov. = function(iX){vcov.lmm(iX, robust = robust, effects = "all")})
                attr(CI.glht$vcov,"robust") <- robust
            }else{
                CI <- NULL
                CI.glht <- NULL
            }

            ## *** test
            iRes <- data.frame("null" = iName.hypo,
                               "statistic" = iStat,
                               "df.num" = iDf[1],
                               "df.denom" = iDf[2],
                               "p.value" = 1 - stats::pf(iStat, df1 = iDf[1], df2 = iDf[2]),
                               stringsAsFactors = FALSE)
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
    if(NROW(objectH0$design$mean)!=NROW(objectH1$design$mean)){
        stop("Mismatch between the design matrices for the mean of the two models - could be due to missing data. \n",
             "Different number of rows: ",NROW(objectH0$design$mean)," vs. ",NROW(objectH1$design$mean),".\n")
    }

    test.X <- identical(objectH0$design$mean[,paramH0[typeH0=="mu"],drop=FALSE], objectH1$design$mean[,paramH0[typeH0=="mu"],drop=FALSE])
    if(test.X==FALSE){
        stop("Mismatch between the design matrices for the mean of the two models - one should be nested in the other. \n")
    }

    paramTest <- setdiff(paramH1,paramH0)
    if(objectH1$method.fit!=objectH1$method.fit){
        stop("The two models should use the same type of objective function for the likelihood ratio test to be valid. \n")
    }
    if(objectH1$method.fit=="REML" && any(typeH1[paramTest]=="mu")){
        objectH0$call$method.fit <- "ML"
        objectH1$call$method.fit <- "ML"
        message("Cannot use a likelihood ratio test to compare mean parameters when the objective function is REML. \n",
                "Will re-estimate the model via ML and re-run the likelihood ratio test. \n")
        return(anova(eval(objectH0$call),eval(objectH1$call)))
    }

    ## ** LRT
    out <- data.frame(null = paste(paste0(paramTest,"==0"), collapse = ", "),
                      logLikH1 = stats::logLik(objectH1),
                      logLikH0 = stats::logLik(objectH0),
                      statistic = NA,
                      df = length(paramTest),
                      p.value = NA,
                      stringsAsFactors = FALSE)
    out$statistic <- 2*(out$logLikH1 - out$logLikH0)
    out$p.value <- 1 - stats::pchisq(out$statistic, df = out$df)
               
    ## ** export
    attr(out, "test") <- "LRT"
    return(out)
}

## * confint.anova_lmm
##' @title Confidence Intervals for Multivariate Wald Tests
##' @description Compute confidence intervals for linear hypothesis tests, possibly with adjustment for multiple comparisons.
##' @name anova
##' 
##' @param object a \code{anova_lmm} object
##' @param level [numeric, 0-1] nominal coverage of the confidence intervals.
##' @param method [character] type of adjustment for multiple comparisons: one of \code{"none"}, \code{"bonferroni"}, \code{"single-step"}.
##' @param parm Not used. For compatibility with the generic method.
##' @param ... Not used. For compatibility with the generic method.
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

## * rbind.anova_lmm
##' @title Linear Hypothesis Testing Across Linear Mixed Models
##' @description Linear hypothesis testing accross linear mixed model.
##'
##' @param model a \code{anova_lmm} object (output of \code{anova} applied to a \code{lmm} object)
##' @param linfct a \code{anova_lmm} object (output of \code{anova} applied to a \code{lmm} object)
##' @param ...  possibly other \code{anova_lmm} objects
##' @param sep [character] character used to separate the outcome and the covariate when naming the tests.
##' 
##' @export
rbind.anova_lmm <- function(model, linfct, ..., sep = ": "){
    default <- LMMstar.options()
    
    ## ** check user input
    if(!inherits(linfct,"anova_lmm")){
        stop("Argument \'anova_lmm\' should inherit from anova_lmm. \n")
    }
    dots <- list(...)
    if(any(sapply(dots,inherits,"anova_lmm")==FALSE)){
        stop("Extra arguments should inherit from anova_lmm. \n")
    }
    ls.object <- c(list(model, linfct),dots)
    if(any(sapply(ls.object,function(iO){"all" %in% names(iO)})==FALSE)){
        stop("All argument should correspond to user specified hypothesis, i.e. call anova with argument linfct. \n")
    }
    if(any(sapply(ls.object,function(iO){!is.null(attr(iO$all,"glht"))})==FALSE)){
        stop("All argument should contain a \"glht\" object, i.e. call anova with argument ci=TRUE. \n")
    }
    ls.glht <- lapply(ls.object, function(iO){attr(iO$all,"glht")[[1]]})
    
    ## ** Extract elements from anova object
    ls.C <- lapply(ls.glht,"[[","linfct")
    ls.rhs <- lapply(ls.glht,"[[","rhs")
    ls.coef <- lapply(ls.glht,"[[","coef")
    ls.lmm <- lapply(ls.glht,"[[","model")
    ls.df <- lapply(ls.glht,"[[","df")
    ls.alternative <- lapply(ls.glht,"[[","alternative")
    if(any(unlist(ls.alternative) != ls.alternative[[1]])){
        stop("Element \'alternative\' should take the same value for all glht objects. \n")
    }
    ls.robust <- lapply(ls.glht,function(iO){attr(iO$vcov,"robust")})
    if(any(unlist(ls.robust) != ls.robust[[1]])){
        stop("Element \'robust\' should take the same value for all glht objects. \n")
    }
    robust <- unique(unlist(ls.robust))
    ls.transform.sigma <- lapply(ls.object[[1]],function(iO){if(is.null(iO$call$transform.sigma)){default$transform.sigma}else{iO$call$transform.sigma}})
    if(any(unlist(ls.transform.sigma) != ls.transform.sigma[[1]])){
        stop("Element \'transform.sigma\' should take the same value for all glht objects. \n")
    }
    transform.sigma <- unique(unlist(ls.transform.sigma))
    ls.transform.k <- lapply(ls.object[[1]],function(iO){if(is.null(iO$call$transform.k)){default$transform.k}else{iO$call$transform.k}})
    if(any(unlist(ls.transform.k) != ls.transform.k[[1]])){
        stop("Element \'transform.k\' should take the same value for all glht objects. \n")
    }
    transform.k <- unique(unlist(ls.transform.k))
    ls.transform.rho <- lapply(ls.object[[1]],function(iO){if(is.null(iO$call$transform.rho)){default$transform.rho}else{iO$call$transform.rho}})
    if(any(unlist(ls.transform.rho) != ls.transform.rho[[1]])){
        stop("Element \'transform.rho\' should take the same value for all glht objects. \n")
    }
    transform.rho <- unique(unlist(ls.transform.rho))
    ls.method.fit <- lapply(ls.lmm,"[[","method.fit")
    if(any(unlist(ls.method.fit) != ls.method.fit[[1]])){
        stop("Element \'ls.method.fit\' should take the same value for all glht objects. \n")
    }
    method.fit <- unique(unlist(ls.method.fit))
    vec.outcome <- unname(unlist(sapply(ls.lmm,"[[","outcome")))

    names(ls.C) <- vec.outcome
    names(ls.coef) <- vec.outcome
    names(ls.lmm) <- vec.outcome

    ## ** extract iid
    ls.iid <- lapply(ls.lmm, function(iO){
        iid(iO, effects = if(method.fit=="REML"){"mean"}else{"all"}, robust = robust,
            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
        
    })
    names(ls.iid) <- vec.outcome
    
    ## ** build glht object
    out <- list(model = ls.lmm, robust = robust)

    out$linfct <- as.matrix(do.call(Matrix::bdiag, ls.C))
    rownames(out$linfct) <- unlist(lapply(vec.outcome, function(iO){paste0(iO,sep,rownames(ls.C[[iO]]))}))
    colnames(out$linfct) <- unlist(lapply(vec.outcome, function(iO){paste0(iO,sep,colnames(ls.C[[iO]]))}))
    
    out$rhs <- unlist(ls.rhs)
    out$coef <- unlist(lapply(vec.outcome, function(iO){stats::setNames(ls.coef[[iO]],paste0(iO,sep,names(ls.coef[[iO]])))}))

    out$vcov <- crossprod(do.call(cbind,ls.iid))
    rownames(out$vcov) <- unlist(lapply(vec.outcome, function(iO){paste0(iO,sep,colnames(ls.iid[[iO]]))}))
    colnames(out$vcov) <- unlist(lapply(vec.outcome, function(iO){paste0(iO,sep,colnames(ls.iid[[iO]]))}))

    if(method.fit=="REML"){
        if(any(abs(out$linfct[,setdiff(colnames(out$linfct),colnames(out$vcov))]>1e-10))){
            stop("Cannot test covariance structure across models when using REML. \n",
                 "Consider setting argument \'method.fit\' to \"ML\" when calling lmm. \n")
        }
        out$linfct <- out$linfct[,colnames(out$vcov),drop=FALSE]
        out$coef <- out$coef[colnames(out$vcov)]
    }
    
    out$df <- unlist(ls.df)
    out$alternative <- ls.alternative[[1]]

    ## ** export
    class(out) <- c("Manova_lmm","glht")
    return(out)
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

## * .simplifyContrast
## remove contrasts making the contrast matrix singular
.simplifyContrast <- function(object, rhs, tol = 1e-10, trace = TRUE){
    object.eigen <- eigen(tcrossprod(object))
    n.zero <- sum(abs(object.eigen$values) < tol)

    if(n.zero==0){return(list(C = object, rhs = rhs, dim = NROW(object), rm = 0))}
    
    keep.lines <- 1:NROW(object)
    for(iLine in NROW(object):1){ ## iLine <- 3
        iN.zero <- sum(abs(eigen(tcrossprod(object[setdiff(keep.lines,iLine),,drop=FALSE]))$values) < tol)
        if(iN.zero<n.zero){
            keep.lines <- setdiff(keep.lines,iLine)
            n.zero <- iN.zero
        }
        if(n.zero==0){
            
            if(trace){
                name.rm <- rownames(object)[-keep.lines]
                if(length(name.rm)==1){
                    message("Singular contrast matrix: contrast \"",name.rm,"\" has been removed. \n")
                }else{
                    message("Singular contrast matrix: contrasts \"",paste(name.rm,collapse= "\" \""),"\" have been removed. \n")
                }
            }

            return(list(C = object[keep.lines,,drop=FALSE], rhs = rhs[keep.lines], dim = length(keep.lines), rm = NROW(object)-length(keep.lines)))
        }
    }

    ## n.zero>0 so failure
    return(NULL)
}

##----------------------------------------------------------------------
### anova.R ends here
