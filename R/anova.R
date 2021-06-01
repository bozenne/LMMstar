### anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:38) 
## Version: 
## Last-Updated: jun  1 2021 (16:32) 
##           By: Brice Ozenne
##     Update #: 354
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * anova.lmm (documentation)
##' @title Multivariate Wald Tests For Linear Mixed Models
##' @description Perform a Wald test testing simultaneously several null hypotheses corresponding to linear combinations of the model paramaters. 
##' @name anova
##' 
##' @param object a \code{lmm} object.
##' @param effects [character] Should the Wald test be computed for all variables (\code{"all"}),
##' or only variables relative to the mean (\code{"mean"}),
##' or only variables relative to the variance structure (\code{"variance"}),
##' or only variables relative to the correlation structure (\code{"correlation"}).
##' Can also be use to specify linear combinations of coefficients, similarly to the \code{linfct} argument of the \code{multcomp::glht} function.
##' @param ci [logical] Should a confidence interval be output for each hypothesis?
##' @param conf.level [numeric, 0-1] nominal coverage of the confidence intervals.
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.
##' @param df [logical] Should a F-distribution be used to model the distribution of the Wald statistic. Otherwise a chi-squared distribution is used.
##' @param type.information,transform.sigma,transform.k,transform.rho,transform.names are passed to the \code{vcov} method. See details section in \code{\link{coef.lmm}}.
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
##' ## fit mixed model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##' 
##' ## chi-2 test
##' anova(eUN.lmm)
##' 
##' anova(eUN.lmm, effects = c("X1=0","X2+X5=10"), ci = TRUE)
##' 
##' ## F-test
##' \dontrun{
##' anova(eUN.lmm, df = TRUE)
##' }
##' 

## * anova.lmm (code)
##' @rdname anova
##' @export
anova.lmm <- function(object, effects = "all", df = !is.null(object$df), ci = FALSE, 
                      type.object = "lmm",
                      transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){
    
    options <- LMMstar.options()

    ## ** normalized user input    
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }

    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, options = options,
                            x.transform.sigma = NULL, x.transform.k = NULL, x.transform.rho = NULL)
    
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho

    if(all(tolower(effects) %in% c("mean","variance","correlation"))){        
        if(transform.k %in% c("sd","var","logsd","logvar")){
            stop("Cannot use \'transform.rho\' equal \"sd\", \"var\", \"logsd\", or \"logvar\". \n",
                 "anova does not handle tests where the null hypothesis is at a boundary of the support of a random variable. \n")
        }
        if(transform.rho %in% c("cov")){
            stop("Cannot use \'transform.rho\' equal \"cov\". \n",
                 "anova does not handle tests where the null hypothesis is at a boundary of the support of a random variable. \n")
        }
        
        effects <- match.arg(effects, c("mean","variance","correlation"), several.ok = TRUE)
        out <- list(mean=NULL,
                    variance=NULL,
                    correlation=NULL)
        ls.assign <- list(mean = attr(object$design$X.mean,"assign"),
                          variance = attr(object$design$X.var$var,"assign"),
                          correlation = attr(object$design$X.var$cor,"assign"))
        ls.nameTerms <- list(mean = attr(stats::terms(object$formula$mean.design),"term.labels"),
                             variance = if(!is.null(object$formula$var.design)){attr(stats::terms(object$formula$var.design),"term.labels")}else{NULL},
                             correlation = if(!is.null(ls.assign$correlation)){object$time$var}else{NULL})
        ls.nameTerms.num <- lapply(ls.nameTerms, function(iName){as.numeric(factor(iName, levels = iName))})
        ls.contrast  <- list(mean = NULL, variance = NULL, correlation = NULL)

        null.mean <- 0
        null.variance <- switch(transform.k,
                                "none" = 1,
                                "square" = 1,
                                "log" = 0,
                                "logsquare" = 0)
        null.correlation <- switch(transform.rho,
                                   "none" = 0,
                                   "atanh" = 0,
                                   "cov" = 0)
        ls.null  <- list(mean = rep(null.mean,length(ls.nameTerms$mean)),
                         variance = rep(null.variance,length(ls.nameTerms$variance)),
                         correlation = rep(null.correlation,length(ls.nameTerms$correlation)))
    }else if(all(grepl("=",effects)==FALSE)){
        stop("Incorrect argument \'effects\': can be \"mean\", \"variance\", \"correlation\", \"all\", \n",
             "or something compatible with the argument \'linfct\' of multcomp::glht. \n ")
    }else{
        out.glht <- try(multcomp::glht(object, linfct = effects), silent = TRUE)
        if(inherits(out.glht,"try-error")){
            stop("Possible mispecification of the argument \'effects\' as running mulcomp::glht lead to the following error: \n",
                 out.glht)
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
    for(iType in type){

        ## skip empty type
        if(length(ls.nameTerms.num[[iType]])==0 || (is.null(ls.contrast[[iType]]) && all(ls.assign[[iType]]==0))){ next }

        iParam <- coef(object, effects = iType,
                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        name.iParam <- names(iParam)

        iLs <- lapply(ls.nameTerms.num[[iType]], function(iTerm){ ## iTerm <- 1
            ## *** contrast matrix
            if(is.null(ls.contrast[[iType]])){
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
            }else{
                iC <- ls.contrast[[iType]]
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
                    iC2 <- iC
                    colnames(iC2) <- name.param[match(colnames(iC),newname)]
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
                CI.glht <- multcomp::glht(object, linfct = iC, rhs = iNull, df = ceiling(max(ci.df)))
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
        if(!is.null(ls.nameTerms[[iType]])){
            names(iLs) <- ls.nameTerms[[iType]]
        }
        out[[iType]] <- do.call(rbind, iLs)
        attr(out[[iType]], "contrast") <- lapply(iLs,attr,"contrast")
        attr(out[[iType]], "CI") <- lapply(iLs,attr,"CI")
        attr(out[[iType]], "glht") <- lapply(iLs,attr,"glht")
        
    }
    
    ## ** export
    class(out) <- append("anova_lmm",class(out))
    return(out)
}

## * confint.anova_lmm
##' @export
confint.anova_lmm <- function(object, parm, level = 0.95, method = "single-step", ...){

    ## ** normalize user input
    if(!missing(parm)){
        stop("Argument \'parm\' is not used - only there for compatibility with the generic method. \n")
    }
    dots <- list(...)
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
            
            if(method == "single-step"){
                iGlht <- attr(iO,"glht")[[iTest]]
                iCi <- confint(iGlht)
                iOut[[iTest]]$lower <- iCi$confint[,"lwr"]
                iOut[[iTest]]$upper <- iCi$confint[,"upr"]
                iOut[[iTest]]$p.value <- summary(iGlht, test = adjusted("single-step"))$test$pvalues
                iOut[[iTest]]$df <- iGlht$df
            }else if(method == "none"){
                iOut[[iTest]]$lower <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * qt(alpha/2, df = iOut[[iTest]]$df)
                iOut[[iTest]]$upper <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * qt(1-alpha/2, df = iOut[[iTest]]$df)
                iOut[[iTest]]$p.value <- 2*(1-stats::pt( abs((iOut[[iTest]]$estimate-iOut[[iTest]]$null) / iOut[[iTest]]$se), df = iOut[[iTest]]$df))
            }else if(method == "bonferroni"){
                p <- NROW(iOut[[iTest]])
                iOut[[iTest]]$lower <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * qt(alpha/(2*p), df = iOut[[iTest]]$df)
                iOut[[iTest]]$upper <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * qt(1-alpha/(2*p), df = iOut[[iTest]]$df)
                iOut[[iTest]]$p.value <- pmin(1,2*p*(1-pt( abs((iOut[[iTest]]$estimate-iOut[[iTest]]$null) / iOut[[iTest]]$se), df = iOut[[iTest]]$df)))
            }
        }
    return(iOut)
    })
return(out)
}

## * print.anova_lmm
##' @title Print Multivariate Wald Tests For Linear Mixed Models
##' 
##' @param x an object of class \code{anova_lmm}, i.e. output of anova function on a \code{lmm} object.
##' @param print.null [logical] should the null hypotheses be printed in the console?
##' 
##' @export
print.anova_lmm <- function(x, print.null = FALSE, ...){

    type <- names(x)
    ci <- confint(x)
    cat("\n")
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
            cat("\n - P-values and confidence interval (adjusted for multiplicity within each global test) \n", sep="")
            print(do.call(rbind,unname(ci[[iType]])))
        }
        cat("\n")
    }

    return(invisible(NULL))
}
## * dfSigma
##' @title Degree of Freedom for the Chi-Square Test
##' @description Computation of the degrees of freedom of the chi-squared distribution
##' relative to the model-based variance. Copied of lavaSearch2:::dfSigmaRobust.
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
