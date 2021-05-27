### anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:38) 
## Version: 
## Last-Updated: May 27 2021 (16:32) 
##           By: Brice Ozenne
##     Update #: 217
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
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.
##' @param df [logical] Should a F-distribution be used to model the distribution of the Wald statistic. Otherwise a chi-squared distribution is used.
##' @param print [logical] should the output be printed in the console.
##' @param print.null [logical] should the null hypotheses be printed in the console.
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
##' anova(eUN.lmm, effects = c("X1=0","X2+X5=0"))
##' 
##' ## F-test
##' \dontrun{
##' anova(eUN.lmm, df = TRUE)
##' }
##' 

## * anova.lmm (code)
##' @export
anova.lmm <- function(object, effects = "all", df = !is.null(object$df), print = TRUE, print.null = FALSE, type.object = "lmm",
                      transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){
    
    options <- LMMstar.options()

    ## ** normalized user input    
    if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }

    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, options = options,
                            x.transform.sigma = NULL, x.transform.k = NULL, x.transform.rho = NULL,
                            backtransform.sigma = NULL, backtransform.k = NULL, backtransform.rho = NULL)
    
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
                          variance = attr(object$design$X.var$var,"assign"))
        ls.nameTerms <- list(mean = attr(terms(object$formula$mean.design),"term.labels"),
                             variance = if(!is.null(object$formula$var.design)){attr(terms(object$formula$var.design),"term.labels")}else{NULL})
        ls.nameTerms.num <- lapply(ls.nameTerms, function(iName){as.numeric(factor(iName, levels = iName))})
        ls.contrast  <- list(mean = NULL, variance = NULL)

        null.mean <- 0
        null.variance <- switch(transform.k,
                                "none" = 1,
                                "square" = 1,
                                "log" = 0,
                                "logsquare" = 0)
        ls.null  <- list(mean = rep(null.mean,length(ls.nameTerms$mean)), variance = rep(null.variance,length(ls.nameTerms$variance)))
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

            ## *** test
            iRes <- data.frame("null" = iName.hypo,
                               "statistic" = iStat/iN.hypo,
                               "df.num" = iN.hypo,
                               "df.denom" = iDf,
                               "p.value" = 1 - pf(iStat/iN.hypo, df1 = iN.hypo, df2 = iDf))
            attr(iRes, "contrast") <- iC
            return(iRes)
            
        })
        if(!is.null(ls.nameTerms[[iType]])){
            names(iLs) <- ls.nameTerms[[iType]]
        }
        out[[iType]] <- do.call(rbind, iLs)
            
        if(print){
            printout <- out
            if(!print.null){
                printout[[iType]][["null"]] <- NULL
            }
            txt.test <- ifelse(df>0,"F-test","Chi-square test")
            if(iType == "all"){
                cat(txt.test," for user-specified linear hypotheses \n", sep="")
                print(printout[[iType]], row.names = FALSE)
            }else{
                cat(txt.test," for the ",iType," coefficients \n", sep="")
                print(printout[[iType]])
            }
            cat("\n")
        }
    }
    
    ## ** export
    return(invisible(out))
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
