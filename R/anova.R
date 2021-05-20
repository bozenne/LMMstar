### anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:38) 
## Version: 
## Last-Updated: May 19 2021 (12:44) 
##           By: Brice Ozenne
##     Update #: 158
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * anova.lmm (code)
anova.lmm <- function(object, effects = "all", df = TRUE, print = TRUE, print.null = FALSE,
                      transform.sigma = "log", transform.k = "log", transform.rho = "atanh", transform.names = TRUE, ...){
    
    

    ## ** normalized user input    
    if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
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
                          variance = attr(object$design$X.var,"assign"))
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
    }else{
        out.glht <- try(glht(object, linfct = effects), silent = TRUE)
        if(inherits(out.glht,"try-error")){
            stop("Incorrect argument \'effects\': can be \"mean\", \"variance\", \"correlation\", \"all\", \n",
                 "or something compatible with the argument \'linfct\' of multcomp::glht. \n ")
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
                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
    name.param <- names(param)
    n.param <- length(param)

    vcov.param <- vcov(object, df = df*2, effects = "all",
                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
    dVcov.param <- attr(vcov.param,"dVcov")
    if(type.information != "observed"){
        warning("when using REML with expected information, the degree of freedom of the F-statistic may depend on the parametrisation of the variance parameters. \n")
    }

    ## ** F-tests
    type <- names(out)
    for(iType in type){

        ## skip empty type
        if(length(ls.nameTerms.num[[iType]])==0){ next }
        iParam <- coef(object, effects = iType,
                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        name.iParam <- names(iParam)

        iLs <- lapply(ls.nameTerms.num[[iType]], function(iTerm){
            ## *** contrast matrix
            if(is.null(ls.contrast[[iType]])){
                iIndex.param <- which(ls.assign[[iType]]==iTerm)
                iN.hypo <- length(iIndex.param)
                iNull <- rep(ls.null[[iType]][iTerm],iN.hypo)
                iName.hypo <- paste(paste0(name.iParam[iIndex.param],"==",iNull), collapse = ", ")
                iC <- matrix(0, nrow = iN.hypo, ncol = n.param, dimnames = list(name.iParam[iIndex.param], name.param))
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

                iNu_m <- dfSigma(contrast = sqrt(D.svd) %*% t(P.svd) %*% iC,
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
            if(iType == "all"){
                cat("F-test for user-specified linear hypotheses \n", sep="")
                print(printout[[iType]], row.names = FALSE)
            }else{
                cat("F-test for the ",iType," coefficients \n", sep="")
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
