### anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:38) 
## Version: 
## Last-Updated: May 14 2021 (15:53) 
##           By: Brice Ozenne
##     Update #: 97
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * anova.lmm (code)
anova.lmm <- function(object, effects = "all", df = TRUE, print = TRUE, print.null = FALSE, ...){
    
    param <- coef(object, effects = "all")
    name.param <- names(param)
    p.param <- length(name.param)
    if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    if(all(tolower(effects) %in% c("mean","variance","correlation"))){
        effects <- match.arg(effects, c("mean","variance","correlation"), several.ok = TRUE)
        out <- list(mean=NULL,
                    variance=NULL,
                    correlation=NULL)
        ls.assign <- list(mean = attr(object$design$X.mean,"assign"),
                          variance = attr(object$design$X.var,"assign"))
        ls.nameTerms <- list(mean = attr(terms(object$formula$mean.design),"term.labels"),
                             variance = if(!is.null(object$formula$variance.design)){attr(terms(object$formula$variance.design),"term.labels")}else{NULL})
        ls.nameTerms.num <- lapply(ls.nameTerms, function(iName){as.numeric(factor(iName, levels = iName))})
        ls.contrast  <- list(mean = NULL, variance = NULL)
        ls.null  <- list(mean = rep(0,length(ls.nameTerms$mean)), variance = rep(1,length(ls.nameTerms$variance)))
    }else{
        out.glht <- try(glht(object, linfct = effects), silent = TRUE)
        if(inherits(out.glht,"try-error")){
            stop("Incorrect argument \'effects\': can be \"mean\", \"variance\", \"correlation\", \"all\", \n",
                 "or something compatible with the argument \'linfct\' of multcomp::glht. \n ")
        }
        out <- list(custom=NULL)
        ls.nameTerms <- list(custom = NULL)
        ls.nameTerms.num <- list(custom = 1)
        ls.contrast <- list(custom = out.glht$linfct)
        ls.null  <- list(custom = out.glht$rhs)
        
    }
    type.information <- attr(object$information,"type.information")
    

    ## ** prepare
    vcov.param <- vcov(object, effects = "all")
    dVcov.param <- attr(object$df,"dVcov")
    
    if(type.information != "observed"){
        warning("when using REML with expected information, the degree of freedom of the F-statistic may depend on the parametrisation of the variance parameters. \n")
    }

    ## ** F-tests
    for(iType in names(out)){

        ## skip empty type
        if(length(ls.nameTerms.num[[iType]])==0){ next }
        iLs <- lapply(ls.nameTerms.num[[iType]], function(iTerm){
            ## *** contrast matrix
            if(is.null(ls.contrast[[iType]])){
                iIndex.param <- which(ls.assign[[iType]]==iTerm)
                iN.hypo <- length(iIndex.param)
                iNull <- rep(ls.null[[iType]][iTerm],iN.hypo)
                iName.hypo <- paste(paste0(name.param[iIndex.param],"==",iNull), collapse = ", ")

                iC <- matrix(0, nrow = iN.hypo, ncol = p.param, dimnames = list(name.param[iIndex.param], name.param))
                if(length(iIndex.param)==1){
                    iC[name.param[iIndex.param],name.param[iIndex.param]] <- 1
                }else{
                    diag(iC[name.param[iIndex.param],name.param[iIndex.param]]) <- 1
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
            if(iType == "custom"){
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
