### anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:38) 
## Version: 
## Last-Updated: maj 31 2022 (10:12) 
##           By: Brice Ozenne
##     Update #: 851
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
##' @description Simultaneous tests of linear combinations of the model paramaters using Wald tests. 
##' @name anova
##' 
##' @param object a \code{lmm} object. Only relevant for the anova function.
##' @param effects [character] Should the Wald test be computed for all variables (\code{"all"}),
##' or only variables relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only variables relative to the variance structure (\code{"variance"}),
##' or only variables relative to the correlation structure (\code{"correlation"}).
##' Can also be use to specify linear combinations of coefficients, similarly to the \code{linfct} argument of the \code{multcomp::glht} function.
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors. 
##' @param rhs [numeric vector] the right hand side of the hypothesis. Only used when the argument effects is a matrix.
##' @param df [logical] Should a F-distribution be used to model the distribution of the Wald statistic. Otherwise a chi-squared distribution is used.
##' @param ci [logical] Should a confidence interval be output for each hypothesis?
##' @param transform.sigma,transform.k,transform.rho,transform.names are passed to the \code{vcov} method. See details section in \code{\link{coef.lmm}}.
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
##' @seealso
##' \code{\link{summary.anova_lmm}} for a summary of the results. \cr
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
##' summary(anova(eUN.lmm, df = FALSE))
##' 
##' ## F-test
##' anova(eUN.lmm)
##' summary(anova(eUN.lmm, effects = "all"))
##' anova(eUN.lmm, effects = c("X1=0","X2+X5=10"))
##'
##' ## another example
##' if(require(multcomp)){
##' amod <- lmm(breaks ~ tension, data = warpbreaks)
##' e.glht <- glht(amod, linfct = mcp(tension = "Tukey"))
##' summary(e.glht, test = Chisqtest()) ## 0.000742
##'
##' e.amod <- anova(amod, effect = mcp(tension = "Tukey"))
##' summary(e.amod)
##' 
##' }

## * anova.lmm (code)
##' @rdname anova
##' @export
anova.lmm <- function(object, effects = NULL, robust = FALSE, rhs = NULL, df = !is.null(object$df), ci = TRUE, 
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
    attr(out,"call") <- match.call()
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
        }else if(grepl("^cor_",effects) && !is.null(object$design$vcov$X$cor.pairwise)){
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
        name.effects <- NULL
       
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
                                       vcov. = function(iX){vcov.lmm(iX, robust = FALSE, effects = "all")}),
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
        name.effects <- rownames(effects)

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
        name.effects <- NULL
        
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
            ls.assign$correlation <- rep(1,sum(object$design$param$type=="rho"))
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
                                       vcov. = function(iX){vcov.lmm(iX, robust = FALSE, effects = "all")}),
                        silent = TRUE)
        newname.coef <- names(stats::coef(object, effects = "all"))
        
        if(inherits(out.glht,"try-error")){

            ## maybe the error is due to white space in the name of the coefficients
            if(any(grepl(" ",effects))){
                lhs <- sapply(strsplit(effects, split = "=",fixed = TRUE),"[",1)
                lhs.term <- unlist(strsplit(unlist(strsplit(lhs, split = "+", fixed = TRUE)), split = "-", fixed = TRUE))
                lhs.variable <- trimws(lapply(strsplit(lhs.term, split = "*", fixed = TRUE), utils::tail, 1), which = "both")
                if(any(grepl(" ",lhs.variable))){
                    stop("Possible mispecification of the argument \'effects\' as running mulcomp::glht lead to the following error: \n",
                         out.glht,
                         "There seems to be whitespace(s) in the name of certain coefficients which is likely to be the cause of the error. \n",
                         "Consider using a contrast matrix to specific the argument \'effects\' \n or renaming the values of categorical variables without whitespace and refiting the lmm object. \n"
                         )
                }
            }
            
            test.reparametrize <- grepl("log", c(object$reparametrize$transform.sigma,object$reparametrize$transform.k)) || grep("atanh", object$reparametrize$transform.rho)
            
            ## restaure untransformed parametrization (glht does not handle log(k). or atanh(cor))
            if(test.reparametrize){
                object2 <- object
                index.var <- which(object$design$param$type %in% c("sigma","k","rho"))
                object2$reparametrize <- .reparametrize(p = object$param[index.var],
                                                        type = object$design$param$type[index.var], 
                                                        level = object$design$param$level[index.var], 
                                                        sigma = object$design$param$sigma[index.var], 
                                                        k.x = object$design$param$k.x[index.var], 
                                                        k.y = object$design$param$k.y[index.var], 
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
        name.effects <- names(effects)
        if(!is.null(name.effects)){
            rownames(ls.contrast$all) <- name.effects
        }
    }
    type.information <- attr(object$information,"type.information")    

    ## ** prepare
    if(robust && object$method.fit=="REML"){
        name.mean <- names(coef(object, effects = "mean"))
        if(any(names(which(colSums(abs(do.call(rbind,ls.contrast)))>0)) %in% name.mean == FALSE)){
            stop("Cannot test variance, covariance, or correlation parameters using robust standard errors with REML. \n")
        }
        effects <- "mean"
        ls.contrast <- lapply(ls.contrast, function(iC){iC[,name.mean,drop=FALSE]})
    }else{
        effects <- "all"
    }
    param <- coef(object, effects = effects,
                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
    newname <- names(coef(object, effects = effects,
                          transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names))
    name.param <- names(param)
    n.param <- length(param)
    vcov.param <- vcov(object, df = df*2, effects = effects, robust = robust,
                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
    dVcov.param <- attr(vcov.param,"dVcov")
    if(df>0 && object$method.fit=="REML" && type.information == "expected"){
        warning("when using REML with expected information, the degree of freedom of the F-statistic may depend on the parametrisation of the variance parameters. \n")
    }
    if(any(sapply(ls.contrast, function(iC){is.null(iC) || identical(colnames(iC), names(param))}) == FALSE)){
        warning("Names of the columns of the contrast matrix do not match the names of the model coefficients. \n")
    }

    ## ** F-tests
    type <- names(out)
    for(iType in type){ ## iType <- "correlation"
        ## skip empty type
        if(length(ls.nameTerms.num[[iType]])==0 || (is.null(ls.contrast[[iType]]) && (all(ls.assign[[iType]]==0)))){ next }

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
                iDf <- c(outSimp$dim,Inf)
                attr(iStat,"error") <- "\n  Could not invert the covariance matrix for the proposed contrast."
            }else{
                iStat <- as.double(t(outSimp$C %*% param - outSimp$rhs) %*% iC.vcov.C_M1 %*% (outSimp$C %*% param - outSimp$rhs))/outSimp$dim 
                iDf <- c(outSimp$dim,Inf)
                if(outSimp$rm>0){
                    iName.hypo <- paste(paste0(rownames(outSimp$C),"==",outSimp$rhs),collapse=", ")
                }
            }

            ## *** degree of freedom
            if(df>0 && !inherits(iC.vcov.C_M1,"try-error")){

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
                if(df>0){
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
                                 partial.R = NA,
                                 p.value = NA,
                                 stringsAsFactors = FALSE)
                CI$statistic <- (CI$estimate-iNull)/CI$se
                rownames(CI) <- rownames(iC)
                CI$partial.R <- sign(CI$statistic)*sqrt(CI$statistic^2 / (CI$df + CI$statistic^2))
                if(!is.null(names(effects)) && !inherits(effects,"mcp")){
                    indexName <- intersect(which(names(effects)!=""),which(!is.na(names(effects))))
                    rownames(CI)[indexName] <- names(effects)[indexName]
                }
                CI.glht <- multcomp::glht(object, linfct = iC, rhs = iNull, df = ceiling(stats::median(ci.df)),
                                          coef. = function(iX){coef.lmm(iX, effects = effects)},
                                          vcov. = function(iX){vcov.lmm(iX, robust = robust, effects = effects)})
                attr(CI.glht$vcov,"robust") <- robust
            }else{
                CI <- NULL
                CI.glht <- NULL
            }

            ## *** test
            iRes <- data.frame("null" = if(!is.null(name.effects)){paste(name.effects,collapse=", ")}else{iName.hypo},
                               "statistic" = iStat,
                               "df.num" = iDf[1],
                               "df.denom" = iDf[2],
                               "partial.R2"  =  iDf[1] * iStat / (iDf[2] + iDf[1] * iStat),
                               "p.value" = 1 - stats::pf(iStat, df1 = iDf[1], df2 = iDf[2]),
                               stringsAsFactors = FALSE)
            ## R2 calculation from
            ## "An R2 statistic for fixed effects in the linear mixed model" by Lloyd J. Edwards et al. 2008 (Statistic in medicine)
            ## Equation 19
            ## DOI: 10.1002/sim.3429

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
        attr(out[[iType]], "CI") <- lapply(iLs,attr,"CI")
        attr(out[[iType]], "glht") <- lapply(iLs,attr,"glht")
        
    }
    
    ## ** export
    attr(out, "df") <- df
    attr(out, "test") <- "Wald"
    attr(out, "robust") <- robust
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
    paramH1 <-  objectH1$design$param$name
    typeH1 <-  objectH1$design$param$type
    paramH0 <-  objectH0$design$param$name
    typeH0 <-  objectH0$design$param$type
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
