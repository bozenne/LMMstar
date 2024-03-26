### anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:38) 
## Version: 
## Last-Updated: Mar 24 2024 (21:19) 
##           By: Brice Ozenne
##     Update #: 1429
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * anova.lmm (documentation)
##' @title Multivariate Tests For Linear Mixed Model
##' @description Simultaneous tests of linear combinations of the model paramaters using Wald tests or Likelihood Ratio Test (LRT). 
##' 
##' @param object a \code{lmm} object. Only relevant for the anova function.
##' @param effects [character or numeric matrix] Should the Wald test be computed for all variables (\code{"all"}),
##' or only variables relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only variables relative to the variance structure (\code{"variance"}),
##' or only variables relative to the correlation structure (\code{"correlation"}).
##' or average counterfactual outcome with respect to the value of a covariate X at each repetition  (\code{"ACO_X"}),
##' or the contrast between average counterfactual outcome for any pair of value of a covariate X (\code{"ATE_X"}).
##' Can also be use to specify linear combinations of coefficients or a contrast matrix, similarly to the \code{linfct} argument of the \code{multcomp::glht} function.
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors. 
##' @param multivariate [logical] Should all hypotheses be simultaneously tested using a multivariate Wald test?
##' @param rhs [numeric vector] the right hand side of the hypothesis. Only used when the argument effects is a matrix.
##' @param df [logical] Should a F-distribution be used to model the distribution of the Wald statistic. Otherwise a chi-squared distribution is used.
##' @param ci [logical] Should an estimate, standard error, confidence interval, and p-value be output for each hypothesis?
##' @param transform.sigma,transform.k,transform.rho,transform.names are passed to the \code{vcov} method. See details section in \code{\link{coef.lmm}}.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A data.frame (LRT) or a list of containing the following elements (Wald):\itemize{
##' \item \code{multivariate}: data.frame containing the multivariate Wald test.
##' The column \code{df.num} refers to the degrees of freedom for the numerator (i.e. number of hypotheses)
##' wherease the column \code{df.denum} refers to the degrees of freedom for the denominator (i.e. Satterthwaite approximation).
##' \item \code{univariate}: data.frame containing each univariate Wald test.
##' \item \code{glht}: used internally to call functions from the multcomp package.
##' \item \code{object}: list containing key information about the linear mixed model.
##' \item \code{vcov}: variance-covariance matrix associated to each parameter of interest (i.e. hypothesis).
##' \item \code{iid}: matrix containing the influence function relative to each parameter of interest (i.e. hypothesis).
##' \item \code{args}: list containing argument values from the function call.
##' }
##' 
##' @details By default adjustment of confidence intervals and p-values for multiple comparisons is based on the distribution of the maximum-statistic.
##' This is refered to as a single-step Dunnett multiple testing procedures in table II of Dmitrienko et al. (2013).
##' It is performed using the multcomp package with the option \code{test = adjusted("single-step")} with equal degrees of freedom
##' or by simulation using a Student's t copula with unequal degrees of freedom (more in the note of the details section of \code{\link{confint.Wald_lmm}}).
##' 
##' @seealso
##' \code{\link{summary.Wald_lmm}} or \code{\link{confint.Wald_lmm}} for a summary of the results. \cr
##' \code{\link{autoplot.Wald_lmm}} for a graphical display of the results. \cr
##' \code{\link{rbind.Wald_lmm}} for combining result across models and adjust for multiple comparisons. \cr
##' 
##' @references Dmitrienko, A. and D'Agostino, R., Sr (2013), Traditional multiplicity adjustment methods in clinical trials. Statist. Med., 32: 5172-5218. https://doi.org/10.1002/sim.5990.
##'  
##' @keywords htest
##' 
##' @examples
##' #### simulate data in the long format ####
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' #### fit Linear Mixed Model ####
##' eUN.lmm <- lmm(Y ~ visit + X1 + X2 + X5,
##'                repetition = ~visit|id, structure = "UN", data = dL)
##' 
##' #### Multivariate Wald test ####
##' ## F-tests
##' anova(eUN.lmm)
##' anova(eUN.lmm, effects = "all")
##' anova(eUN.lmm, robust = TRUE, df = FALSE)
##' summary(anova(eUN.lmm), method = "bonferroni")
##' 
##' ## user defined F-test
##' summary(anova(eUN.lmm, effects = c("X1=0","X2+X5=10")))
##' 
##' ## chi2-tests
##' anova(eUN.lmm, df = FALSE)
##' 
##' ## with standard contrast
##' if(require(multcomp)){
##' amod <- lmm(breaks ~ tension, data = warpbreaks)
##' e.amod <- anova(amod, effect = mcp(tension = "Tukey"))
##' summary(e.amod)
##' }
##' 
##' #### Likelihood ratio test ####
##' eUN0.lmm <- lmm(Y ~ X1 + X2, repetition = ~visit|id, structure = "UN", data = dL)
##' anova(eUN.lmm, eUN0.lmm) 
##' 
##' eCS.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "CS", data = dL)
##' anova(eUN.lmm, eCS.lmm)

## * anova.lmm (code)
##' @export
anova.lmm <- function(object, effects = NULL, robust = FALSE, multivariate = TRUE, rhs = NULL, df = !is.null(object$df), ci = TRUE, 
                      transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

    call <- match.call()

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
        attr(out,"call") <- call
        if(!inherits(out,"LRT_lmm")){  ## in the case of re-estimation of the model via ML no need to re-assign the class
            class(out) <- append("LRT_lmm",class(out))
        }
    }else{ ## Wald test
        attr(robust, "call") <- "robust" %in% names(call)
        out <- .anova_Wald(object, effects = effects, robust = robust, multivariate = multivariate, rhs = rhs, df = df, ci = ci, 
                           transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        attr(out,"call") <- call
        class(out) <- append("Wald_lmm",class(out))
    }
    ## ** export
    return(out)
}

## * .anova_Wald
.anova_Wald <- function(object, effects, robust, multivariate, rhs, df, ci, 
                        transform.sigma, transform.k, transform.rho, transform.names){

    options <- LMMstar.options()
    type.information <- object$args$type.information
    original.transform.sigma <- transform.sigma
    original.transform.k <- transform.k
    original.transform.rho <- transform.rho
    
    ## ** normalized user input
    terms.mean <- attr(stats::terms(object$formula$mean.design),"term.labels")
    subeffect <- NULL
    
    init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho)
    
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho

    name.coef <- names(stats::coef(object, effects = "all"))
    out <- list(multivariate = NULL, univariate = NULL, glht = NULL)

    if(inherits(effects,"mcp")){
        out.glht <- try(multcomp::glht(object, linfct = effects), ## only used for generating contrast matrix
                        silent = TRUE)
        if(inherits(out.glht,"try-error")){
            stop("Possible mispecification of the argument \'effects\' as running mulcomp::glht lead to the following error: \n",
                 out.glht)
        }
        type <- "all"
        if(length(names(effects))==1){
            ls.nameTerms <- list(all = names(effects))
        }else{
            ls.nameTerms <- list(all = NULL)
        }
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
                stop("Argument \'effect\' should have column names matching the coefficient names when a matrix. \n")
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
        ls.nameTerms <- list(all = NULL)
        ls.nameTerms.num <- list(all = 1)
        ls.contrast <- list(all = out.glht$linfct)
        ls.null  <- list(all = out.glht$rhs)        
        name.effects <- rownames(effects)
        type <- "all"

    }else if(all(tolower(effects) %in% c("all","mean","fixed","variance","correlation"))){

        if("all" %in% effects){
            effects <- c("mean","variance","correlation")
        }
        if(transform.k %in% c("sd","var","logsd","logvar")){
            stop("Cannot use \'transform.rho\' equal \"sd\", \"var\", \"logsd\", or \"logvar\". \n",
                 "anova does not handle tests where the null hypothesis is at a boundary of the support of a random variable. \n")
        }
        if(transform.rho %in% c("cov")){
            stop("Cannot use \'transform.rho\' equal \"cov\". \n",
                 "anova does not handle tests where the null hypothesis is at a boundary of the support of a random variable. \n")
        }
        
        effects <- match.arg(effects, c("all","mean","fixed","variance","correlation"), several.ok = TRUE)
        effects[effects=="fixed"] <- "mean"
        name.effects <- NULL
        
        type <- NULL
        ls.assign <- list()
        ls.nameTerms <- list()
        ls.contrast <- list()
        ls.null <- list()
        if("mean" %in% effects){
            type <- c(type, "mean")
            ls.assign$mean <- attr(object$design$mean,"assign")
            ls.nameTerms$mean <- attr(stats::terms(object$formula$mean.design),"term.labels")
            ls.contrast <- c(ls.contrast,list(mean = NULL))
            null.mean <- 0
            ls.null$mean <- rep(null.mean,length(ls.nameTerms$mean))            
        }
        if("variance" %in% effects){
            type <- c(type, "variance")
            ls.assign$variance <- attr(object$design$vcov$var$X,"assign")            
            ls.nameTerms$variance <- attr(stats::terms(object$formula$var),"term.labels")
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
            type <- c(type, "correlation")
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
        
    }else{

        if(all(grepl("=",effects)==FALSE)){
            X <- object$design$mean
            if(all(effects %in% attr(X,"variable"))){
                assign <- which(attr(X,"variable") %in% effects)
                labels <- colnames(X)[attr(X,"assign") %in% assign]
                if(is.null(rhs)){
                    rhs <- rep(0, NROW(effects))
                }
                effects <- paste0(labels,"=",rhs)
            }else{
                stop("Incorrect argument \'effects\': can be \"mean\", \"variance\", \"correlation\", \"all\", \n",
                     "or an equation such compatible with the argument \'linfct\' of multcomp::glht \n ",
                     "or a contrast matrix. \n",
                     "or covariate names \n ")
            }
        }

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
            test.reparametrize <- any(grepl("log", c(object$reparametrize$transform.sigma,object$reparametrize$transform.k))) || grep("atanh", object$reparametrize$transform.rho)
            
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
        type <- "all"
        ls.contrast <- list(all = out.glht$linfct)
        ls.null  <- list(all = out.glht$rhs)
        name.effects <- names(effects)
        if(!is.null(name.effects)){
            rownames(ls.contrast$all) <- name.effects            
        }
        if(length(effects)==1){
            ls.nameTerms <- list(all = name.effects)
        }else{
            ls.nameTerms <- list(all = NULL)
        }
        ls.nameTerms.num <- list(all = 1)
    }

    ## ** prepare
    if(robust && object$args$method.fit=="REML"){
        name.mean <- names(coef(object, effects = "mean"))
        if(all(effects %in% c("mean","variance","correlation"))){
            if(all(effects %in% c("variance","correlation"))){
                stop("Cannot test variance, covariance, or correlation parameters using robust standard errors with REML. \n")
            }
        }else if(any(names(which(colSums(abs(do.call(rbind,ls.contrast)))>0)) %in% name.mean == FALSE)){
            stop("Cannot test variance, covariance, or correlation parameters using robust standard errors with REML. \n")
        }
        ls.contrast <- lapply(ls.contrast, function(iC){
            if(!is.null(iC)){iC[,name.mean,drop=FALSE]}else{iC}
        })
        effects <- "mean"        
    }else{
        effects <- "all"
    }
    type.param <- stats::setNames(object$design$param$type,object$design$param$name)
    param <- coef(object, effects = effects,
                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
    newparam <- coef(object, effects = effects,
                     transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
    newname <- names(newparam)
    name.param <- names(param)
    name.paramSigma <- names(type.param)[type.param=="sigma"]
    name.paramK <- names(type.param)[type.param=="k"]
    name.paramRho <- names(type.param)[type.param=="rho"]
    n.param <- length(param)
    vcov.param <- vcov(object, df = df*2, effects = effects, robust = robust,
                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
    dVcov.param <- attr(vcov.param,"dVcov")
    if(df>0 && object$args$method.fit=="REML" && type.information == "expected"){
        warning("when using REML with expected information, the degree of freedom of the F-statistic may depend on the parametrisation of the variance parameters. \n")
    }
    if(any(sapply(ls.contrast, function(iC){is.null(iC) || identical(colnames(iC), names(param))}) == FALSE)){
        warning("Names of the columns of the contrast matrix do not match the names of the model coefficients. \n")
    }

    ## ** F-tests
    out$glht <- stats::setNames(vector(mode = "list", length = length(type)), type)

    for(iType in type){ ## iType <- "variance"
        ## skip empty type
        if(length(ls.nameTerms.num[[iType]])==0 || (is.null(ls.contrast[[iType]]) && (all(ls.assign[[iType]]==0)))){ next }

        iParam <- coef(object, effects = iType,
                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
        name.iParam <- names(iParam)
        
        if(is.null(ls.nameTerms[[iType]])){
            ls.nameTerms[[iType]] <- ls.nameTerms.num[[iType]]
        }
        out$glht[[iType]] <- stats::setNames(vector(mode = "list", length = length(ls.nameTerms.num[[iType]])), ls.nameTerms[[iType]])

        for(iTerm in ls.nameTerms.num[[iType]]){ ## iTerm <- 1
            iNameTerm <- ls.nameTerms[[iType]][[which(iTerm == ls.nameTerms.num[[iType]])]]

            ## *** contrast matrix
            if(is.null(ls.contrast[[iType]])){
                if(all(ls.assign[[iType]]!=iTerm)){next}
                if(!is.null(subeffect) && ls.nameTerms$mean[iTerm]!=subeffect){next}
                iIndex.param <- which(ls.assign[[iType]]==iTerm)
                iN.hypo <- length(iIndex.param)
                iNull <- rep(ls.null[[iType]][iTerm],iN.hypo)
                iName.hypo <- paste(paste0(name.iParam[iIndex.param],"=",iNull), collapse = ", ")

                iC <- matrix(0, nrow = iN.hypo, ncol = n.param, dimnames = list(name.iParam[iIndex.param], name.param))
                if(length(iIndex.param)==1){
                    iC[name.iParam[iIndex.param],name.iParam[iIndex.param]] <- 1
                }else{
                    diag(iC[name.iParam[iIndex.param],name.iParam[iIndex.param]]) <- 1
                }
                iC.uni <- iC
            }else{
                iC <- ls.contrast[[iType]]
                iC.uni <- iC
                iN.hypo <- NROW(iC)
                iNull <- ls.null[[iType]]
                iName.hypo  <- paste(paste0(rownames(iC),"=",iNull),collapse=", ")
            }

            ## *** statistic
            if(multivariate){
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
                        iName.hypo <- paste(paste0(rownames(outSimp$C),"=",outSimp$rhs),collapse=", ")
                    }
                }
            }else{
                iStat <- NA
                iDf <- c(NA,NA)
            }

            ## *** degree of freedom
            if(multivariate && df>0 && !inherits(iC.vcov.C_M1,"try-error")){

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
            if(ci>0){
                if(df>0){
                    ci.df <-  .dfX(X.beta = iC.uni, vcov.param = vcov.param, dVcov.param = dVcov.param, return.vcov = ci>0.5)
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
                CI$statistic <- (CI$estimate-iNull)/CI$se
                rownames(CI) <- rownames(iC)
                if(!is.null(names(effects)) && !inherits(effects,"mcp")){                    
                    indexName <- intersect(which(names(effects)!=""),which(!is.na(names(effects))))
                    rownames(CI)[indexName] <- names(effects)[indexName]
                }
                CI.glht <- multcomp::glht(object, linfct = iC, rhs = iNull, df = ceiling(stats::median(ci.df)),
                                          coef. = function(iX){coef.lmm(iX, effects = effects)},
                                          vcov. = function(iX){vcov.lmm(iX, robust = robust, effects = effects)})
                CI.glht$model <- NULL
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
            ## R2 calculation from
            ## "An R2 statistic for fixed effects in the linear mixed model" by Lloyd J. Edwards et al. 2008 (Statistic in medicine)
            ## Equation 19
            ## DOI: 10.1002/sim.3429
            if(iType=="all"){
                iType2 <- apply(iC, MARGIN = 1, FUN = function(iRow){
                    ## type of parameter for the given contrast
                    iType <- unique(type.param[names(iRow)[which(abs(iRow)>0)]])
                    if(length(iType)==0){ ## if the contrast includes no parameter uses the type of parameter for the whole contrast matrix
                        iType <- unique(type.param[names(which(colSums(abs(iC)>0)>0))])
                    }
                    if(length(iType)==0){ ## use by effects.lmm when asking for the treatment effect on the change from baseline in a lmm without interaction: no coefficient to test.
                        return("all")
                    }else if(length(iType)>1){
                        return("all")
                    }else{
                        return(iType)
                    }
                })
            }else{
                iType2 <- switch(iType,
                                 "mean"="mu",
                                 "variance"="k",
                                 "correlation"="rho")
            }
            if(length(unique(iType2))==1){
                out$multivariate <- rbind(out$multivariate, cbind(type = unname(iType2[1]), test = iNameTerm, iRes))
            }else{
                out$multivariate <- rbind(out$multivariate, cbind(type = "all", test = iNameTerm, iRes))
            }
            if(ci){
                out$univariate <- rbind(out$univariate, cbind(type = iType2, test = iNameTerm, CI))
            }
            out$glht[[iType]][[iNameTerm]] <- CI.glht
        }
    }

    ## ** save some of the objects for possible use of rbind.Wald_lmm
    if(ci > 0.5){
        out$object <- list(outcome = object$outcome$var,
                           method.fit = object$args$method.fit,
                           type.information = type.information,
                           cluster.var = object$cluster$var,
                           structure = c(type = object$design$vcov$type, class = object$design$vcov$class))
        if(!is.na(attr(object$cluster$var,"original"))){
            out$object$cluster <- object$cluster$levels
        }
                           
        globalC <- do.call(rbind,lapply(unlist(out$glht, recursive=FALSE),"[[","linfct"))

        ## default: non-robust vcov and robust for iid
        if(attr(robust,"call")){
            out$vcov <- globalC %*% vcov.param %*% t(globalC)
            robust2 <- robust
        }else{
            out$vcov <- globalC %*% vcov(object, df = FALSE, effects = "all", robust = FALSE,
                                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE) %*% t(globalC)
            robust2 <- TRUE
        }
        if(object$args$method.fit == "ML"){
            out$iid <- iid(object, effects = effects, robust = robust2, 
                           transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names) %*% t(globalC)
        }else if(object$args$method.fit == "REML" && all(type.param[which(colSums(globalC!=0)>0)]=="mu")){
            globalC <- globalC[,name.param[type.param=="mu"],drop=FALSE]
            out$iid <- iid(object, effects = "mean", robust = robust2, 
                           transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names) %*% t(globalC)
        }
    }

    ## ** prepare for back-transformation
    out$args <- data.frame(type = NA, robust = robust, df = df, ci = ci,
                           transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                           transform.names = transform.names)

    if(all(name.paramSigma %in% colnames(globalC) == FALSE) || all(globalC[,name.paramSigma,drop=FALSE]==0)){
        out$args$transform.sigma <- NA
    }
    if(all(name.paramK %in% colnames(globalC) == FALSE) || all(globalC[,name.paramK,drop=FALSE]==0)){
        out$args$transform.k <- NA
    }
    if(all(name.paramRho %in% colnames(globalC) == FALSE) || all(globalC[,name.paramRho,drop=FALSE]==0)){
        out$args$transform.rho <- NA
    }

    test.original <- is.null(original.transform.sigma) && is.null(original.transform.k) && is.null(original.transform.rho)
    test.backtransform <- stats::na.omit(c(sigma = out$args$transform.sigma, k = out$args$transform.k, rho = out$args$transform.rho))
    ## should back transformation be performed when latter calling coef/confint?
    ## - no if no sigma/k/rho parameters in the linear hypotheses
    ## - no for the p-value, yes for the estimate, and NA for CI when evaluating a contrast on a tranformed scaled
    out$args$backtransform <- length(test.backtransform[test.backtransform != "none"])>0
    
    if(test.original && out$args$backtransform){ 

        ## find parameters subject to backtransformation
        param.with.backtransform <- object$design$param[object$design$param$type %in% names(test.backtransform)[test.backtransform!="none"],"name"]

        ## restrict contrast matrix to hypotheses with backtransformation
        M.allContrast <- do.call(rbind,lapply(out$glht[lengths(out$glht)>0], function(iLs){do.call(rbind,lapply(iLs,"[[","linfct"))}))
        col.with.backtransform <- intersect(colnames(M.allContrast),param.with.backtransform)
        row.with.backtransform <- rowSums(M.allContrast[,col.with.backtransform,drop=FALSE]!=0)>0
        Mback.allContrast <- M.allContrast[row.with.backtransform,,drop=FALSE]

        ## sanity check
        if(any(rowSums(M.allContrast[row.with.backtransform,setdiff(colnames(M.allContrast),param.with.backtransform),drop=FALSE]!=0)>0)){
            warning("Hypothesis with mean and variance-covariance parameter(s): possible misleading results due to transformation of the variance-covariance parameter(s). \n")
        }

        ## evaluate un-transformed value
        paramMback.allContrast <- rowSums(Mback.allContrast!=0)
        if(any(paramMback.allContrast>1)){
            param.untransformed <- coef(object, effects = effects, transform.sigma = "none", transform.k = "none", transform.rho = "none")
            attr(out$univariate,"backtransform") <- (Mback.allContrast[paramMback.allContrast>1,,drop=FALSE] %*% param.untransformed)[,1]
        }        
    }

    ## nice naming
    if(transform.names && is.null(names(effects))){
        backtransform.names <- names(coef(object,
                                          effects = effects, 
                                          transform.sigma = gsub("log","",transform.sigma),
                                          transform.k = gsub("log","",transform.k),
                                          transform.rho = gsub("atanh","",transform.rho),
                                          transform.names = transform.names))
        out$args$backtransform.names <- list(stats::setNames(rownames(out$univariate),rownames(out$univariate)))
        index <- match(newname,out$args$backtransform.names[[1]])
        out$args$backtransform.names[[1]][stats::na.omit(index)] <- backtransform.names[which(!is.na(index))]
    }
    if(identical(type,"all")){
        out$args$type <- list("all")
    }else{
        out$args$type <- list(c("mu"="mean", "k"="variance", "rho"="correlation")[out$multivariate$type])
    }

    ## ** export
    return(out)
}

## * .anova_LRT
.anova_LRT <- function(object1,object2){
    tol <- 1e-10
    
    ## ** normalize user input
    ## *** re-order models
    logLik1 <- logLik(object1)
    logLik2 <- logLik(object2)
    if(is.na(logLik1) || is.na(logLik2)){
        stop("Cannot perform a likelihood ratio test when the log-likelihood is NA for one of the models.\n")
    }else if(logLik2>=logLik1){
        type <- "2-1"
        objectH0 <- object1
        objectH1 <- object2
    }else if(logLik1>=logLik2){
        type <- "1-2"
        objectH0 <- object2
        objectH1 <- object1
    }

    ## *** extract coefficients
    name.paramH0 <- names(coef(objectH0, effects = "all"))   
    name.paramH1 <- names(coef(objectH1, effects = "all"))
    table.paramH0 <- objectH0$design$param
    table.paramH1 <- objectH1$design$param
    type.paramH0 <- stats::setNames(table.paramH0[match(name.paramH0, table.paramH0$name),"type"], name.paramH0)
    type.paramH1 <- stats::setNames(table.paramH1[match(name.paramH1, table.paramH1$name),"type"], name.paramH1)
    mismatchH0 <- stats::setNames(name.paramH0 %in% name.paramH1 == FALSE, name.paramH0)
    mismatchH1 <- stats::setNames(name.paramH1 %in% name.paramH0 == FALSE, name.paramH1)
    rhs <- stats::setNames(rep("0", sum(mismatchH1)), names(which(mismatchH1)))
    current.mismatchH0 <- mismatchH0
    current.mismatchH1 <- mismatchH1

    ## *** objective function
    if(object1$args$method.fit!=object2$args$method.fit){
        stop("The two models should use the same type of objective function for the likelihood ratio test to be valid. \n")
    }
    if(objectH1$args$method.fit=="REML" && (any(type.paramH1[mismatchH1]=="mu") || any(type.paramH0[mismatchH0]=="mu"))){
        objectH0$call$method.fit <- "ML"
        objectH1$call$method.fit <- "ML"
        message("Cannot use a likelihood ratio test to compare mean parameters when the objective function is REML. \n",
                "Will re-estimate the model via ML and re-run the likelihood ratio test. \n")
        out <- anova(eval(objectH0$call),eval(objectH1$call))
        attr(out,"type") <- type
        return(out)
    }

    ## *** number of observations
    nobsH0 <- nobs(objectH0)
    nobsH1 <- nobs(objectH1)
    if(any(nobsH0 != nobsH1)){
        if(nobsH0["missing"]!=nobsH1["missing"]){
            stop("Mismatch between the number of observations between the two models - could be due to missing data. \n",
                 "H0: ",paste(paste(names(nobsH0),"=",nobsH0), collapse = ", "),".\n",
                 "H1: ",paste(paste(names(nobsH1),"=",nobsH0), collapse = ", "),".\n")
        }else{
            stop("Mismatch between the number of observations between the two models. \n",
                 "H0: ",paste(paste(names(nobsH0),"=",nobsH0), collapse = ", "),".\n",
                 "H1: ",paste(paste(names(nobsH1),"=",nobsH0), collapse = ", "),".\n")
        }
    }

    ## *** check nesting
    if(any(abs(objectH1$design$Y-objectH0$design$Y)>tol)){
            stop("Mismatch in outcome between the two models. \n")
    }

    diff.strata <- setdiff(objectH0$design$vcov$name$strata,objectH0$design$vcov$name$strata)
    if(length(diff.strata)){
        stop("Cannot perform a likelihood ratio test when the model are not nested. \n",
             "Variance-covariance parameters are stratified with respect to \"",diff.strata,"\" in the model with the smallest likelihood but not in the model with the largest likelihood. \n",sep="")
    }

    if(any(mismatchH0)){
        ## name may not match even with nested model when using different covariance patterns
        ## e.g. CS: rho while UN gives rho(1,2), rho(1,3), ...
        response <- objectH1$design$Y
        
        if(length(name.paramH0)==length(name.paramH1)){
            stop("Cannot perform a likelihood ratio test when the model are not nested. \n",
                 "The two models have the same number of parameters. \n")
        }
        if(length(name.paramH0)>length(name.paramH1)){
            stop("Cannot perform a likelihood ratio test when the model are not nested. \n",
                 "The model with the largest likelihood had less parameters than the model with the smallest likelihood. \n")
        }
        
        ## mean
        if("mu" %in% type.paramH0[mismatchH0]){
            mismatchH0.mu <- names(type.paramH0[mismatchH0])[type.paramH0[mismatchH0] %in% c("mu")]
            mismatchH1.mu <- names(type.paramH1[mismatchH1])[type.paramH1[mismatchH1] %in% c("mu")]
            if(length(mismatchH0.mu)>length(mismatchH1.mu)){
                stop("Cannot perform a likelihood ratio test when the model are not nested. \n",
                     "The model with the largest likelihood has less mean parameters than the model with the smallest likelihood. \n")
            }

            XH0mu <- model.matrix(objectH0, effects = "mean")
            XH1mu <- model.matrix(objectH1, effects = "mean")
            H.X0mu <- XH0mu %*% solve(crossprod(XH0mu)) %*% t(XH0mu)
            H.X1mu <- XH1mu %*% solve(crossprod(XH1mu)) %*% t(XH1mu)
            
            ## attempt to check whether the design matrices are equivalent
            if(any(abs(H.X0mu-H.X1mu)>tol)){
                stop("Do not understand how the two models are nested. \n",
                     "The model with the lowest likelihood had mean parameters \"",paste(mismatchH0.mu, collapse="\" \""),"\" that were not present in the model with the largest likelihood. \n",
                     "Consider reparametrizing the mean structure so the nesting is more obvious.\n")
            }else{
                current.mismatchH0[mismatchH0.mu] <- FALSE
                current.mismatchH1[mismatchH1.mu] <- FALSE
                rhs <- rhs[names(current.mismatchH1)[current.mismatchH1]]
            }            
        }

        ## variance
        if("sigma" %in% type.paramH0[mismatchH0] || "k" %in% type.paramH0[mismatchH0]){
            mismatchH0.ksigma <- names(type.paramH0[mismatchH0])[type.paramH0[mismatchH0] %in% c("sigma","k")]
            mismatchH1.ksigma <- names(type.paramH1[mismatchH1])[type.paramH1[mismatchH1] %in% c("sigma","k")]
            if(length(mismatchH0.ksigma)>length(mismatchH1.ksigma)){
                stop("Cannot perform a likelihood ratio test when the model are not nested. \n",
                     "The model with the largest likelihood has less variance parameters than the model with the smallest likelihood. \n")
            }
             
            if(identical(objectH1$design$vcov$name$strata,objectH0$design$vcov$name$strata)){
                
                XH0ksigma <- model.matrix(objectH0, effects = "variance")$var$X
                XH1ksigma <- model.matrix(objectH1, effects = "variance")$var$X
                H.X0ksigma <- XH0ksigma %*% solve(crossprod(XH0ksigma)) %*% t(XH0ksigma)
                H.X1ksigma <- XH1ksigma %*% solve(crossprod(XH1ksigma)) %*% t(XH1ksigma)

                ## attempt to check whether the design matrices are equivalent
                if(any(abs(H.X0ksigma-H.X1ksigma)>tol)){
                    stop("Do not understand how the two models are nested. \n",
                         "The model with the lowest likelihood had variance parameters \"",paste(mismatchH0.ksigma, collapse="\" \""),"\" that were not present in the model with the largest likelihood. \n",
                         "Consider reparametrizing the mean structure so the nesting is more obvious.\n")
                }else{
                    current.mismatchH0[mismatchH0.ksigma] <- FALSE
                    current.mismatchH1[mismatchH1.ksigma] <- FALSE
                    rhs <- rhs[names(current.mismatchH1)[current.mismatchH1]]
                }
            }else if(all(is.na(objectH0$design$vcov$name$strata))){
                mismatchH1.ksigma2 <- sapply(strsplit(mismatchH1.ksigma,":",fixed = TRUE),"[[",1)
                if(all(mismatchH1.ksigma2 %in% mismatchH0.ksigma)){
                    current.mismatchH0[mismatchH0.ksigma] <- FALSE
                    rhs[mismatchH1.ksigma[mismatchH1.ksigma2 %in% mismatchH0.ksigma]] <- mismatchH1.ksigma2[mismatchH1.ksigma2 %in% mismatchH0.ksigma]
                }else{
                    stop("Do not understand how the two models are nested. \n",
                         "The model with the lowest likelihood had variance parameters \"",paste(mismatchH0.ksigma, collapse="\" \""),"\" that were not present in the model with the largest likelihood. \n",
                         "Consider reparametrizing the mean structure so the nesting is more obvious.\n")
                }
            }else{
                stop("Do not understand how the two models are nested. \n",
                     "The model with the lowest likelihood had variance parameters \"",paste(mismatchH0.ksigma, collapse="\" \""),"\" that were not present in the model with the largest likelihood. \n",
                     "Consider reparametrizing the mean structure so the nesting is more obvious.\n")
            }
        }

        ## correlation
        if("rho" %in% type.paramH0[mismatchH0]){
            mismatchH0.rho <- names(type.paramH0[mismatchH0])[type.paramH0[mismatchH0] %in% c("rho")]
            mismatchH1.rho <- names(type.paramH1[mismatchH1])[type.paramH1[mismatchH1] %in% c("rho")]
            if(length(mismatchH0.rho)>length(mismatchH1.rho)){
                stop("Cannot perform a likelihood ratio test when the model are not nested. \n",
                     "The model with the largest likelihood has less correlation parameters than the model with the smallest likelihood. \n")
            }

            if(identical(objectH1$design$vcov$name$strata,objectH0$design$vcov$name$strata)){
                ## no strata or same strata: effect of the type of  structure
                n.strata <- objectH1$strata$n
                test.nested1 <- objectH1$design$vcov$class=="UN" && objectH0$design$vcov$class=="CS" && is.na(objectH0$design$vcov$name$var) && is.na(objectH0$design$vcov$name$cor)
                if(test.nested1){
                    mismatchH0 <- setdiff(mismatchH0,mismatchH0.rho)                
                    for(iS in 1:n.strata){
                        iRhoH0 <- table.paramH0[which((table.paramH0$type=="rho")*(table.paramH0$index.strata==1)==1),"name"]
                        iRhoH1 <- table.paramH1[which((table.paramH1$type=="rho")*(table.paramH1$index.strata==1)==1),"name"]
                        rhs[iRhoH1] <- iRhoH0
                    }
                }else{
                    stop("Do not understand how the two models are nested. \n",
                         "The model with the lowest likelihood had correlation parameters \"",paste(mismatchH0.rho, collapse="\" \""),"\" that were not present in the model with the largest likelihood. \n")
                }
            }else{
                ## effect of strata
                mismatchH1.rho2 <- sapply(strsplit(mismatchH1.rho,":",fixed = TRUE),"[[",1)
                if(all(mismatchH1.rho2 %in% mismatchH0.rho)){
                    current.mismatchH0[mismatchH0.rho] <- FALSE
                    rhs[mismatchH1.rho[mismatchH1.rho2 %in% mismatchH0.rho]] <- mismatchH1.rho2[mismatchH1.rho2 %in% mismatchH0.rho]
                }else{
                    stop("Do not understand how the two models are nested. \n",
                         "The model with the lowest likelihood had correlation parameters \"",paste(mismatchH0.rho, collapse="\" \""),"\" that were not present in the model with the largest likelihood. \n")
                }
            }
        }
    }

    n.paramTest <- length(name.paramH1)-length(name.paramH0)

    ## ** LRT
    out <- data.frame(null = paste(paste0(names(rhs),"==",rhs), collapse = ", "),
                      logLikH1 = stats::logLik(objectH1),
                      logLikH0 = stats::logLik(objectH0),
                      statistic = NA,
                      df = n.paramTest,
                      p.value = NA,
                      stringsAsFactors = FALSE)
    out$statistic <- 2*(out$logLikH1 - out$logLikH0)
    out$p.value <- 1 - stats::pchisq(out$statistic, df = out$df)

    ## ** export
    attr(out,"type") <- type
    return(out)
}

## * anova.mlmm
##' @export
anova.mlmm <- function(object, effects = NULL, rhs = NULL, ...){

    ## ** normalize argument
    if(is.null(effects)){
        class(object) <- setdiff(class(object), c("mlmm"))
        return(object)
    }
    if(inherits(effects, "mcp")){
        if(length(effects)!=1){
            stop("Argument \'effects\' must specify a single hypothesis test when being of class \"mcp\". \n",
                 "Something like mcp(group = \"Dunnett\") or mcp(group = \"Tukey\") \n")
        }
        effects.save <- effects
        constraint <- effects.save[[1]]
        effects <- names(effects.save)
        if(!grepl("=",effects)){
            effects <- paste0(effects,"=0")
        }
    }else{
        constraint <- NULL
    }

    ## ** test linear combinations
    robust <- object$args$robust
    df <- object$args$df 
    ci <- object$args$ci

    transform.sigma <- if(is.na(object$args$transform.sigma)){NULL}else{object$args$transform.sigma}
    transform.k <- if(is.na(object$args$transform.k)){NULL}else{object$args$transform.k}
    transform.rho <- if(is.na(object$args$transform.rho)){NULL}else{object$args$transform.rho}

    ls.lmm <- object$model
    name.lmm <- names(ls.lmm)
    ls.anova <- stats::setNames(lapply(name.lmm, function(iName){ ## iName <- name.lmm[1]
        anova(ls.lmm[[iName]], effects = effects, rhs = rhs, df = df, ci = ci, robust = robust,
              transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
    }), name.lmm)

    ## ** regenerate a new mlmm object
    out <- do.call("rbind.Wald_lmm",
                   args = c(list(model = ls.anova[[1]], effects = constraint, rhs = rhs, name = names(object$model), sep = object$args$sep), unname(ls.anova[-1]))
                   )
    
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
    if(length(rhs)==1 && NROW(object)>1){
        rhs <- rep(rhs, NROW(object))
    }

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
