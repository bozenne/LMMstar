### anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:38) 
## Version: 
## Last-Updated: Jul 22 2024 (20:22) 
##           By: Brice Ozenne
##     Update #: 1676
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
##' Can also be use to specify linear combinations of coefficients or a contrast matrix, similarly to the \code{linfct} argument of the \code{multcomp::glht} function.
##' @param rhs [numeric vector] the right hand side of the hypothesis. Only used when the argument \code{effects} is a matrix.
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors. 
##' @param df [logical] Should degrees of freedom be estimated using a Satterthwaite approximation?
##' If yes F-distribution (multivariate) and Student's t-distribution (univariate) are used.
##' Other chi-squared distribution and normal distribution are used.
##' @param univariate [logical] Should an estimate, standard error, confidence interval, and p-value be output for each hypothesis?
##' @param multivariate [logical] Should all hypotheses be simultaneously tested using a multivariate Wald test?
##' @param transform.sigma,transform.k,transform.rho are passed to the \code{vcov} method. See details section in \code{\link{coef.lmm}}.
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
anova.lmm <- function(object, effects = NULL, rhs = NULL, robust = FALSE, df = !is.null(object$df),
                      univariate = TRUE, multivariate = TRUE, 
                      transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, ...){

    call <- match.call()
    options <- LMMstar.options()

    ## ** normalized user input    
    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** effects
    if(is.null(effects)){
        effects <- options$effects
    }


    ## ** run test
    if(inherits(effects,"lmm")){
        ## *** Likelihood Ratio Test (LRT)
        out <- .anova_LRT(object1 = object, object2 = effects)
    }else{  
        ## *** Wald test

        #n# extract from object
        object.coef <- stats::model.tables(object, effects = "param")
        name.coef <- object.coef$name
        n.coef <- length(name.coef)
        type.coef <- stats::setNames(object.coef$type, name.coef)

        ## initialize tranformation
        init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                                x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho,
                                simplify = FALSE)
        transform.sigma <- init$transform.sigma
        transform.k <- init$transform.k
        transform.rho <- init$transform.rho
        
        ## prepare output
        ls.contrast <- list()
        ls.null <- list()
        backtransform <- TRUE

        ## generate contrast matrix
        if(all(tolower(effects) %in% c("all","mean","fixed","variance","correlation"))){

            ## further extract from object
            object.X <- model.matrix(object, effects = "all", simplify = 0.5)

            ## further normalize input
            effects <- tolower(effects)
            if("all" %in% effects){
                if(length(effects) == 1){
                    effects <- c("mean","variance","correlation")
                }else{
                    stop("When argument \'effect\' contains \"all\" it should be of length 1. \n")
                }
            }else if("fixed" %in% effects){
                effects[effects=="fixed"] <- "mean"
            }
            if(init$transform.k %in% c("sd","var","logsd","logvar")){
                stop("Cannot use \'transform.rho\' equal \"sd\", \"var\", \"logsd\", or \"logvar\". \n",
                     "anova does not handle tests where the null hypothesis is at a boundary of the support of a random variable. \n")
            }
            if(all(attr(object.X$mean,"assign")>0)){
                effects <- setdiff(effects,"mean")
            }
            if("k" %in% type.coef == FALSE){
                effects <- setdiff(effects,"variance")
            }
            if("rho" %in% type.coef == FALSE){
                effects <- setdiff(effects,"correlation")
            }

            ## name of the sigma coefficient
            if("variance" %in% effects || "correlation" %in% effects){
                name.coef.sigma <- name.coef[type.coef == "sigma"]
                name.strata.sigma <- names(attr(object.X$var,"ls.level"))[match(name.coef.sigma,colnames(object.X$vcov$var$X))]

            if("mean" %in% effects & any(attr(object.X$mean,"assign")>0)){
           
                ## names of the terms in the design matrix
                terms.mean <- attr(object.X$mean,"term.labels")

                ## contrast matrix
                ls.contrast$mu <- stats::setNames(lapply(1:length(terms.mean), function(iT){ ## iT <- 1
                    iIndex <- which(attr(object.X$mean,"assign")==iT)
                    iCoef <- colnames(object.X$mean)[iIndex]
                    iC <- matrix(0, nrow = length(iIndex), ncol = n.coef,
                                 dimnames = list(paste0(iCoef,"=0"), name.coef))
                    iC[,iCoef] <- diag(1, nrow = length(iIndex))
                    return(iC)                
                }), terms.mean)
                ls.null$mu <- lapply(ls.contrast$mu, function(iC){stats::setNames(rep(0, NROW(iC)),rownames(iC))})
            }
            if("variance" %in% effects & any("k" %in% type.coef)){
                ## name of the sigma coefficient
                name.coef.sigma <- name.coef[type.coef == "sigma"]

                ## names of the k coefficients
                name.coef.k <- name.coef[type.coef == "k"]
                n.coef.k <- length(name.coef.k)

                ## terms
                if(object$strata$n==1){
                    terms.var <- unique(attr(object.X$var,"term.labels") [match(name.coef.k, colnames(object.X$var))])
                }else{
                    terms.var <- paste(names(attr(object.X$var,"ls.level")[match(name.coef.sigma, colnames(object.X$var))]),
                                       paste(attr(object.X$var,"variable")[-1],collapse = ":"), sep = ":")
                }

                ## null hypothesis
                null.variance <- switch(init$transform.k,
                                        "none" = 1,
                                        "square" = 1,
                                        "log" = 0,
                                        "logsquare" = 0)

                ## contrast matrix 
                ls.contrast$k <- stats::setNames(lapply(name.coef[type.coef=="sigma"], function(iSigma){ ## iSigma <- name.coef[type.coef=="sigma"][1]
                    iCoef <- intersect(name.coef.k,name.coef[attr(object.coef,"sigma") == iSigma])
                    iC <- matrix(0, nrow = length(iCoef), ncol = n.coef,
                                 dimnames = list(paste0(iCoef,"=0"), name.coef))
                    iC[,iCoef] <- diag(1, nrow = length(iCoef))
                    return(iC)                
                }),terms.var)
                ls.null$k <- lapply(ls.contrast$k, function(iC){stats::setNames(rep(null.variance, NROW(iC)),rownames(iC))})
            }
            if("correlation" %in% effects & NCOL(object.X$cor)>0){

                
                ## names of the rho coefficients
                name.coef.rho <- name.coef[type.coef == "rho"]
                n.coef.rho <- length(name.coef.rho)
                terms.cor <-  paste(name.strata.sigma,paste(object.X$vcov$name$cor[[1]],collapse=":"),sep=":")

                ## null hypothesis
                null.correlation <- switch(init$transform.rho,
                                           "none" = 0,
                                           "atanh" = 0)

                ## contrast matrix 
                ls.contrast$rho <- stats::setNames(lapply(name.coef[type.coef=="sigma"], function(iSigma){ ## iSigma <- name.coef[type.coef=="sigma"][1]
                    iCoef <- intersect(name.coef.rho,name.coef[attr(object.coef,"sigma") == iSigma])
                    iC <- matrix(0, nrow = length(iCoef), ncol = n.coef,
                                 dimnames = list(paste0(iCoef,"=0"), name.coef))
                    iC[,iCoef] <- diag(1, nrow = length(iCoef))
                    return(iC)                
                }),terms.cor)
                ls.null$rho <- lapply(ls.contrast$rho, function(iC){stats::setNames(rep(null.correlation, NROW(iC)),rownames(iC))})
            }
        
        }else if(inherits(effects,"mcp")){
            browser()
            out.glht <- try(multcomp::glht(object, linfct = effects), silent = TRUE)
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
            simplify <- FALSE ## keep vcov and iid
       
        }else if(is.matrix(effects)){
            browser()
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
            simplify <- FALSE ## keep vcov and iid
    if(any(sapply(ls.contrast, function(iC){is.null(iC) || identical(colnames(iC), names(param))}) == FALSE)){
        warning("Names of the columns of the contrast matrix do not match the names of the model coefficients. \n")
    }

    }else if(is.character(effects)){
        ## normalized user input (transform)
        original.transform.sigma <- transform.sigma
        original.transform.k <- transform.k
        original.transform.rho <- transform.rho
        init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                                x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho)
        transform.sigma <- init$transform.sigma
        transform.k <- init$transform.k
        transform.rho <- init$transform.rho

        object.X <- model.matrix(object, effects = "all", simplify = FALSE)
        object.coef <- names(coef(object, effects = "all", transform.names = FALSE))
        attr(object.coef, "rescue") <- names(coef(object, effects = "all",
                                                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = TRUE))
        
        out.eq2c <- equation2contrast(effects,
                                      name.coef = object.coef,
                                      X = list(mean = object.X$mean, var = object.X$vcov$var$X),
                                      name.arg = "effects")

        ls.contrast <- list(user = list("1" = out.eq2c$contrast))
        ls.null  <- list(user = list("1" = out.eq2c$rhs))
        if(out.eq2c$rescue==TRUE && (!is.null(original.transform.sigma) || !is.null(original.transform.k) || !is.null(original.transform.rho))){
            ## the user specifically requests the transformed scale
            colnames(ls.contrast$user) <- object.coef
            backtransform <- FALSE 
        }

    }else{
        stop("Incorrect argument 'effects': can be \"mean\", \"variance\", \"correlation\", \"all\", \n", 
             "or an equation such compatible with the argument 'linfct' of multcomp::glht \n ", 
             "or a contrast matrix. \n", "or covariate names \n ")
    }

        attr(robust, "call") <- "robust" %in% names(call)
        out <- .anova_Wald(object, contrast = ls.contrast, null = ls.null, robust = robust, df = df,
                           multivariate = multivariate, univariate = univariate, 
                           transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, backtransform = backtransform)
    }

    ## ** export
    attr(out,"call") <- call
    return(out)
}

## * .anova_Wald
.anova_Wald <- function(object, contrast, null, robust, df,
                        univariate, multivariate, 
                        transform.sigma, transform.k, transform.rho, backtransform){

    ## ** prepare
    ## *** all parameters with their type
    theta <- coef(object, effects = "all", transform.names = FALSE, simplify = FALSE)
    name.theta <- names(theta)
    type.theta <- stats::setNames(attr(theta,"type"),name.theta)
    if(transform.k %in% c("sd","logsd","var","logvar")){
        type.theta[type.theta=="k"] <- "sigma"
    }

    ## *** transformation (ignore sd->var or cor->cov)
    transform2.sigma <- switch(transform.sigma,
                              "log" = "log",
                              "logsquare" = "log",
                              "one" = "one",
                              "none")
    transform2.k <- switch(transform.k,
                          "log" = "log",
                          "logsquare" = "log",
                          "logsd" = "log",
                          "logvar" = "log",
                          "none")
    transform2.rho <- switch(transform.rho,
                             "atanh" = "atanh",
                             "none")

    ## *** gather all (multivariate) tests
    ls.grid <- lapply(names(contrast), function(iName){ ## iName <- names(contrast)[1]
        iDf <- data.frame(type.original = iName, type = NA, term = names(contrast[[iName]]), n.test = sapply(contrast[[iName]],NROW))
        if(iName=="user"){
            attr(iDf,"ls.type") <- apply(contrast[[iName]][[1]], MARGIN = 1, FUN = function(iRow){type.theta[names(which(iRow!=0))]}, simplify = FALSE)
            iUtype <- unique(unlist(attr(iDf,"ls.type")))
            if(length(iUtype)==1){
                iDf$type <- iUtype
            }else{
                iDf$type <- "all"
            }
        }else{
            iDf$type <- iName
        }
        return(iDf)
    })
    grid <- do.call(rbind,ls.grid)
    if(NROW(grid)>1){
        rownames(grid) <- paste(grid$type, grid$term, sep = "_")
    }
    n.grid <- NROW(grid)

    if(length(unique(grid$type))==1){
        effects <- switch(grid$type[1],
                          "mu" = "mean",
                          "k" = "variance",
                          "rho" = "correlation",
                          "all" = "all")
    }else{
        effects <- "all"
    }

    ## *** extract model coefficient and uncertainty
    param <- coef(object, effects = effects,
                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
    vcov.param <- vcov(object, df = df*2, effects = effects, robust = robust,
                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
    dVcov.param <- attr(vcov.param,"dVcov")

    ## *** output
    out <- list(multivariate = NULL,
                univariate = NULL,
                glht = NULL,
                object = NULL,
                vcov = NULL,
                object = list(outcome = object$outcome$var,
                              method.fit = object$args$method.fit,
                              type.information = object$args$type.information,
                              cluster.var = object$cluster$var,
                              cluster = object$cluster$levels),
                args = data.frame(type = effects, robust = robust, df = df,
                                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
                )

    if(multivariate){
        out$multivariate <- data.frame(matrix(NA, nrow = n.grid, ncol = 8,
                                              dimnames = list(NULL,c("type","term","null","statistic","df.num","df.denom","p.value","message"))))
        out$multivariate$type <- grid$type
        out$multivariate$test <- grid$term
        out$multivariate$null <- unname(sapply(unlist(contrast, recursive = FALSE), FUN = function(iRow){paste(rownames(iRow), collapse = ", ")}))
        rownames(out$multivariate) <- rownames(grid)
    }
    if(univariate){
        ## name will differ from term when testing mean, variance, and correlation parameters which may all refer to the same variable, say time
        out$univariate <- data.frame(matrix(NA, nrow = sum(grid$n.test), ncol = 13,
                                            dimnames = list(NULL,c("type","term","name","estimate","se","df","statistic","lower","upper","null","p.value","transform","backtransform"))))
        
        out$univariate$type <- unname(unlist(mapply(iType = grid$type, iN = grid$n.test, FUN = function(iType,iN){rep(iType,iN)}, SIMPLIFY = FALSE)))
        out$univariate$term <- unname(unlist(mapply(iTerm = grid$term, iN = grid$n.test, FUN = function(iTerm,iN){rep(iTerm,iN)}, SIMPLIFY = FALSE)))
        out$univariate$name <- rownames(grid)[merge(x = out$univariate, y = cbind(index = 1:n.grid, grid), by = c("type","term"), sort = FALSE)$index]
        rownames(out$univariate) <- unname(unlist(sapply(unlist(contrast, recursive = FALSE), rownames)))
        
        ## back-transformation
        if(all(grid$type == "mu") || (transform2.sigma=="none" && transform2.k=="none" && transform2.rho=="none")){ ## also include names(contrast) is "mean"
            ## no need for transformation for all linear combinations
            out$univariate[out$univariate$type %in% c("sigma","k","rho"),"transform"] <- FALSE
            out$univariate[out$univariate$type %in% c("sigma","k","rho"),"backtransform"] <- FALSE
        }else if("variance" %in% names(contrast) || "correlation" %in% names(effects)){
            ## potential log-transformation/exp-backtransform for linear combinations including sigma parameters
            out$univariate[out$univariate$type %in% "sigma","transform"] <- (transform2.sigma=="log")
            out$univariate[out$univariate$type %in% "sigma","backtransform"] <- (transform2.sigma=="log")
            ## potential log-transformation/exp-backtransform for linear combinations including k parameters
            out$univariate[out$univariate$type %in% "k","transform"] <- (transform2.k=="log")
            out$univariate[out$univariate$type %in% "k","backtransform"] <- (transform2.k=="log")
            ## potential atanh-transformation/tanh-backtransform for linear combinations including rho parameters
            out$univariate[out$univariate$type %in% "rho","transform"] <- (transform2.rho=="atanh")
            out$univariate[out$univariate$type %in% "rho","backtransform"] <- (transform2.rho=="atanh")
        }else{ ## only remains names(contrast) is "user"
            dfType <-  do.call(rbind,lapply(attr(ls.grid[[1]],"ls.type"), function(iType){
                if(length(unique(iType))>1){
                    return(data.frame(type = "all", p = length(iType)))
                }else{
                    return(data.frame(type = iType[1], p = length(iType)))
                }
            }))
            ## potential log-transformation/exp-backtransform for linear combinations including ONLY sigma parameters or ONLY k parameters
            ## NOTE: log(k2) - log(k1) --> log(k2/k1) --> k2/k1 after back-transformation
            out$univariate[dfType$type == "sigma","transform"] <- (transform2.sigma=="log")
            out$univariate[dfType$type == "sigma","backtransform"] <- backtransform && (transform2.sigma=="log")
            out$univariate[dfType$type == "k","transform"] <- (transform2.k=="log")
            out$univariate[dfType$type == "k","backtransform"] <- backtransform && (transform2.k=="log")
            ## potential atanh-transformation/tanh-backtransform for linear combinations including ONLY rho parameters
            if(any("rho" %in% out$univariate$type)){
                ## single rho parameter per linear combination
                out$univariate[dfType$type == "rho" & dfType$p == 1,"transform"] <- (transform2.rho=="atanh")
                out$univariate[dfType$type == "rho" & dfType$p == 1,"backtransform"] <- backtransform && (transform2.rho=="atanh")
                ## multiple rho parameter per linear combination: only output p-value on the transformed scale
                out$univariate[dfType$type == "rho" & dfType$p > 1,"transform"] <- FALSE
                out$univariate[dfType$type == "rho" & dfType$p > 1,"backtransform"] <- FALSE
            }
            if(any(out$univariate=="all")){
                if(backtransform){
                    stop("Results may be misleading with contrasts mixing parameters using distinct transformations. \n",
                         "Consider setting the arguments \'transform.sigma\', \'transform.k\', and \'transform.rho\' to \"none\". \n",
                         "Or expression the argument \'effects\' as a function of the transformed parameters. \n")
                }else{
                    out$univariate[dfType$type == "all","transform"] <- TRUE
                    out$univariate[dfType$type == "all","backtransform"] <- FALSE
                }
            }
        }
        out$glht <- stats::setNames(vector(mode = "list", length = n.grid), rownames(grid))
    }

    ## ** Wald tests
    for(iG in 1:NROW(grid)){ ## iG <- 1
        iType <- grid[iG,"type"]
        iTerm <- grid[iG,"term"]
        
        iContrast <- contrast[[iType]][[iTerm]]
        iN.hypo <- NROW(iContrast)
        iNull <- null[[iType]][[iTerm]]
        
        ## *** Multivariate Wald test
        if(multivariate){
            
            iSimplify <- simplifyContrast(iContrast, iNull) ## remove extra lines
            out$multivariate[iG,"df.num"] <- iSimplify$dim
            
            iC.vcov.C_M1 <- try(solve(iSimplify$C %*% vcov.param %*% t(iSimplify$C)), silent = TRUE)
                
            if(inherits(iC.vcov.C_M1,"try-error")){
                
                out$multivariate[iG,"message"] <- "Could not invert the covariance matrix for the proposed contrast."
                
            }else{
                
                out$multivariate[iG,"statistic"] <- as.double(t(iSimplify$C %*% param - iSimplify$rhs) %*% iC.vcov.C_M1 %*% (iSimplify$C %*% param - iSimplify$rhs))/iSimplify$dim 

                ## degree of freedom
                if(df>0){

                    iSVD <- eigen(iC.vcov.C_M1)
                    iSVD.D <- diag(iSVD$values, nrow = iSimplify$dim, ncol = iSimplify$dim)
                    iSVD.P <- iSVD$vectors
                    iSVD.contrast <- sqrt(iSVD.D) %*% t(iSVD.P) %*% iSimplify$C
                    colnames(iSVD.contrast) <- colnames(iSimplify$C)

                    iNu_m <- dfSigma(contrast = iSVD.contrast,
                                     vcov = vcov.param,
                                     dVcov = dVcov.param,
                                     keep.param = colnames(iSimplify$C))
                
                    iEQ <- sum(iNu_m/(iNu_m - 2))
                    out$multivariate[iG,"df.denom"] <- 2 * iEQ/(iEQ - iSimplify$dim)
                    
                }

            }
        }

        ## *** Univariate Wald test
        if(univariate){

            iG.univariate <- which(out$univariate$name == rownames(out$multivariate)[iG])

            if(df>0){
                out$univariate[iG.univariate,"df"] <- as.double(.dfX(X.beta = iContrast, vcov.param = vcov.param, dVcov.param = dVcov.param))
            }else{
                out$univariate[iG.univariate,"df"] <- rep(Inf, length(iG.univariate))
            }
            out$univariate[iG.univariate,"estimate"] <- as.double(iContrast %*% param)
            out$univariate[iG.univariate,"se"] <- sqrt(diag(iContrast %*% vcov.param %*% t(iContrast)))
            out$univariate[iG.univariate,"null"] <- iNull
    
            out$glht[[iG]] <- multcomp::glht(object, linfct = iContrast, rhs = iNull, df = ceiling(stats::median(out$univariate[iG.univariate,"df"])),
                                             coef. = function(iX){coef.lmm(iX, effects = effects)},
                                             vcov. = function(iX){vcov.lmm(iX, robust = robust, effects = effects)})
            out$glht[[iG]]$model <- NULL
        }

    }

    if(multivariate){
        out$multivariate$p.value <- 1 - stats::pf(out$multivariate$statistic, df1 = out$multivariate$df.num, df2 = out$multivariate$df.denom)
    }
    if(univariate){
        out$univariate$statistic <- (out$univariate$estimate-out$univariate$null)/out$univariate$se
    }
browser()    
    ## ** save some of the objects for possible use of rbind.Wald_lmm
    if(!simplify){
        globalC <- do.call(rbind, lapply(unlist(out$glht, recursive = FALSE), "[[", "linfct"))
        
        if(attr(robust,"call")){
            out$vcov <- globalC %*% vcov.param %*% t(globalC)
            robust2 <- robust
        }else{
            ## default: non-robust vcov and robust for iid
            if(robust){
                out$vcov <- globalC %*% vcov(object, df = FALSE, effects = "all", robust = FALSE,
                                             transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE) %*% t(globalC)
            }else{
                out$vcov <- globalC %*% vcov.param %*% t(globalC)
            }
            robust2 <- TRUE
        }
        out$iid <- iid(object, effects = effects, robust = robust2, 
                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE) %*% t(globalC)
    }


    ## ** export
    class(out) <- append("Wald_lmm",class(out))
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

    ## ** check nesting
    testEqual <- .checkNesting(objectH0, objectH1)
    rhs <- attr(testEqual,"rhs")

    ## ** objective function
    if(objectH0$args$method.fit!=objectH1$args$method.fit){
        stop("The two models should use the same type of objective function for the likelihood ratio test to be valid. \n")
    }
     if(objectH1$args$method.fit=="REML" && (testEqual["mean"]==FALSE)){
        objectH0$call$method.fit <- "ML"
        objectH1$call$method.fit <- "ML"
        if(testEqual["var"] && testEqual["cor"]){
            message("Cannot use a likelihood ratio test to compare mean parameters when the objective function is REML. \n",
                    "Will re-estimate the model via ML and re-run the likelihood ratio test. \n")
        }else{
            message("Cannot use a likelihood ratio test to compare mean parameters when the objective function is REML. \n",
                    "Will re-estimate the model via ML and re-run the likelihood ratio test. \n",
                    "This will affect the estimation of the variance and correlation parameters. \n")
        }
        out <- .anova_LRT(eval(objectH0$call),eval(objectH1$call))
        attr(out,"type") <- type
        return(out)
    }

    ## ** LRT
    name.paramH0 <- names(coef(objectH0, effects = "all"))   
    name.paramH1 <- names(coef(objectH1, effects = "all"))
    n.paramTest <- length(name.paramH1)-length(name.paramH0)

    if(is.null(rhs)){
        null <- ""
    }else{
        null <- paste(paste0(names(rhs),"==",rhs), collapse = "\n                   ")
    }
    out <- data.frame(null = null,
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
    class(out) <- append("LRT_lmm",class(out))
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


##----------------------------------------------------------------------
### anova.R ends here
