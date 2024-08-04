### anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:38) 
## Version: 
## Last-Updated: Aug  4 2024 (16:40) 
##           By: Brice Ozenne
##     Update #: 2024
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
##' Can also be \code{2} compute the degrees of freedom w.r.t. robust standard errors instead of w.r.t. model-based standard errors.
##' @param df [logical] Should degrees of freedom be estimated using a Satterthwaite approximation?
##' If yes F-distribution (multivariate) and Student's t-distribution (univariate) are used.
##' Other chi-squared distribution and normal distribution are used.
##' @param univariate [logical] Should an estimate, standard error, confidence interval, and p-value be output for each hypothesis?
##' @param multivariate [logical] Should all hypotheses be simultaneously tested using a multivariate Wald test?
##' @param transform.sigma,transform.k,transform.rho are passed to the \code{vcov} method. See details section in \code{\link{coef.lmm}}.
##' @param simplify [logical] should only the estimates, variance-covariance with its gradient, and degrees of freedom relative to the parameters involved in the Wald test be stored (TRUE)
##' or relative to all model parameters (0.5) along with their iid decomposition (FALSE).
##' 
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A data.frame (LRT) or a list of containing the following elements (Wald):\itemize{
##' \item \code{multivariate}: data.frame containing the multivariate Wald test.
##' The column \code{df.num} refers to the degrees of freedom for the numerator (i.e. number of hypotheses)
##' wherease the column \code{df.denum} refers to the degrees of freedom for the denominator (i.e. Satterthwaite approximation).
##' \item \code{univariate}: data.frame containing each univariate Wald test.
##' \item \code{glht}: used internally to call functions from the multcomp package.
##' \item \code{object}: list containing key information about the linear mixed model.
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
anova.lmm <- function(object, effects = NULL, rhs = NULL, robust = NULL, df = NULL,
                      univariate = TRUE, multivariate = TRUE, 
                      transform.sigma = NULL, transform.k = NULL, transform.rho = NULL,
                      simplify = NULL, ...){

    call <- match.call()
    options <- LMMstar.options()

    ## ** normalized user input    
    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    if(!is.matrix(effects) && !is.null(rhs)){
        message("Argument \'rhs\' is ignored unless argument \"effect\" is a contrast matrix. \n")
    }

    ## ** run test
    if(inherits(effects,"lmm")){
        
        ## *** Likelihood Ratio Test (LRT)
        out <- .anova_LRT(object1 = object, object2 = effects)
        
    }else{
        
        ## *** Wald test
        if(!univariate & !multivariate){
            stop("Argument \'univariate\' and \'multivariate\' should not be both FALSE. \n")
        }

        ## initialize effects
        if(is.null(effects)){
            effects <- options$effects
        }

        ## initialize robust
        if(is.null(robust)){
            robust <- FALSE
        }
        
        ## initialize df
        if(is.null(df)){
            df <- !is.null(object$df)
        }

        ## initialize simplify
        if(is.null(simplify)){
            if((df == FALSE) || (is.character(effects) && all(effects %in% c("all","mean","fixed","variance","correlation")))){
                simplify <- TRUE
            }else{
                simplify <- 0.5
            }
        }else{
            if(!is.numeric(simplify) && !is.logical(simplify)){
                stop("Argument \'simplify\' must be numeric or logical. \n")
            }
            if(length(simplify)!=1){
                stop("Argument \'simplify\' must have length 1. \n")
            }
            if(simplify %in% c(0,0.5,1) == FALSE){
                stop("Argument \'simplify\' must be TRUE/1, 0.5, or FALSE/0. \n")
            }
        }
        
        ## initialize tranformation
        init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                                x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho,
                                simplify = FALSE)
        transform.sigma <- init$transform.sigma
        transform.k <- init$transform.k
        transform.rho <- init$transform.rho
        
        ## generate contrast matrix
        ls.e2c <- effects2contrast(object, effects = effects, rhs = rhs,
                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)

        ## extract model coefficient and uncertainty
        object.param <- model.tables(object, effects = "param")
        type.param <- stats::setNames(object.param$type,object.param$name)
        if(transform.k %in% c("sd","logsd","var","logvar")){
            type.param[type.param=="sigma"] <- "k"
        }
        param <- coef(object, effects = "all",
                      transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = TRUE)
        param.notrans <- coef(object, effects = "all",
                              transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = FALSE)
        vcov.param <- vcov(object, df = df, effects = list("all",c("all","gradient"))[[df+1]], robust = robust,
                           transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
        dVcov.param <- attr(vcov.param,"gradient")
        attr(vcov.param,"gradient") <- NULL
        attr(vcov.param,"df") <- NULL
        if(simplify == FALSE){
            iid.param <- iid(object, effects = "all", robust = robust, 
                             transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
        }else{
            iid.param <- NULL
        }

        ## call Wald test
        out <- .anova_Wald(param = param, param.notrans = param.notrans, vcov.param = vcov.param, dVcov.param = dVcov.param, iid.param = iid.param, type.param = type.param,
                           contrast = ls.e2c$contrast, null = ls.e2c$null, robust = robust, df = df,
                           multivariate = multivariate, univariate = univariate, 
                           transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, backtransform = ls.e2c$backtransform,
                           simplify = simplify)

        ## add extra information to object that may be using by rbind
        out$object <- list(outcome = object$outcome$var,
                           method.fit = object$args$method.fit,
                           type.information = object$args$type.information,
                           cluster.var = object$cluster$var,
                           cluster = object$cluster$levels)

    }

    ## ** export
    attr(out,"call") <- call
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
        objectH0 <- eval(objectH0$call)
        objectH1 <- eval(objectH1$call)
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

## * .anova_Wald
.anova_Wald <- function(param, param.notrans, vcov.param, iid.param, dVcov.param, type.param,
                        contrast, null, robust, df,
                        univariate, multivariate, 
                        transform.sigma, transform.k, transform.rho, backtransform,
                        simplify){

    ## ** prepare
    ## *** robust not equal to 1: use current vcov matrix (when =1: vcov.param robust whereas attr(.,"vcov") model based)
    if(df && is.null(attr(dVcov.param,"vcov"))){
        df_vcov.param <- vcov.param
    }else{
        df_vcov.param <- attr(dVcov.param,"vcov")
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
    ## grid combining type (mean,variance,covariance,user) and terms (e.g. time, group, time:group, age)
    ls.grid <- lapply(names(contrast), function(iName){ ## iName <- names(contrast)[1]
        iDf <- data.frame(type.original = iName, type = NA, term = names(contrast[[iName]]), n.test = sapply(contrast[[iName]],NROW))
        return(iDf)
    })
    grid <- do.call(rbind,ls.grid)
    if(NROW(grid)>1){
        rownames(grid) <- paste(grid$type.original, grid$term, sep = "_")
    }
    n.grid <- NROW(grid)

    if(all(grid$type.original=="user")){
        ## single contrast matrix
        M.type <- do.call(rbind,apply(contrast[[1]][[1]], MARGIN = 1, FUN = function(iRow){ ## iRow <- contrast[[1]][[1]][1,]
            table(factor(type.param[names(which(iRow!=0))], levels = c("mu","sigma","k","rho")))
        }, simplify = FALSE))
        table.type <- cbind(as.data.frame(M.type), overall = apply(M.type, MARGIN = 1, function(iRow){ifelse(sum(iRow!=0)==1,names(which(iRow!=0)),"all")}))
        if(length(unique(table.type$overall))==1){
            grid$type <- table.type$overall[1]
        }else{
            grid$type <- "all"
        }
        names(contrast) <- grid$type
        names(null) <- grid$type
    }else{ 
        grid$type <- grid$type.original
    }

    ## *** update null according to transformation
    if(all(grid$type.original=="user")){

        message.error <- c("Consider not using any transformation (setting \'transform.sigma\', \'transform.k\', and \'transform.rho\' to \"none\"), \n",
                           "or only considering the transformed scale (specify \'transform.sigma\', \'transform.k\', \'transform.rho\' and explicit the contrast in term of the transformed parameters). \n")

        ## when back-transformed is FALSE, the user has specified the transformation and the transformed parameters so the null should also be on the right scale
        if(backtransform && transform2.sigma=="log" && any(table.type$overall=="sigma")){

            nt.sigma <- .null_transform(contrast = contrast[[1]]$user[table.type$overall=="sigma",,drop=FALSE],
                                        null = null[[1]]$user[table.type$overall=="sigma"],
                                        transformation = "log",
                                        n.param = table.type[table.type$overall=="sigma","sigma"],
                                        message.error = message.error, arg.rhs = "rhs")

            contrast[[1]]$user[table.type$overall=="sigma",] <- nt.sigma$contrast
            null[[1]]$user[table.type$overall=="sigma"] <- nt.sigma$null

        }
        if(backtransform && transform2.k=="log" && any(table.type$overall=="k")){

            nt.k <- .null_transform(contrast = contrast[[1]]$user[table.type$overall=="k",,drop=FALSE],
                                    null = null[[1]]$user[table.type$overall=="k"],
                                    transformation = "log",
                                    n.param = table.type[table.type$overall=="k","k"],
                                    message.error = message.error, arg.rhs = "rhs")

            contrast[[1]]$user[table.type$overall=="k",] <- nt.k$contrast
            null[[1]]$user[table.type$overall=="k"] <- nt.k$null

        }
        if(backtransform && transform2.rho=="atanh" && any(table.type$overall=="atanh")){
    
            nt.rho <- .null_transform(contrast = contrast[[1]]$user[table.type$overall=="rho",,drop=FALSE],
                                      null = null[[1]]$user[table.type$overall=="rho"],
                                      transformation = "atanh",
                                      n.param = table.type[table.type$overall=="rho","rho"],
                                      message.error = message.error, arg.rhs = "rhs")

            contrast[[1]]$user[table.type$overall=="rho",] <- nt.rho$contrast
            null[[1]]$user[table.type$overall=="rho"] <- nt.rho$null

        }

        if(any(table.type$overall=="all")){

            type.all <- names(which(colSums(table.type[table.type$overall=="all",c("mu","sigma","k","rho"),drop=FALSE])>0))
            transform2.all <- unique(c(mu = "none", sigma = transform.sigma, k = transform.k, rho = transform.rho)[type.all])
            if(backtransform && length(transform2.all)>1){
                stop("Cannot move from the original scale to the transformed scale when testing linear hypothesis involving parameters with different transformations. \n",
                     message.error)
            }else if(backtransform && transform2.all != "none"){
                nt.all <- .null_transform(contrast = contrast[[1]]$user[table.type$overall=="all",,drop=FALSE],
                                          null = null[[1]]$user[table.type$overall=="all"],
                                          transformation = transform2.all,
                                          n.param = rowSums(table.type[table.type$overall=="all",c("mu","sigma","k","rho")]),
                                          message.error = message.error, arg.rhs = "rhs")
                contrast[[1]]$user[table.type$overall=="all",] <- nt.all$contrast
                null[[1]]$user[table.type$overall=="all"] <- nt.all$null
            }

        }else{

            transform2.all <- "none"

        }

    }else{
        
        if(transform2.k=="log"){
            null$k <- lapply(null$k, log)
        }
        if(transform2.rho=="atanh"){ ## should have no effect since atanh(0)=0
            null$rho <- lapply(null$rho, atanh)
        }
        transform2.all <- "none"

    }

    ## *** output
    out <- list(args = data.frame(type = ifelse(all(grid$type.original=="user"),"user","auto"), robust = robust, df = df,
                                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.all = transform2.all,
                                  univariate = univariate, multivariate = multivariate, simplify = as.numeric(simplify))
                )

    if(multivariate){
        out$multivariate <- data.frame(matrix(NA, nrow = n.grid, ncol = 7,
                                              dimnames = list(NULL,c("type","term","null","statistic","df.num","df.denom","p.value"))))
        out$multivariate$type <- grid$type
        out$multivariate$term <- grid$term
        out$multivariate$null <- unname(sapply(unlist(contrast, recursive = FALSE), FUN = function(iRow){paste(rownames(iRow), collapse = ", ")}))
        rownames(out$multivariate) <- rownames(grid)
    }
    if(univariate){

        out$glht <- stats::setNames(vector(mode = "list", length = n.grid), rownames(grid))
        ## name (e.g. mean_visit, var_visit) may differ from term (visit, visit) when testing mean, variance, and correlation parameters which may all refer to the same variable, say time
        out$univariate <- data.frame(matrix(NA, nrow = sum(grid$n.test), ncol = 15,
                                            dimnames = list(NULL,c("type","term","name","n.param","estimate","se","df","statistic","null","p.value","quantile","lower","upper",
                                                                   "transformed","tobacktransform"))))

        out$univariate$term <- unname(unlist(mapply(iTerm = grid$term, iN = grid$n.test, FUN = function(iTerm,iN){rep(iTerm,iN)}, SIMPLIFY = FALSE)))

        if(all(grid$type.original=="user")){
            out$univariate$type <- table.type$overall
            out$univariate$n.param <- rowSums(table.type[c("mu","sigma","k","rho")])
            out$univariate$name <- rownames(grid)
        }else{
            out$univariate$type <- unname(unlist(mapply(iType = grid$type, iN = grid$n.test, FUN = function(iType,iN){rep(iType,iN)}, SIMPLIFY = FALSE)))
            out$univariate$n.param <- 1
            out$univariate$name <- rownames(grid)[merge(x = out$univariate[c("type","term")],
                                                        y = cbind(index = 1:n.grid, grid[c("type","term")]),
                                                        by = c("type","term"), sort = FALSE)$index]
        }
        rownames(out$univariate) <- unname(unlist(sapply(unlist(contrast, recursive = FALSE), rownames)))

        ## back-transformation
        ## no need when mean parameters
        out$univariate[out$univariate$type == "mu",c("transformed","tobacktransform")] <- FALSE
        ## potential log-transformation/exp-backtransform for linear combinations including only sigma parameters
        out$univariate[out$univariate$type == "sigma", "transformed"] <- (transform.sigma!="none")
        out$univariate[out$univariate$type == "sigma", "tobacktransform"] <- backtransform & (transform2.sigma=="log")
        ## potential log-transformation/exp-backtransform for linear combinations including only k parameters
        out$univariate[out$univariate$type == "k", "transformed"] <- (transform.k!="none")
        out$univariate[out$univariate$type == "k", "tobacktransform"] <- backtransform & (transform2.k=="log")
        ## potential atanh-transformation/tanh-backtransform for linear combinations including only rho parameters
        out$univariate[out$univariate$type == "rho","transformed"] <- ((out$univariate[out$univariate$type == "rho","n.param"] == 1) | !backtransform) & (transform.rho!="none")
        out$univariate[out$univariate$type == "rho","tobacktransform"] <- (out$univariate[out$univariate$type == "rho","n.param"] == 1) & backtransform & (transform2.rho=="atanh")
        ## it is checked previously that all transformations are the same (i.e. k sigma parameters)
        out$univariate[out$univariate$type == "all","transformed"] <- (transform2.all!="none")
        out$univariate[out$univariate$type == "all","tobacktransform"] <- backtransform & (transform2.all!="none")
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
                attr(out$multivariate,"statistic") <- "Could not invert the covariance matrix for the proposed contrast."
            }else{
                out$multivariate[iG,"statistic"] <- as.double(t(iSimplify$C %*% param - iSimplify$rhs) %*% iC.vcov.C_M1 %*% (iSimplify$C %*% param - iSimplify$rhs))/iSimplify$dim
            }

            ## degree of freedom
            if(df>0 && is.null(attr(dVcov.param,"vcov"))){
                df_iC.vcov.C_M1 <- iC.vcov.C_M1
            }else if(df>0 && !is.null(attr(dVcov.param,"vcov"))){
                ## case of robust s.e. but df evaluated on model-based s.e.
                df_iC.vcov.C_M1 <- try(solve(iSimplify$C %*% df_vcov.param %*% t(iSimplify$C)), silent = TRUE)
            }

            if(df==0 || inherits(df_iC.vcov.C_M1,"try-error")){
                out$multivariate[iG,"df.denom"] <- Inf
            }else{
                iSVD <- eigen(df_iC.vcov.C_M1, symmetric = TRUE)
                iSVD.D <- diag(iSVD$values, nrow = iSimplify$dim, ncol = iSimplify$dim)
                iSVD.P <- iSVD$vectors
                iSVD.contrast <- sqrt(iSVD.D) %*% t(iSVD.P) %*% iSimplify$C
                colnames(iSVD.contrast) <- colnames(iSimplify$C)

                iNu_m <- .df_contrast(contrast = iSVD.contrast[,dimnames(dVcov.param)[[1]],drop=FALSE], ## subset for rbind when dVcov.param is not full in the first 2 dimensions
                                      vcov.param = df_vcov.param,
                                      dVcov.param = dVcov.param)
                
                    iEQ <- sum(iNu_m/(iNu_m - 2))
                    out$multivariate[iG,"df.denom"] <- 2 * iEQ/(iEQ - iSimplify$dim)
            }
        }

        ## *** Univariate Wald test
        if(univariate){

            iG.univariate <- which(out$univariate$name == rownames(grid)[iG])

            if(df>0){
                
                out$univariate[iG.univariate,"df"] <- as.double(.df_contrast(contrast = iContrast[,dimnames(dVcov.param)[[1]],drop=FALSE], ## subset for rbind when dVcov.param is not full
                                                                             vcov.param = df_vcov.param,
                                                                             dVcov.param = dVcov.param))
                if(multivariate){
                    if(out$multivariate[iG,"df.denom"]<1){
                        ## warning("Suspicious degree of freedom was found for the F-statistic (",out$multivariate[iG,"df.denom"],"). \n",
                        ##         "It has been set to the average degree of freedom of the univariate Wald tests (",mean(out$univariate[iG.univariate,"df"]),") instead. \n")
                        attr(out$multivariate,"df") <- "average degree of freedom of the univariate Wald tests"
                        out$multivariate[iG,"df.denom"] <- mean(out$univariate[iG.univariate,"df"])
                    }
                }
            }else{
                out$univariate[iG.univariate,"df"] <- rep(Inf, length(iG.univariate))
            }
            ## handle the case of linear combinations of correlation where no transformation should be made for the estimate
            iTrans.univariate <- intersect(which(out$univariate$transformed),iG.univariate)
            names(iTrans.univariate) <- rownames(out$univariate)[iTrans.univariate]
            iNtrans.univariate <- intersect(which(!out$univariate$transformed),iG.univariate)
            names(iNtrans.univariate) <- rownames(out$univariate)[iNtrans.univariate]

            out$univariate[iTrans.univariate,"estimate"] <- as.double(iContrast[names(iTrans.univariate),,drop=FALSE] %*% param)
            out$univariate[iNtrans.univariate,"estimate"] <- as.double(iContrast[names(iNtrans.univariate),,drop=FALSE] %*% param.notrans)
            out$univariate[iG.univariate,"se"] <- sqrt(diag(iContrast %*% vcov.param %*% t(iContrast)))
            out$univariate[iG.univariate,"null"] <- iNull

            ## create glht object
            out$glht[[iG]] <- list(model = NULL,
                                   linfct = iContrast,
                                   rhs = iNull,
                                   coef = stats::setNames(param, names(param.notrans)),
                                   vcov = vcov.param,
                                   df = round(stats::median(out$univariate[iG.univariate,"df"])),
                                   alternative = "two.sided",
                                   coef.name = names(param), ## possibly transformed named (coef will always be named without transformed names)
                                   coef.notrans = param.notrans, ## values without transformation
                                   coef.type = type.param,
                                   dVcov = dVcov.param,
                                   iid = iid.param)

            if(simplify == TRUE){  
                
                index.simplify <- which(colSums(iContrast!=0)>0)
                out$glht[[iG]]$linfct <- out$glht[[iG]]$linfct[,index.simplify,drop=FALSE]
                out$glht[[iG]]$coef <- out$glht[[iG]]$coef[index.simplify] ## transformation.name = FALSE (use transformed values but not transformed name)
                out$glht[[iG]]$vcov <- out$glht[[iG]]$vcov[index.simplify,index.simplify,drop=FALSE]
                out$glht[[iG]]$coef.name <- out$glht[[iG]]$coef.name[index.simplify] ## transformation.name = TRUE (use transformed values and transformed name)
                out$glht[[iG]]$coef.notrans <- out$glht[[iG]]$coef.notrans[index.simplify] ## transformation.sigma = "none", ... (use original value and original name)
                out$glht[[iG]]$coef.type <- out$glht[[iG]]$coef.type[index.simplify]
                out$glht[[iG]]$dVcov <- NULL ## useless to have dVcov if we do not have the full vcov for df computation
                out$glht[[iG]]$iid <- out$glht[[iG]]$iid[,index.simplify,drop=FALSE]

            }else if(simplify == 0.5){

                name.simplify <- names(which(colSums(iContrast!=0)>0))
                out$glht[[iG]]$dVcov <- out$glht[[iG]]$dVcov[name.simplify,name.simplify,,drop=FALSE]

            }

            class(out$glht[[iG]]) <- "glht"
        }
    }

    
    if(multivariate){
        out$multivariate$p.value <- 1 - stats::pf(out$multivariate$statistic, df1 = out$multivariate$df.num, df2 = out$multivariate$df.denom)
    }
    if(univariate){
        out$univariate$statistic <- (out$univariate$estimate-out$univariate$null)/out$univariate$se
    }

    ## ** export
    class(out) <- append("Wald_lmm",class(out))
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

##----------------------------------------------------------------------
### anova.R ends here
