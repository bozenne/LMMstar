### anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:38) 
## Version: 
## Last-Updated: okt 24 2025 (17:33) 
##           By: Brice Ozenne
##     Update #: 2234
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * anova.lmm (documentation)
##' @title Multivariate Tests For a Linear Mixed Model
##' @description Test linear combinations of parameters from a linear mixed model
##' using Wald test or Likelihood Ratio Test (LRT). 
##' 
##' @param object a \code{lmm} object. Only relevant for the anova function.
##' @param effects [character or numeric matrix] Should the Wald test be computed for all variables (\code{"all"}),
##' or only variables relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only variables relative to the variance structure (\code{"variance"}),
##' or only variables relative to the correlation structure (\code{"correlation"}).
##' Can also be use to specify linear combinations of coefficients or a contrast matrix, similarly to the \code{linfct} argument of the \code{multcomp::glht} function.
##' @param rhs [numeric vector] the right hand side of the hypothesis. Only used when the argument \code{effects} is a matrix.
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors. 
##' Can also be \code{2} compute the degrees-of-freedom w.r.t. robust standard errors instead of w.r.t. model-based standard errors.
##' @param type.information [character] Should the expected information be used  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param df [logical] Should degrees-of-freedom be estimated using a Satterthwaite approximation?
##' If yes F-distribution (multivariate) and Student's t-distribution (univariate) are used.
##' Other chi-squared distribution and normal distribution are used.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the variance-covariance matrix. Only relevant if differs from the fitted values.
##' @param univariate [logical] Should an estimate, standard error, confidence interval, and p-value be output for each hypothesis?
##' @param multivariate [logical] Should all hypotheses be simultaneously tested using a multivariate Wald test?
##' @param transform.sigma,transform.k,transform.rho are passed to the \code{vcov} method. See details section in \code{\link{coef.lmm}}.
##' @param simplify [logical] when \code{FALSE} the lmm object is stored in the output.
##' 
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A data.frame (LRT) or a list of containing the following elements (Wald):\itemize{
##' \item \code{args}: list containing argument values from the function call.
##' \item \code{multivariate}: data.frame containing the multivariate Wald test.
##' The column \code{df.num} refers to the degrees-of-freedom for the numerator (i.e. number of hypotheses)
##' wherease the column \code{df.denum} refers to the degrees-of-freedom for the denominator (i.e. Satterthwaite approximation).
##' \item \code{univariate}: data.frame containing each univariate Wald test.
##' \item \code{glht}: used internally to call functions from the multcomp package.
##' \item \code{model} (optional): linear mixed model used to generate the Wald tests.
##' }
##' 
##' @details By default adjustment of confidence intervals and p-values for multiple comparisons is based on the distribution of the maximum-statistic.
##' This is refered to as a single-step Dunnett multiple testing procedures in table II of Dmitrienko et al. (2013).
##' It is performed using the multcomp package with the option \code{test = adjusted("single-step")} with equal degrees-of-freedom
##' or by simulation using a Student's t copula with unequal degrees-of-freedom (more in the note of the details section of \code{\link{confint.Wald_lmm}}).
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
anova.lmm <- function(object, effects = NULL, rhs = NULL, type.information = NULL, robust = NULL, df = NULL, 
                      univariate = TRUE, multivariate = TRUE, p = NULL,
                      transform.sigma = NULL, transform.k = NULL, transform.rho = NULL,
                      simplify = NULL, ...){

    mycall <- match.call()

    ## ** special case (Likelihood Ratio Test - LRT)
    if(inherits(effects,"lmm")){
        ## hidden argument force
        out <- .anova_LRT(object1 = object, object2 = effects, ...)
        attr(out,"call") <- mycall
        return(out)
    }        

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** effects
    if(is.null(effects)){
        effects <- options$effects
    }else if(!is.matrix(effects) && !is.null(rhs)){
        message("Argument \'rhs\' is ignored unless argument \"effect\" is a contrast matrix. \n")
    }

    ## *** type.information
    if(is.null(type.information)){
        type.information <- object$args$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    ## ***  robust
    if(is.null(robust)){
        robust <- FALSE
    }
    
    ## ***  df
    if(is.null(df)){
        df <- !is.null(object$df)
    }

    ## *** univariate and multivariate
    if(!univariate & !multivariate){
        stop("Argument \'univariate\' and \'multivariate\' should not be both FALSE. \n")
    }

    ## *** simplify
    if(is.null(simplify)){
        if((is.character(effects) && all(effects %in% c("all","mean","fixed","variance","correlation")))){ ## automatic contrast
            simplify <- TRUE
        }else{ ## default may change in the future?
            simplify <- TRUE
        }
    }else{
        if(!is.numeric(simplify) && !is.logical(simplify)){
            stop("Argument \'simplify\' must be numeric or logical. \n")
        }
        if(length(simplify)!=1){
            stop("Argument \'simplify\' must have length 1. \n")
        }
        if(!is.logical(simplify) && simplify %in% c(0,1) == FALSE){
            stop("Argument \'simplify\' must be TRUE or FALSE. \n")
        }
    }
    
    ## *** transformation & p
    init <- .init_transform(p = p, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho,
                            table.param = object$design$param)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform
    theta <- init$p

    ## ** generate contrast matrix
    ls.e2c <- effects2contrast(object, effects = effects, rhs = rhs, options = options)

    ## ** extract model coefficient and uncertainty
    table.param <- stats::model.tables(object, effects = "param", options = options)
    type.param <- stats::setNames(table.param$type,table.param$name)
    if(transform.k %in% c("sd","logsd","var","logvar")){
        type.param[type.param=="sigma"] <- "k"
    }
    param <- stats::coef(object, p = theta, effects = "all",
                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = TRUE, options = options)
    param.notrans <- stats::coef(object, p = theta, effects = "all",
                                 transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = FALSE, options = options)
    vcov.param <- stats::vcov(object, p = theta, df = df, effects = list("all",c("all","gradient"))[[df+1]], type.information = type.information, robust = robust, 
                              transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE, options = options)
    dVcov.param <- attr(vcov.param,"gradient")
    attr(vcov.param,"gradient") <- NULL

    ## ** run Wald test
    out <- .anova_Wald(param = param, param.notrans = param.notrans, vcov.param = vcov.param, dVcov.param = dVcov.param, type.param = type.param,
                       contrast = ls.e2c$contrast, null = ls.e2c$null, df = df,
                       multivariate = multivariate, univariate = univariate, 
                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, backtransform = ls.e2c$backtransform,
                       simplify = simplify)
    
    ## ** add extra information to object, e.g. to retrieve the original contrast matrix
    out$args$p.null <- is.null(p) ## original lmm estimates
    out$args$method.fit <- object$args$method.fit
    out$args$robust <- robust
    out$args$type.information <- type.information
    if(simplify == FALSE){
        out$model <- object
    }else{
        ## out$model <- parent.env(environment()) ## save pointer to the environment
    }

    out$param = cbind(trans.value = param, value = param.notrans, table.param[c("trans.name","name","type")])
    rownames(out$param) <- NULL

    ## ** export
    out$call <- mycall
    return(out)
}

## * .anova_LRT
.anova_LRT <- function(object1, object2, force = FALSE){
    tol <- 1e-10
    
    ## ** normalize user input
    ## *** re-order models
    logLik1 <- stats::logLik(object1)
    logLik2 <- stats::logLik(object2)
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
    if(force){
        testEqual <- NULL
        rhs <- NULL
    }else{
        testEqual <- .checkNesting(objectH0, objectH1)
        rhs <- attr(testEqual,"rhs")
    }

    ## ** objective function
    if(objectH0$args$method.fit!=objectH1$args$method.fit){
        stop("The two models should use the same type of objective function for the likelihood ratio test to be valid. \n")
    }
    if(!force && objectH1$args$method.fit=="REML" && (testEqual["mean"]==FALSE)){
        if(testEqual["var"] && testEqual["cor"]){
            message("Cannot use a likelihood ratio test to compare mean parameters when the objective function is REML. \n",
                    "Will re-estimate the model via ML and re-run the likelihood ratio test. \n")
        }else{
            message("Cannot use a likelihood ratio test to compare mean parameters when the objective function is REML. \n",
                    "Will re-estimate the model via ML and re-run the likelihood ratio test. \n",
                    "This will affect the estimation of the variance and correlation parameters. \n")
        }
        objectH0 <- stats::update(objectH0, method.fit = "ML")
        objectH1 <- stats::update(objectH1, method.fit = "ML")
    }

    ## ** LRT
    name.paramH0 <- stats::model.tables(objectH0, effects = "param")$name   
    name.paramH1 <- stats::model.tables(objectH1, effects = "param")$name
    n.paramTest <- length(name.paramH1)-length(name.paramH0)

    if(is.null(rhs) || force){
        null <- ""
    }else{
        null <- paste(paste0(names(rhs),"==",rhs), collapse = paste("\n",  strrep(" ",nchar("  Null hypothesis: "))))
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
.anova_Wald <- function(param, param.notrans, vcov.param, dVcov.param, type.param,
                        contrast, null, df,
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
        if(backtransform && transform2.rho=="atanh" && any(table.type$overall=="rho")){
            
            nt.rho <- .null_transform(contrast = contrast[[1]]$user[table.type$overall=="rho",,drop=FALSE],
                                      null = null[[1]]$user[table.type$overall=="rho"],
                                      transformation = "atanh",
                                      n.param = table.type[table.type$overall=="rho","rho"],
                                      message.error = message.error, arg.rhs = "rhs")

            contrast[[1]]$user[table.type$overall=="rho",] <- nt.rho$contrast
            null[[1]]$user[table.type$overall=="rho"] <- nt.rho$null

        }

        if(any(table.type$overall=="all")){
            type.sum <- colSums(table.type[table.type$overall=="all",c("mu","sigma","k","rho"),drop=FALSE])
            if(all(type.sum==0)){
                type.all <- "mu"
            }else{
                type.all <- names(which(type.sum>0))
            }
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
    out <- list(args = data.frame(type = ifelse(all(grid$type.original=="user"),"user","auto"), method.fit = NA, type.information = NA, robust = NA, df = df,
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
                                            dimnames = list(NULL,c("type","term","name","n.param","estimate","se","df","quantile","lower","upper","statistic","null","p.value",
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
        rownames(out$univariate) <- unname(unlist(lapply(unlist(contrast, recursive = FALSE), rownames)))
        
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
        out$univariate[out$univariate$type == "rho","transformed"] <- (transform.rho!="none")
        out$univariate[out$univariate$type == "rho","tobacktransform"] <- backtransform & (transform2.rho=="atanh") & (out$univariate[out$univariate$type == "rho","n.param"] == 1)
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

            ## degree-of-freedom
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

                iNu_m <- .df_contrast(contrast = iSVD.contrast,
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

                out$univariate[iG.univariate,"df"] <- as.double(.df_contrast(contrast = iContrast,
                                                                             vcov.param = df_vcov.param,
                                                                             dVcov.param = dVcov.param))
                if(multivariate){
                    if(is.na(out$multivariate[iG,"df.denom"])){
                        if(any(iSVD$value <= 0)){
                            warning("Variance-covariance matrix relative to the linear contrasts is not positive definite. \n")
                        }
                    }else if(out$multivariate[iG,"df.denom"]<1){
                        ## warning("Suspicious degrees-of-freedom was found for the F-statistic (",out$multivariate[iG,"df.denom"],"). \n",
                        ##         "It has been set to the average degree-of-freedom of the univariate Wald tests (",mean(out$univariate[iG.univariate,"df"]),") instead. \n")
                        attr(out$multivariate,"df") <- "average degree-of-freedom of the univariate Wald tests"
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

            if(length(iTrans.univariate)>0){
                out$univariate[iTrans.univariate,"estimate"] <- as.double(iContrast[names(iTrans.univariate),,drop=FALSE] %*% param)
            }
            if(length(iNtrans.univariate)>0){
                out$univariate[iNtrans.univariate,"estimate"] <- as.double(iContrast[names(iNtrans.univariate),,drop=FALSE] %*% param.notrans)
            }
            out$univariate[iG.univariate,"se"] <- sqrt(diag(iContrast %*% vcov.param %*% t(iContrast)))
            out$univariate[iG.univariate,"null"] <- iNull

            ## create glht object
            ## dVcov is not stored as it is only useful to retrieve df when having the full vcov (not only the contrast vcov)            
            index.simplify <- which(colSums(iContrast!=0)>0)
            out$glht[[iG]] <- list(model = NULL,
                                   linfct = iContrast[,index.simplify,drop=FALSE],
                                   rhs = iNull,
                                   coef = stats::setNames(param, names(param.notrans))[index.simplify], ## use transformed values but not transformed name
                                   vcov = vcov.param[index.simplify,index.simplify,drop=FALSE],
                                   df = attr(vcov.param, "df")[names(index.simplify)],
                                   alternative = "two.sided")
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

    transform.sigma <- if(is.na(object$args$transform.sigma)){NULL}else{object$args$transform.sigma}
    transform.k <- if(is.na(object$args$transform.k)){NULL}else{object$args$transform.k}
    transform.rho <- if(is.na(object$args$transform.rho)){NULL}else{object$args$transform.rho}

    ls.lmm <- object$model
    name.lmm <- names(ls.lmm)

    ls.anova <- stats::setNames(lapply(name.lmm, function(iName){ ## iName <- name.lmm[1]
        anova(ls.lmm[[iName]], effects = effects, rhs = rhs, df = df, robust = robust,
              transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, simplify = FALSE)
    }), name.lmm)

    ## ** regenerate a new mlmm object
    out <- do.call("rbind.Wald_lmm",
                   args = c(list(model = ls.anova[[1]], effects = constraint, rhs = rhs, name = names(object$model), sep = ":"), unname(ls.anova[-1]))
                   )
    
    return(out)
    
}

##----------------------------------------------------------------------
### anova.R ends here
