### lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: May 12 2024 (19:18) 
##           By: Brice Ozenne
##     Update #: 3070
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * lmm (documentation)
##' @title Fit Linear Mixed Model
##' @description Fit a linear mixed model defined by a mean and a covariance structure.
##'
##' @param formula [formula] Specify the model for the mean.
##' On the left hand side the outcome and on the right hand side the covariates affecting the mean value.
##' E.g. Y ~ Gender + Gene.
##' @param repetition [formula] Specify the structure of the data: the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' @param structure [character] type of covariance structure, either \code{"CS"} (compound symmetry) or \code{"UN"} (unstructured).
##' @param data [data.frame] dataset (in the long format) containing the observations.
##' @param method.fit [character] Should Restricted Maximum Likelihoood (\code{"REML"}) or Maximum Likelihoood (\code{"ML"}) be used to estimate the model parameters?
##' @param type.information [character] Should the expected information be computed  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param df [logical] Should the degree of freedom be computed using a Satterthwaite approximation?
##' @param weights [formula or character] variable in the dataset used to weight the log-likelihood and its derivative. Should be constant within cluster.
##' @param scale.Omega [formula or character] variable in the dataset used to rescale the residual variance-covariance matrix. Should be constant within cluster.
##' @param trace [interger, >0] Show the progress of the execution of the function.
##' @param control [list] Control values for the optimization method.
##' The element \code{optimizer} indicates which optimizer to use and additional argument will be pass to the optimizer.
##' 
##'
##' @details \bold{Computation time} the \code{lmm} has not been developped to be a fast function as, by default, it uses REML estimation with the observed information matrix and uses a Satterthwaite approximation to compute degrees of freedom (this require to compute the third derivative of the log-likelihood which is done by numerical differentiation). The computation time can be substantially reduced by using ML estimation with the expected information matrix and no calculation of degrees of freedom: arguments \code{method.fit="ML"}, \code{type.information="expected"}, \code{df=FALSE}. This will, however, lead to less accurate p-values and confidence intervals in small samples.
##'
##' By default, the estimation of the model parameters will be made using a Newton Raphson algorithm.
##' This algorithm does not ensure that the residual covariance matrix is positive definite and therefore may sometimes fail.
##' See argument optimizer in \code{\link{LMMstar.options}}.
##'
##' \bold{Argument control:} when using the optimizer \code{"FS"}, the following elements can be used
##' \itemize{
##' \item \code{init}: starting values for the model parameters.
##' \item \code{n.iter}: maximum number of interations of the optimization algorithm.
##' \item \code{tol.score}: score value below which convergence has been reached.
##' \item \code{tol.param}: difference in estimated parameters from two successive iterations below which convergence has been reached.
##' \item \code{trace}: display progress of the optimization procedure.
##' }
##' 
##' \bold{Argument repetition:} when numeric, it will be converted into a factor variable, possibly adding a leading 0 to preserve the ordering.
##' This transformation may cause inconsistency when combining results between different \code{lmm} object. 
##' This is why the grouping variable should preferably be of type character or factor.
##' 
##' @seealso
##' \code{\link{summary.lmm}} for a summary of the model fit. \cr
##' \code{\link{model.tables.lmm}} for a data.frame containing estimates with their uncertainty. \cr
##' \code{\link{plot.lmm}} for a graphical display of the model fit or diagnostic plots. \cr
##' \code{\link{levels.lmm}} to display the reference level. \cr
##' \code{\link{anova.lmm}} for testing linear combinations of coefficients (F-test, multiple Wald tests) \cr
##' \code{\link{effects.lmm}} for evaluating average marginal or counterfactual effects \cr
##' \code{\link{sigma.lmm}} for extracting estimated residual variance-covariance matrices. \cr
##' \code{\link{residuals.lmm}} for extracting residuals or creating residual plots (e.g. qqplots). \cr
##' \code{\link{predict.lmm}} for evaluating mean and variance of the outcome conditional on covariates or other outcome values.

##' @return an object of class \code{lmm} containing the estimated parameter values, the residuals, and relevant derivatives of the likelihood.
##'
##' @keywords models

## * lmm (examples)
##' @examples
##' #### 1- simulate data in the long format ####
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' dL$X1 <- as.factor(dL$X1)
##' dL$X2 <- as.factor(dL$X2)
##' 
##' #### 2- fit Linear Mixed Model ####
##' eCS.lmm <- lmm(Y ~ X1 * X2 + X5, repetition = ~visit|id, structure = "CS", data = dL)
##' 
##' logLik(eCS.lmm) ## -670.9439
##' summary(eCS.lmm)
##'
##'
##' #### 3- estimates ####
##' ## reference level
##' levels(eCS.lmm)$reference
##' ## mean parameters
##' coef(eCS.lmm)
##' model.tables(eCS.lmm)
##' confint(eCS.lmm)
##'
##' ## all parameters
##' coef(eCS.lmm, effects = "all")
##' model.tables(eCS.lmm, effects = "all")
##' confint(eCS.lmm, effects = "all")
##'
##' ## variance-covariance structure
##' sigma(eCS.lmm)
##' 
##' #### 4- diagnostic plots ####
##' quantile(residuals(eCS.lmm))
##' quantile(residuals(eCS.lmm, type = "normalized"))
##'
##' \dontrun{
##' if(require(ggplot2)){
##'   ## investigate misspecification of the mean structure
##'   plot(eCS.lmm, type = "scatterplot")
##'   ## investigate misspecification of the variance structure
##'   plot(eCS.lmm, type = "scatterplot2")
##'   ## investigate misspecification of the correlation structure
##'   plot(eCS.lmm, type = "correlation")
##'   ## investigate misspecification of the residual distribution
##'   plot(eCS.lmm, type = "qqplot")
##' }
##' }
##' 
##' #### 5- statistical inference ####
##' anova(eCS.lmm) ## effect of each variable
##' anova(eCS.lmm, effects = "X11-X21=0") ## test specific coefficient
##' ## test several hypothese with adjustment for multiple comparisons
##' summary(anova(eCS.lmm, effects = c("X11=0","X21=0")))
##'
##' #### 6- prediction ####
##' ## conditional on covariates
##' newd <- dL[1:3,]
##' predict(eCS.lmm, newdata = newd, keep.data = TRUE)
##' ## conditional on covariates and outcome
##' newd <- dL[1:3,]
##' newd$Y[3] <- NA
##' predict(eCS.lmm, newdata = newd, type = "dynamic", keep.data = TRUE)
##'
##' #### EXTRA ####
##' if(require(mvtnorm)){
##' ## model for the average over m replicates
##' ## (only works with independent replicates)
##' Sigma1 <- diag(1,1,1); Sigma5 <- diag(1,5,5)
##' n <- 100
##' dfW <- rbind(data.frame(id = 1:n, rep = 5, Y = rowMeans(rmvnorm(n, sigma = Sigma5))),
##'              data.frame(id = (n+1):(2*n), rep = 1, Y = rmvnorm(n, sigma = Sigma1)))
##' 
##' e.lmmW <- lmm(Y~1, data = dfW, scale.Omega=~rep, control = list(optimizer = "FS"))
##' e.lmm0 <- lmm(Y~1, data = dfW, control = list(optimizer = "FS"))
##' model.tables(e.lmmW, effects = "all")
##' model.tables(e.lmm0, effects = "all")
##' ## TRUE standard error is 1
##'
##' }
##' 

## * lmm (code)
##' @export
lmm <- function(formula, repetition, structure, data,
                weights = NULL, scale.Omega = NULL,
                method.fit = NULL, df = NULL, type.information = NULL, trace = NULL, control = NULL){

    name.out <- c("args", "call", "cluster", "df", "d2Omega", "data", "data.original", 
                  "design", "dOmega", "dVcov", "fitted", "formula", "information", 
                  "logLik", "Omega", "OmegaM1", "opt", "outcome", "param", "reparametrize", 
                  "residuals", "score", "strata", "time", "vcov", "weights", "xfactor"
                  )
    out <- stats::setNames(vector(mode = "list", length = length(name.out)), name.out)

    ## ** 0. extract call and default from package
    out$call <- match.call()
    out$data.original <- data

    options <- LMMstar.options()
    if(is.null(trace)){
        trace <- options$trace
    }else if(trace %in% 0:10 == FALSE){
        stop("Argument \'trace\' should be 0 (no progression displayed in the console) \n",
             "                             1 (display start and end of key steps in the console) \n",
             "                             2 or more (display the execution of the optimization in the console). \n")
    }

    ## handle the case where structure is defined as a function
    if(!missing(structure) && !is.character(structure) && !inherits(structure,"structure")){
        if(inherits(out$call$structure,"name")){ ## instead of a character
            structure <- deparse(out$call$structure)
        }else if(inherits(structure,"function")){ ## instead of a structure object
            structure <- do.call(structure, args = list(~1))
        }
    }
    
    ## ** 1. check and normalize user input
    if(trace>=1){cat("1. Check and normalize user input")}

    ## *** Check all arguments, initialize all arguments but data and structure
    if(trace>=2){cat("\n- normalize argument")}
    outArgs <- .lmmNormalizeArgs(formula = formula, repetition = repetition, structure = structure, data = data,
                                 weights = weights, scale.Omega = scale.Omega,
                                 method.fit = method.fit, df = df, type.information = type.information, trace = trace, control = control,
                                 options = options)    
    var.outcome <- outArgs$var.outcome
    if(trace>=2){cat("\n")}

    ## *** initialize data, e.g. add integer version of cluster/time/strata (XXcluster.indexXX, XXtime.indexXX, ...)
    if(trace>=2){cat("- normalize data")}
    var.all <- c(var.outcome,outArgs$var.X,outArgs$var.cluster,outArgs$var.time,outArgs$var.strata,outArgs$ranef$var,outArgs$var.weights,outArgs$var.scale.Omega)
    if(!missing(structure) && inherits(structure,"structure")){
        var.all <- c(var.all, unlist(structure$name))
    }
    outData <- .lmmNormalizeData(as.data.frame(data)[unique(stats::na.omit(var.all))],
                                 var.outcome = var.outcome, 
                                 var.cluster = outArgs$var.cluster,
                                 var.time = outArgs$var.time,
                                 var.strata = outArgs$var.strata,                         
                                 droplevels = TRUE,
                                 initialize.cluster = outArgs$ranef$crossed,
                                 initialize.time = setdiff(outArgs$ranef$vars, outArgs$var.cluster),
                                 na.rm = TRUE)
    data <- outData$data    
    var.cluster <- outData$var.cluster
    var.time <- outData$var.time
    var.time.original <- attr(var.time,"original")
    var.strata <- outData$var.strata
    var.strata.original <- attr(var.strata,"original")
    index.na <- outData$index.na
    if(trace>=2){cat("\n")}

    ## *** identify cluster/time/strata
    ## cluster
    U.cluster <- levels(data$XXclusterXX) 
    n.cluster <- max(data$XXcluster.indexXX) ## may not match U.cluster in presence of missing values

    ## time
    U.time <- levels(data$XXtimeXX)
    if(length(var.time.original)>=1 && any(!is.na(var.time.original))){
        attr(U.time,"original") <- unique(data[do.call(order,data[var.time.original]),var.time.original,drop=FALSE])
    }
    n.time <- max(data$XXtime.indexXX) ## may not match U.time in presence of missing values

    ## strata
    U.strata <- levels(data$XXstrataXX)
    if(length(var.strata.original)>=1 && any(!is.na(var.strata.original))){
        attr(U.strata,"original") <- unique(data[do.call(order,data[var.strata.original]),var.strata.original,drop=FALSE])
    }
    n.strata <- max(data$XXstrata.indexXX) ## may not match U.strata in presence of missing values

    ## *** normalize cluster/time/strata
    if(trace>=2){cat("- normalize structure")}
    structure <- .lmmNormalizeStructure(structure = structure,
                                        data = data,
                                        ranef = outArgs$ranef,
                                        var.outcome = var.outcome,
                                        var.cluster = var.cluster,
                                        n.cluster = n.cluster,
                                        var.time = var.time,
                                        n.time = n.time,
                                        var.strata = var.strata)

    if(trace>=2){cat("\n")}

    ## *** store results
    out$formula <- list(mean = outArgs$formula,
                        mean.outcome = outArgs$formula.outcome,
                        mean.design = outArgs$formula.design,
                        var = structure$formula$var,
                        cor = structure$formula$cor)
    out$outcome <- list(var = var.outcome)

    dfNNA.Ucluster <- data[!duplicated(data$XXclusterXX),c("XXclusterXX","XXcluster.indexXX")]
    cluster.matchNNA <- match(U.cluster,dfNNA.Ucluster$XXclusterXX)
    out$cluster <- list(n = length(U.cluster), levels = U.cluster, index = dfNNA.Ucluster$XXcluster.indexXX[cluster.matchNNA], var = var.cluster)
    
    dfNNA.Utime <- data[!duplicated(data$XXtimeXX),c("XXtimeXX","XXtime.indexXX")]
    time.matchNNA <- match(U.time,dfNNA.Utime$XXtimeXX)
    out$time <- list(n = length(U.time), levels = U.time, index = dfNNA.Utime$XXtime.indexXX[time.matchNNA], var = var.time)

    dfNNA.Ustrata <- data[!duplicated(data$XXstrataXX),c("XXstrataXX","XXstrata.indexXX")]
    strata.matchNNA <- match(U.strata,dfNNA.Ustrata$XXstrataXX)
    out$strata <- list(n = length(U.strata), levels = U.strata, index = dfNNA.Ustrata$XXstrata.indexXX[strata.matchNNA], var = var.strata)

    out$index.na <- index.na
    out$data <- data ## possible missing data have been removed
    out$weights <-  list(var = c("logLik" = outArgs$var.weights, "Omega" = outArgs$var.scale.Omega))
    out$args <- list(method.fit = outArgs$method.fit,
                     type.information = outArgs$type.information,
                     df = outArgs$df,
                     control =  outArgs$control)
    
    if(trace>=1){cat("\n")}

    ## ** 2. Design matrix and precomputation
    if(trace>=1){cat("2. Design matrix and precomputations \n")}

    ## *** update transformation and precompute moments
    if(structure$class=="CUSTOM"){
        precompute.moments <- FALSE
        if((is.null(structure$d2FCT.sigma) || is.null(structure$d2FCT.rho)) && (outArgs$df || outArgs$method.fit=="REML" || outArgs$type.information=="observed")){
            ## need second derivative but transformation based on dJacobian not implemented!
            options$transform.sigma <- "none"
            options$transform.k <- "none"
            options$transform.rho <- "none"
        }
    }else{
        precompute.moments <- is.na(out$weights$var["Omega"]) && options$precompute.moments        
    }

    ## *** design matrix
    out$design <- .model.matrix.lmm(formula.mean = out$formula$mean.design,
                                    structure = structure,
                                    data = data, var.outcome = out$outcome$var, var.weights = out$weights$var,
                                    precompute.moments = precompute.moments,
                                    drop.X = options$drop.X,
                                    options = options)

    ## *** update xfactor according to factors used in the vcov structure
    ## NOTE: use model.frame to handline splines in the formula
    out$xfactor <- c(list(mean = stats::.getXlevels(stats::terms(out$formula$mean.design),
                                                    stats::model.frame(out$formula$mean.design,data))),
                     out$design$vcov$xfactor)
    out$design$vcov$xfactor <- NULL
    
    if(trace>=1){cat("\n")}

    ## ** 3. Estimate model parameters
    if(trace>=1){cat("3. Estimate model parameters")}

    valid.control <- c("init","n.iter","optimizer","tol.score","tol.param","trace")
    if(any(names(out$args$control) %in% valid.control  == FALSE)){
        stop("Incorrect elements in argument \'control\': \"",paste(names(out$args$control)[names(out$args$control) %in% valid.control  == FALSE], collapse = "\" \""),"\". \n",
             "Valid elements: \"",paste(valid.control, collapse = "\" \""),"\".\n")
    }
    if(identical(out$args$control$init,"lmer")){
        out$args$control$init <- .initializeLMER(formula = out$formula$mean, structure = out$design$vcov, data = data,
                                                 param = out$design$param, method.fit = out$args$method.fit, weights = out$design$weights, scale.Omega = out$design$scale.Omega)
    }else if(inherits(out$design$vcov,"CUSTOM")){
        init.Omega <- .calc_Omega(out$design$vcov, param = c(out$design$vcov$init.sigma,out$design$vcov$init.rho), keep.interim = FALSE)
        out$args$control$init <- init.Omega[[which.max(out$design$vcov$Upattern$n.time)]]
        
    }

    if(trace>0){
        if(out$args$control$trace>0){cat("\n")}
        if(out$args$control$trace>1){cat("\n")}
    }
    outEstimate <- .estimate(design = out$design, time = out$time, method.fit = out$args$method.fit, type.information = out$args$type.information,
                             transform.sigma = options$transform.sigma, transform.k = options$transform.k, transform.rho = options$transform.rho,
                             precompute.moments = precompute.moments, 
                             optimizer = out$args$control$optimizer, init = out$args$control$init, n.iter = out$args$control$n.iter,
                             tol.score = out$args$control$tol.score, tol.param = out$args$control$tol.param, trace = out$args$control$trace)
    param.value <- outEstimate$estimate
    out$opt <- outEstimate[c("cv","n.iter","score","previous.estimate","previous.logLik","control")]
    if((trace==0 && out$args$control$trace>0)){
        cat("\n")
    }
    if(out$opt$cv<=0){
        warning("Convergence issue: no stable solution has been found. \n")
    }
        
    out$param <- param.value

    if(trace>=1){cat("\n")}

    ## ** 4. Compute likelihood derivatives
    if(trace>=1){cat("4. Compute likelihood derivatives \n")}
    outMoments <- .moments.lmm(value = out$param, design = out$design, time = out$time, method.fit = out$args$method.fit, type.information = out$args$type.information,
                               transform.sigma = options$transform.sigma, transform.k = options$transform.k, transform.rho = options$transform.rho,
                               logLik = TRUE, score = TRUE, information = TRUE, vcov = TRUE, df = out$args$df, indiv = FALSE, effects = c("mean","variance","correlation"), robust = FALSE,
                               trace = trace>=2, precompute.moments = precompute.moments, method.numDeriv = options$method.numDeriv, transform.names = FALSE)
    out[names(outMoments)] <- outMoments
    out$fitted <- out$fitted[,1]
    out$residuals <- out$residuals[,1]

    if(trace>=1){cat("\n")}

    ## ** 5. convert to lmm and export
    class(out) <- "lmm"
    return(out)
}

## * .lmmNormalizeArgs 
##' @description Normalize all arguments for lmm
##' @noRd
.lmmNormalizeArgs <- function(formula, repetition, structure, data,
                              weights, scale.Omega,
                              method.fit, df, type.information, trace, control,
                              options){

    ## ** data
    if(!inherits(data,"data.frame")){
        stop("Argument \'data\' should be a \"data.frame\" or inherits from this class. \n")
    }

    ## ** formula
    detail.formula <- formula2var(formula, name.argument = "formula", suggestion = "Something like Y ~ X. \n")
    formula <- detail.formula$formula$regressor ## remove possible random effects

    ## extract variable names
    var.all_meanformula <- detail.formula$vars$all 
    var.X <- detail.formula$vars$regressor
    var.outcome <- detail.formula$vars$response

    ## extract random effects
    if(detail.formula$special=="ranef"){
        if(length(detail.formula$vars$time)>0){
            stop("Incorrect argument \'formula\', \n",
                 "Current version can only handle random intercepts (i.e. no covariates in random effects). \n")
        }
        detail.ranef <- detail.formula$ranef
    }else if(!missing(structure) && inherits(structure,"RE")){
        detail.ranef <- structure$ranef
    }else{
        detail.ranef <- NULL
    }

    ## check
    if(detail.formula$special == "repetition"){
        stop("Random effects in argument \'formula\' should be wrapped into parentheses. \n",
             "Something like Y ~ X1 + (1|id). Otherwise consider using argument \'repetition\'. \n",
             sep = "")
    }
    if(any(var.all_meanformula %in% names(data) == FALSE)){
        invalid <- var.all_meanformula[var.all_meanformula %in% names(data) == FALSE]
        stop("Argument \'formula\' is inconsistent with argument \'data\'. \n",
             "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
             sep = "")
    }
    if(length(var.outcome)!=1){
        stop("Argument \'formula\' must be contain exactly one outcome variable \n")
    }
    if(any(grepl("Intercept",var.X)) ||any(grepl("(Intercept)",var.X))){
        stop("Argument \'formula\' should not contain a variable called \"Intercept\" or \"(Intercept)\". \n")
    }
    if(any(grepl(":",var.X,fixed=TRUE))){
        stop("Argument \'formula\' should not contain a variable whose name contains \":\". \n")
    }
    
    ## ** structure
    missing.structure <- missing(structure) || is.null(structure)
    if(!missing.structure){
        if(inherits(structure,"character")){
            if(detail.formula$special=="ranef" && !identical(structure,"RE")){
                stop("When random effects are specified in the argument \'formula\', argument \'structure\' must be \"RE\". \n")
            }
        }else if(inherits(structure,"structure")){
            if(detail.formula$special=="ranef" && !inherits(structure,"RE")){
                stop("When random effects are specified in the argument \'formula\', argument \'structure\' must be \"RE\". \n")
            }
            if(detail.formula$special=="ranef" && !is.na(structure$name$cor)){
                stop("When random effects are specified in the argument \'formula\', argument \'structure\' should not specify the correlation structure. \n")
            }
            if(detail.formula$special=="ranef" && !is.na(structure$name$strata) && any(detail.formula$ranef$vars %in% structure$name$strata)){
                stop("Strata variable in argument \'structure\' should not correspond to any random effect from argument \'formula\'. \n")
            }
            if(!identical(structure$name$cluster, NA) && structure$name$cluster %in% names(data) == FALSE){
                stop("Cluster variable in structure inconsistent with argument \'data\'. \n",
                     "Variable \"",structure$name$cluster,"\" could not be found in argument \'data\'. \n",
                     sep = "")
            }
            if(!identical(structure$name$time, NA) && structure$name$time %in% names(data) == FALSE){
                stop("Time variable in structure inconsistent with argument \'data\'. \n",
                     "Variable \"",structure$name$time,"\" could not be found in argument \'data\'. \n",
                     sep = "")
            }
            if(!identical(structure$name$strata, NA) && structure$name$strata %in% names(data) == FALSE){
                stop("Strata variable in structure inconsistent with argument \'data\'. \n",
                     "Variable \"",structure$name$strata,"\" could not be found in argument \'data\'. \n",
                     sep = "")
            }
            if(!identical(structure$name$var[[1]], NA) && any(structure$name$var[[1]] %in% names(data) == FALSE)){
                invalid <- structure$name$var[[1]][structure$name$var[[1]] %in% names(data) == FALSE]
                stop("Variance structure inconsistent with argument \'data\'. \n",
                     "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in argument \'data\'. \n",
                     sep = "")
            }
            if(!identical(structure$name$cor[[1]], NA) && any(structure$name$cor[[1]] %in% names(data) == FALSE)){
                invalid <- structure$name$cor[[1]][structure$name$cor[[1]] %in% names(data) == FALSE]
                stop("Correlation structure inconsistent with argument \'data\'. \n",
                     "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in argument \'data\'. \n",
                     sep = "")
            }
            
        }else{
            stop("Argument \'structure\' must be a character or a structure object. \n")
        }
    }

    ## ** repetition
    if(!missing(repetition) && inherits(try(repetition,silent=TRUE),"try-error")){ ## handle typical user mis-communication
        stop("Could not evaluate argument \'repetition\': maybe the symbol \'~\' is missing. \n",
             try(repetition,silent=TRUE))
    }

    missing.repetition <- missing(repetition) || is.null(repetition)
    if(missing.repetition){

        if(detail.formula$special=="ranef"){
            if(detail.ranef$crossed){            
                var.cluster <- NA
                var.time <- detail.ranef$cluster
            }else{
                var.cluster <- detail.ranef$cluster
                var.time <- NA
            }
            if(!missing.structure && inherits(structure,"structure")){
                var.strata <- structure$name$strata
                if(is.na(var.cluster) && !is.na(structure$name$cluster)){
                    var.cluster <- structure$name$cluster
                }else if(!is.na(var.cluster) && !is.na(structure$name$cluster) && !identical(as.character(var.cluster)==as.character(structure$name$cluster))){
                    stop("Inconsistency between the cluster defined via the \'structure\' and the \'formula\' argument. \n",
                         "\"",paste(structure$name$cluster, collapse = "\" \""),"\" vs. \"",paste(var.cluster, collapse = "\" \""),"\" \n")
                }
                if(is.na(var.time) && !is.na(structure$name$time)){
                    var.time <- structure$name$time
                }else if(!is.na(var.time) && !is.na(structure$name$time) && !identical(as.character(var.time)==as.character(structure$name$time))){
                    stop("Inconsistency between the time defined via the \'structure\' and the \'formula\' argument. \n",
                         "\"",paste(structure$name$time, collapse = "\" \""),"\" vs. \"",paste(var.time, collapse = "\" \""),"\" \n")
                }
            }else{
                var.strata <- NA
            }
        }else if(!missing.structure && inherits(structure,"structure")){
            var.cluster <- structure$name$cluster
            var.time <- structure$name$time
            var.strata <- structure$name$strata
        }else{
            var.cluster <- NA
            var.time <- NA
            var.strata <- NA
        }

    }else if(!missing.repetition){
        detail.repetition <- formula2var(repetition, name.argument = "repetition", suggestion = "Something like: ~ time or ~ time|cluster. \n")
        if(detail.repetition$special=="ranef"){
            stop("Incorrect specification of argument \'repetition\'. \n",
                 "Should be something like: ~ time or ~ time|cluster. \n")
        }
        if(any(detail.repetition$vars$all %in% names(data) == FALSE)){
            invalid <- detail.repetition$vars$all[detail.repetition$vars$all %in% names(data) == FALSE]
            if(identical(invalid,"repetition")){
                stop("Argument \'repetition\' is inconsistent with argument \'data\'. \n",
                     "A variable \"repetition\" is used in the \'repetition\' argument, but that may be due to a missing \"=\" sign. \n",
                     sep = "")
            }else{
                stop("Argument \'repetition\' is inconsistent with argument \'data\'. \n",
                     "Variable",if(length(invalid)>1){"s"}," \"",paste(invalid, collapse = "\", \""),"\" could not be found in the dataset. \n",
                     sep = "")
            }
        }
        if(detail.repetition$special == "repetition"){
            var.cluster <- detail.repetition$var$cluster
            var.time <- detail.repetition$var$time
        }else if(detail.repetition$special=="none"){            
            var.cluster <- NA
            var.time <- detail.repetition$var$regressor
        }
        if(!is.null(detail.repetition$var$response)){
            var.strata <- detail.repetition$var$response
        }else{
            var.strata <- NA
        }
        
        ## compatibility structure/repetition
        if(!missing.structure && inherits(structure,"structure")){
            if(!is.na(structure$name$cluster)){
                if(!is.na(var.cluster) && !identical(as.character(structure$name$cluster),as.character(detail.repetition$var$cluster))){
                    stop("Inconsistency between the cluster defined via the \'structure\' and the \'repetition\' argument. \n",
                         "\"",paste(structure$name$cluster, collapse = "\" \""),"\" vs. \"",paste(detail.repetition$var$cluster, collapse = "\" \""),"\" \n")
                }else if(is.na(var.cluster)){
                    var.cluster <- structure$name$cluster
                }
            }
            if(!is.na(structure$name$time)){
                if(!is.na(var.time) && !identical(as.character(structure$name$time),as.character(detail.repetition$var$time))){
                    stop("Inconsistency between the time defined via the \'structure\' and the \'repetition\' argument. \n",
                         "\"",paste(structure$name$time, collapse = "\" \""),"\" vs. \"",paste(detail.repetition$var$time, collapse = "\" \""),"\" \n")
                }else if(is.na(var.time)){
                    var.time <- structure$name$time
                }
            }
            if(!is.na(structure$name$strata)){
                if(!is.na(var.strata) && !identical(as.character(structure$name$strata),as.character(detail.repetition$var$strata))){
                    stop("Inconsistency between the strata defined via the \'structure\' and the \'repetition\' argument. \n",
                         "\"",paste(structure$name$strata, collapse = "\" \""),"\" vs. \"",paste(detail.repetition$var$strata, collapse = "\" \""),"\" \n")
                }else if(is.na(var.strata)){
                    var.strata <- structure$name$strata
                }
            }
            ## Covariates from certain structures contains missing time variables
            test.timecluster <- length(var.time)==1 && !is.na(var.time) && length(var.cluster)==1 && !is.na(var.cluster)
            if(test.timecluster && !inherits(structure,"CUSTOM") && identical(structure$name$var,structure$name$cor) && length(structure$name$var[[1]]==1)){
                if(any(tapply(data[[var.time]],data[[var.cluster]], function(iVec){any(duplicated(iVec))}))){
                    var.time <- c(structure$name$var[[1]], var.time)
                }                
            }
        }
    }

    ## ** weights
    if(!is.null(weights)){
        if(is.character(weights)){
            var.weights <- weights
        }else if(inherits(weights,"formula")){
            var.weights <- all.vars(weights)
        }else {
            stop("Argument \'weights\' should be a character or a formula. \n")
        }
        if(length(var.weights)>1){
            stop("Can only handle a single weights variable. \n")
        }
        if(var.weights %in% names(data)==FALSE){
            stop("Argument \'weights\' is inconsistent with argument \'data\'. \n",
                 "Variable \"",var.weights,"\" could not be found in the dataset. \n",
                 sep = "")
        }
        if(any(data[[var.weights]]<=0)){
            stop("Weights should be strictly positives. \n")
        }
        if(!is.na(var.cluster)){
            if(any(tapply(data[[var.weights]], data[[var.cluster]], function(iW){sum(!duplicated(iW))})>1)){
                stop("Invalid argument \'weights\': weights should be constant within clusters. \n")
            }
        }
    }else{
        var.weights <- NA
    }

    ## ** scale.Omega
    if(!is.null(scale.Omega)){
        if(is.character(scale.Omega)){
            var.scale.Omega <- scale.Omega
        }else if(inherits(scale.Omega,"formula")){
            var.scale.Omega <- all.vars(scale.Omega)
        }else {
            stop("Argument \'scale.Omega\' should be a character or a formula. \n")
        }
        if(length(var.scale.Omega)>1){
            stop("Can only handle a single scale.Omega variable. \n")
        }
        if(var.scale.Omega %in% names(data)==FALSE){
            stop("Argument \'scale.Omega\' is inconsistent with argument \'data\'. \n",
                 "Variable \"",var.scale.Omega,"\" could not be found in the dataset. \n",
                 sep = "")
        }
        if(any(data[[var.scale.Omega]]<=0)){
            stop("Scale.Omega should be strictly positives. \n")
        }
        if(!is.na(var.cluster)){
            if(any(tapply(data[[var.scale.Omega]], data[[var.cluster]], function(iW){sum(!duplicated(iW))})>1)){
                stop("Invalid argument \'scale.Omega\': scale.Omega should be constant within clusters. \n")
            }
        }
    }else{
        var.scale.Omega <- NA
    }

    ## ** method.fit
    if(is.null(method.fit)){
        if(length(var.X)==0 && detail.formula$terms$intercept == FALSE){
            method.fit <- "ML"
        }else{
            method.fit <- options$method.fit
        }
    }else{
        method.fit <- match.arg(method.fit, choices = c("ML","REML"))
        if(length(var.X)==0 && detail.formula$terms$intercept == FALSE && method.fit == "REML"){
            message("Revert back to ML estimation as there is no mean structure. \n")
            method.fit <- "ML"
        }
    }
    
    ## ** degrees of freedom
    if(is.null(df)){
        df <- options$df
    }else if(!is.logical(df)){
        stop("Argument \'df\' should be TRUE or FALSE. \n")
    }
    
    ## ** type of information
    if(is.null(type.information)){
        type.information <- options$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    ## ** control
    if(is.null(control$optimizer)){
        control$optimizer <- options$optimizer
    }else{
        optimx.method <- c("BFGS", "CG", "Nelder-Mead", "nlminb", "bobyqa")
        control$optimizer <- match.arg(control$optimizer, c("FS",optimx.method)) ## FS = fisher scoring
    }
    if(is.null(control$trace)){
        control$trace <- trace-1
    }else{
    }

    ## ** export
    return(list(formula = formula, formula.outcome = stats::update(formula,.~1), formula.design = detail.formula$formula$design, ranef = detail.ranef,
                var.outcome = var.outcome, var.X = var.X, var.cluster = var.cluster, var.time = var.time, var.strata = var.strata,
                var.weights = var.weights, var.scale.Omega = var.scale.Omega,
                method.fit = method.fit,
                df = df,
                type.information = type.information,
                control = control))

}

## * .lmmNormalizeData
##' @description Normalize data argument for lmm
##'
##' @param droplevels [logical] should un-used levels in cluster/time/strata be removed?
##' @param initalize.cluster [logical] when the cluster variable is NA/NULL,
##' should a different cluster per observation be used (0) or the same cluster for all observations (1)
##' @param initalize.time [logical] when the time variable is NA/NULL, ##' how should the observations be ordered
##' @param na.rm [logical] should row containing missing values be removed from the design matrix?
##'
##' @details Argument \code{var.outcome} is NA when the function is called by \code{model.matrix.lmm}
##' 
##' @noRd
.lmmNormalizeData <- function(data,
                              var.outcome, 
                              var.cluster, var.time, var.strata, droplevels,
                              initialize.cluster, initialize.time, na.rm){


    ## ** normalize
    if(is.null(var.outcome)){var.outcome <- NA}
    if(is.null(var.cluster)){var.cluster <- NA}
    if(is.null(var.time)){var.time <- NA}
    if(is.null(var.strata)){var.strata <- NA}

    out <- list(var.cluster = var.cluster,
                var.time = var.time,
                var.strata = var.strata)

    if(is.list(droplevels)){
        refLevel.strata <- droplevels$strata
        refLevel.time <- droplevels$time
        droplevels <- TRUE
    }else{
        refLevel.strata <- NULL
        refLevel.time <- NULL
    }

    ## ** test
    names.data <- names(data)
    if(NROW(data)==0){
        stop("Argument \'data\' has 0 rows. \n")
    }
    if("XXindexXX" %in% names.data){
        stop("Argument \'data\' should not contain a column named \"XXindexXX\" as this name is used internally by the lmm function. \n")
    }
    if("XXtimeXX" %in% names.data){
        stop("Argument \'data\' should not contain a column named \"XXtimeXX\" as this name is used internally by the lmm function. \n")
    }
    if("XXtime.indexXX" %in% names.data){
        stop("Argument \'data\' should not contain a column named \"XXtime.indexXX\" as this name is used internally by the lmm function. \n")
    }
    if("XXclusterXX" %in% names.data){
        stop("Argument \'data\' should not contain a column named \"XXclusterXX\" as this name is used internally by the lmm function. \n")
    }
    if("XXcluster.indexXX" %in% names.data){
        stop("Argument \'data\' should not contain a column named \"XXcluster.indexXX\" as this name is used internally by the lmm function. \n")
    }
    if("XXstrataXX" %in% names.data){
        stop("Argument \'data\' should not contain a column named \"XXstrataXX\" as this name is used internally by the lmm function. \n")
    }
    if("XXstrata.indexXX" %in% names.data){
        stop("Argument \'data\' should not contain a column named \"XXstrata.indexXX\" as this name is used internally by the lmm function. \n")
    }
 
    if(any(!is.na(var.cluster)) && any(var.cluster %in% names.data == FALSE)){
        stop("Argument \'repetition\' is inconsistent with argument \'data\'. \n",
             "Could not find column \"",paste(var.cluster, collapse = "\" \""),"\" indicating the cluster in argument \'data\'. \n", sep="")
    }
    if(any(!is.na(var.time)) && any(var.time %in% names.data == FALSE)){
        stop("Argument \'repetition\' is inconsistent with argument \'data\'. \n",
             "Could not find column \"",paste(var.time, collapse = "\" \""),"\" indicating the time in argument \'data\'. \n", sep="")
    }
    if(any(!is.na(var.strata)) && any(var.strata %in% names.data == FALSE)){
        stop("Argument \'structure\' is inconsistent with argument \'data\'. \n",
             "Could not find column \"",paste(var.strata, collapse = "\" \""),"\" indicating the strata in argument \'data\'. \n",sep="")
    }
    if(length(var.cluster)>1){
        stop("Incorrect specification of argument \'repetition\': too many cluster variables. \n",
             "There should be exactly one variable after the grouping symbol (|), something like: ~ time|cluster or strata ~ time|cluster. \n", sep = "")
    }

    ## ** convert logical into factor
    test.logical <- sapply(data,is.logical)
    if(any(test.logical)){ ## avoid an error when computing partial residuals since the formula interface treat logical as factor 
        data[test.logical] <- lapply(data[test.logical], as.factor)
    }
    
    ## ** index
    data$XXindexXX <- 1:NROW(data)

    ## ** outcome
    ## handle a possible transformation of the outcome
    ## (model.frame excludes NAs - this is avoided by removing NAs before hand)
    ## if(!is.na(var.outcome)){
    ##     indexY.NNA <- which(!is.na(data[[var.outcome]]))
    ##     data[indexY.NNA,var.outcome] <- stats::model.response(stats::model.frame(formula.outcome, data[indexY.NNA,,drop=FALSE]))
    ## }

    ## ** cluster
    if(is.na(var.cluster)){
        if(!is.null(initialize.cluster) && initialize.cluster==1){
            if(any(duplicated(data[,var.time]))){
                stop("Incorrect specification of argument \'repetition\': missing cluster variable. \n",
                     "Cannot guess the cluster variable with non-unique levels for the cross random effects. \n")
            }
            data$XXclusterXX <- as.factor("1")
        }else{
            data$XXclusterXX <- addLeading0(1:NROW(data), as.factor = TRUE, code = attr(var.cluster,"code"))
        }

        out$var.cluster <- "XXclusterXX"
        attr(out$var.cluster, "original") <- NA
    }else{
        if(is.factor(data[[var.cluster]])){
            if(droplevels){
                data$XXclusterXX <- droplevels(data[[var.cluster]])
            }else{
                data$XXclusterXX <- data[[var.cluster]]
            }
        }else if(is.numeric(data[[var.cluster]]) && all(data[[var.cluster]] %% 1 == 0)){
            data$XXclusterXX <- addLeading0(data[[var.cluster]], as.factor = TRUE, code = attr(var.cluster,"code"))
        }else{
            data$XXclusterXX <- factor(data[[var.cluster]], levels = sort(unique(data[[var.cluster]])))
        }
        attr(out$var.cluster, "original") <- var.cluster
    }

    ## ** time
    if(length(var.time)>1){
        ## create new variable summarizing all variables
        data$XXtimeXX <- nlme::collapse(lapply(var.time, function(iX){as.factor(data[[iX]])}), as.factor = is.null(refLevel.time))
        out$var.time <- "XXtimeXX"
        attr(out$var.time, "original") <- var.time     
    }else if(is.na(var.time)){
        if(length(initialize.time)==0){
            iSplit <- do.call(rbind,by(data, data$XXclusterXX, FUN = function(iDF){cbind(XXindexXX = iDF$XXindexXX, XXtimeXX = 1:NROW(iDF))}))
        }else{ ## re-order dataset according to variables defining the variance-covariance pattern to minimize the number of patterns
            iSplit <- do.call(rbind,by(data, data$XXclusterXX, FUN = function(iDF){ ## iDF <- data[data$XXclusterXX==data$XXclusterXX[1],]
                iOrder <- do.call(order,iDF[initialize.time])
                cbind(XXindexXX = iDF[iOrder,"XXindexXX"], XXtimeXX = 1:NROW(iDF))
            }))
        }
        data[iSplit[,"XXindexXX"],"XXtimeXX"] <- addLeading0(iSplit[,"XXtimeXX"], as.factor = TRUE, code = attr(var.time,"code"))
        
        out$var.time <- "XXtimeXX"
        attr(out$var.time, "original") <- NA
    }else{
        if(is.factor(data[[var.time]])){            
            if(droplevels){
                data$XXtimeXX <- droplevels(data[[var.time]])
            }else{
                data$XXtimeXX <- data[[var.time]]
            }
        }else if(is.numeric(data[[var.time]])){
            if(all(data[[var.time]]>=0) && all(data[[var.time]] %% 1 == 0)){
                data$XXtimeXX <- addLeading0(data[[var.time]], as.factor = TRUE, code = attr(var.time,"code"))
            }else{
                data$XXtimeXX <- factor(data[[var.time]], levels = unique(sort(data[[var.time]])))
            }
        }else{
            data$XXtimeXX <- factor(data[[var.time]], levels = sort(unique(data[[var.time]])))
        }
        attr(out$var.time, "original") <- var.time
    }
    if(!is.null(refLevel.time)){  ## from model.matrix to restaure levels used when fitting the lmm
        data$XXtimeXX <- factor(data$XXtimeXX, refLevel.time)
    }

    ## ** strata
    ## NOTE: .formulaStructure makes sure that there is at most 1 strata variable
    if(is.na(var.strata)){
        data$XXstrataXX <- factor(1)

        out$var.strata <- "XXstrata.indexXX"
        attr(out$var.strata, "original") <- NA
    }else{
        if(is.factor(data[[var.strata]])){
            if(droplevels){
                data$XXstrataXX <- droplevels(data[[var.strata]])
            }else{
                data$XXstrataXX <- data[[var.strata]]
            }
        }else if(is.numeric(data[[var.strata]]) && all(data[[var.strata]] %% 1 == 0)){
            data$XXstrataXX <- addLeading0(data[[var.strata]], as.factor = TRUE, code = attr(var.strata,"code"))
        }else{
            data$XXstrataXX <- factor(data[[var.strata]], levels = sort(unique(data[[var.strata]])))
        }
        attr(out$var.strata, "original") <- var.strata
    }
    if(!is.null(refLevel.strata)){ ## from model.matrix to restaure levels used when fitting the lmm
        data$XXstrataXX <- factor(data$XXstrataXX, refLevel.strata)
    }

    ## ** indicators
    data$XXcluster.indexXX <- as.numeric(data$XXclusterXX)
    data$XXtime.indexXX <- as.numeric(data$XXtimeXX)
    data$XXstrata.indexXX <- as.numeric(data$XXstrataXX)

    ## ** check distinct time and unique strata values within clusters
    n.cluster <- max(data$XXcluster.indexXX)
    n.time <- max(data$XXtime.indexXX)
    n.strata <- max(data$XXstrata.indexXX)
    if(n.cluster >= n.time){
        test.duplicated <- tapply(data$XXcluster.indexXX, data$XXtime.indexXX, function(iT){any(duplicated(iT))})
    }else{
        test.duplicated <- tapply(data$XXtime.indexXX, data$XXcluster.indexXX, function(iT){any(duplicated(iT))})
    }    
    if(any(test.duplicated)){
        stop("Argument \'repetition\' is inconsistent with argument \'data\'. \n",
             "The time variable should contain unique values within clusters \n", sep="")
    }

    if(n.strata>1 & n.time >1){
        test.sameStrata <- tapply(data$XXstrata.indexXX, data$XXcluster.indexXX, function(iT){any(iT[1]!=iT[-1])})
        if(any(test.sameStrata)){
            stop("Argument \'repetition\' is inconsistent with argument \'data\'. \n",
                 "When a variable is used to stratify the variance structure, all observations within a cluster must belong to same strata. \n")
        
        }
    }

    ## ** missing data
    index.na <- which(rowSums(is.na(data[names.data]))>0)
    if(na.rm && length(index.na) == NROW(data)){
        var.na <- names.data[colSums(!is.na(data))==0]
        if(length(var.na)==0){
            stop("All observations have at least one missing data. \n")
        }else if(length(var.na)==1){
            stop("Variable \"",var.na,"\" contains only missing data. \n")
        }else{
            stop("Variables \"",paste(var.na, collapse="\" \""),"\" contain only missing data. \n")
        }
    }else if(length(index.na)>0){

        test.naOutcome <- is.na(data[[var.outcome]])
        test.naOther <- rowSums(is.na(data[setdiff(names.data,var.outcome)]))>0
        warning <- FALSE
        text.warning <- NULL
        
        if( any(test.naOther > test.naOutcome) ){
            index.row <- which(test.naOther > test.naOutcome)
            test.naOther2 <- colSums(is.na(data[index.row,setdiff(names.data,var.outcome),drop=FALSE]))
            name.naOther2 <- names(test.naOther2)[test.naOther2>0]

            warning <- TRUE
            nobs.naOther <- sum(test.naOther)
            text.warning <- c(text.warning,
                              paste0("Can only handle missing values in the outcome variable ",var.outcome,". \n",
                                     "  ",nobs.naOther," observation",if(nobs.naOther>1){"s"}," with missing values in \"",paste(name.naOther2, collapse="\" \""),"\" ",ifelse(nobs.naOther==1,"has","have")," been removed. \n", sep = ""))
        }

        attr(index.na, "cluster") <- data[index.na,"XXclusterXX"]
        attr(index.na, "time") <- data[index.na,"XXtimeXX"]

        if(na.rm){
            keep <- list(nlevel.cluster = max(data$XXcluster.indexXX),
                         nlevel.time = max(data$XXtime.indexXX),
                         nlevel.strata = max(data$XXstrata.indexXX))
            data <- data[-index.na,, drop=FALSE]
            data$XXcluster.indexXX <- as.numeric(droplevels(data$XXclusterXX))
            if(droplevels){
                data$XXtime.indexXX <- as.numeric(droplevels(data$XXtimeXX))
                data$XXstrata.indexXX <- as.numeric(droplevels(data$XXstrataXX))
            }

            loss.cluster <- keep$nlevel.cluster - max(data$XXcluster.indexXX)
            loss.time <- keep$nlevel.time - max(data$XXtime.indexXX) 
            loss.strata <- keep$nlevel.strata - max(data$XXstrata.indexXX)
            if(loss.cluster>0){
                warning <- TRUE
                text.warning <- c(text.warning,paste0("  ",loss.cluster," cluster",ifelse(loss.cluster==1," has","s have")," been removed. \n"))
            }
            if(loss.time>0){
                warning <- TRUE
                text.warning <- c(text.warning,paste0("  ",loss.time," timepoint",ifelse(loss.time==1," has","s have")," been removed. \n"))
            }
            if(loss.strata>0){
                warning <- TRUE
                text.warning <- c(text.warning,paste0("  ",loss.strata," strata",ifelse(loss.strata==1," has"," have")," been removed. \n"))
            }
            if(warning){
                warning(text.warning)
            }
        }
    }else{
        index.na <- NULL
    }

    ## ** export
    out$data <- data
    out$index.na <- index.na
    return(out)
}

## * .lmmNormalizeStructure
##' @description Normalize structure argument for lmm
##' @noRd
.lmmNormalizeStructure <- function(structure, data, ranef,
                                   var.outcome,
                                   var.cluster, n.cluster,
                                   var.time, n.time,
                                   var.strata){

    var.cluster.original <- attr(var.cluster,"original")
    var.time.original <- attr(var.time,"original")
    var.strata.original <- attr(var.strata,"original")
    structure1time2ID <- c("CS","RE","TOEPLITZ","UN","EXP") ## simplify structures to ID when a single timepoint
    
    ## ** initialize structure when not specified
    if(missing(structure) || is.null(structure)){
        if(!is.null(ranef$formula)){
            structure <- "RE"
        }else if(any(duplicated(data[["XXclusterXX"]]))){
            if(all(is.na(var.time.original))){
                structure <- "CS"
            }else{
                structure <- "UN"
            }
        }else if(n.time>1){
            structure <- "IND"
        }else{
            structure <- "ID"
        }
    }            

    ## ** check that time and cluster are appropriately defined with respect to structure
    if(is.character(structure)){
        if(all(is.na(var.cluster.original)) && structure %in% c("CS","TOEPLITZ","UN")){
            stop("Incorrect specification of argument \'repetition\': missing cluster variable. \n",
                 "Should have exactly one variable after the grouping symbol (|), something like: ~ time|cluster or strata ~ time|cluster. \n")
        }
        if(all(is.na(var.time.original)) && structure %in% c("IND","TOEPLITZ","UN")){
            stop("Incorrect specification of argument \'repetition\': missing time variable. \n",
                 "Should have exactly one variable before the grouping symbol (|), something like: ~ time|cluster or strata ~ time|cluster. \n")
        }
    }else if(inherits(structure,"structure")){
        if(all(is.na(var.cluster.original)) && structure$class %in% c("CS","TOEPLITZ","UN")){
            stop("Incorrect specification of argument \'repetition\': missing cluster variable. \n",
                 "Should have exactly one variable after the grouping symbol (|), something like: ~ time|cluster or strata ~ time|cluster. \n")
        }
        if(all(is.na(var.time.original)) && structure$class %in% c("IND","TOEPLITZ","UN") && all(is.na(structure$name$var[[1]]))){
            stop("Incorrect specification of argument \'repetition\': missing time variable. \n",
                 "Should have exactly one variable before the grouping symbol (|), something like: ~ time|cluster or strata ~ time|cluster. \n")
        }
    }

    ## ** update structure with cluster/time/strata variables
    var.clusterS <- "XXcluster.indexXX"
    attr(var.clusterS,"original") <- var.cluster.original

    var.timeS <- "XXtime.indexXX"
    attr(var.timeS,"original") <- var.time.original

    ## random effect not defined in formula or structure but can be guessed from repetition
    if((inherits(structure,"RE") && is.null(structure$ranef) || identical(structure,"RE")) && is.null(ranef) && !is.null(var.cluster.original)){
        ranef <- formula2var(stats::as.formula(paste0("~(1|",var.cluster.original,")")))$ranef
    }
    ## special case where no formula in structure and only type as an argument of structure
    if(inherits(structure,"structure")){
        if(is.null(structure$formula$var) && is.null(structure$formula$cor) && identical("type", setdiff(names(structure$call),""))){
            type <- structure$type
            structure <- structure$class
        }else{
            type <- NULL
        }
    }else{
        type <- NULL
    }

    if(inherits(structure,"structure")){

        ## exclude useless strata variables
        if(any(!is.na(var.strata.original))){
            test.uniqueStrata <- sapply(var.strata.original,function(iName){length(unique(data[[iName]]))})
            var.strata.original <- var.strata.original[test.uniqueStrata>1]
            if(all(test.uniqueStrata==1)){
                message("Single strata: move from stratified to non-stratified ",structure$class," structure. \n")
            }else if(any(test.uniqueStrata==1)){
                message("Remove strata variable \"",paste(var.strata.original[test.uniqueStrata==1], collapse = "\", \""),"\" as it takes a single distinct value. \n")
            }
        }
        if(n.time==1 && (structure$class %in% structure1time2ID || (structure$class == "IND" && (all(is.na(structure$name$var)) || all(structure$name$var[[1]] %in% var.time.original))))){
            message("Single timepoint: move from ",structure$class," to ID structure. \n")
            structure$call[[1]] <- parse(text="ID")
        }

        if(inherits(structure,"RE")){
            ## special case with random effects
            structure <- stats::update(structure, var.cluster = var.clusterS, var.time = var.timeS, var.strata = var.strata.original, ranef = ranef)
        }else{
            structure <- stats::update(structure, var.cluster = var.clusterS, var.time = var.timeS, var.strata = var.strata.original, n.time = n.time)
        }

    }else if(is.character(structure)){
        args.structure <- list(var.cluster = var.clusterS,
                               var.time = var.timeS)

        if(!is.null(type)){
            args.structure$type <- type
        }

        if(is.na(var.strata.original)){
            args.structure$formula <- ~1
        }else{ ## exclude useless strata variables
            test.uniqueStrata <- sapply(var.strata.original,function(iName){length(unique(data[[iName]]))})
            if(all(test.uniqueStrata==1)){
                message("Single strata: move from stratified to non-stratified ",structure," structure. \n")
                args.structure$formula <- ~1
            }else if(any(test.uniqueStrata==1)){
                message("Remove strata variable \"",paste(var.strata.original[test.uniqueStrata==1], collapse = "\", \""),"\" as it takes a single distinct value. \n")
                args.structure$formula <- stats::as.formula(paste(var.strata.original[test.uniqueStrata!=1],"~1"))
            }else{
                args.structure$formula <- stats::as.formula(paste(var.strata.original,"~1"))
            }
        }
        if(n.time==1 && (structure %in% structure1time2ID || (structure == "IND" && all(all.vars(args.structure$formula) %in% var.time.original)))){
            message("Single timepoint: move from ",structure," to ID structure. \n")
            structure <- "ID"
        }else if(structure == "UN" && all(table(data[[var.cluster]])<2)){
            message("Single repetition per cluster: move from ",structure," to IND structure. \n")
            structure <- "IND"
        }

        if(structure %in% c("UN","EXP","TOEPLITZ") || (n.time>1 && structure == "IND")){
            args.structure$add.time <- var.time.original
        }else if(structure %in% "RE"){
            args.structure$ranef <- ranef
        }
        structure <- do.call(structure, args = args.structure)    
    }

    ## ** Sanity checks

    ## *** Variability within cluster (for clusters with more than a single value)
    if(!is.null(structure$formula$cor)){
        if(all(table(data$XXcluster.indexXX)<2)){
            warning("Single observation per cluster: will not be able to estimate correlation parameters. \n")
        }else{
            check <- FALSE
            for(iC in 1:n.cluster){ ## iC <- 1
                iData <- data[data$XXcluster.indexXX==iC,var.outcome]
                if(length(iData)==1 || all(is.na(iData))){
                    next
                }else if(max(iData, na.rm = TRUE)-min(iData, na.rm = TRUE)>1e-12){
                    check <- TRUE
                    break
                }
            }
            if(check == FALSE){
                warning("Constant outcome values within cluster. \n")
            }
        }
    }

    ## ** export
    return(structure)
}



######################################################################
### lmm.R ends here
