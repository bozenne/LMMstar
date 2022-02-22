### lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: feb 22 2022 (17:34) 
##           By: Brice Ozenne
##     Update #: 1461
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
##'g
##' @param formula [formula] Specify the model for the mean.
##' On the left hand side the outcome and on the right hand side the covariates affecting the mean value.
##' E.g. Y ~ Gender + Gene.
##' @param repetition [formula] Specify the model for the covariance.
##' On the right hand side the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' On the left hand side, a possible stratification variable, e.g. group ~ time|id. In that case the mean structure should only be stratified on this variable using interactions.
##' @param structure [character] type of covariance structure, either \code{"CS"} (compound symmetry) or \code{"UN"} (unstructured).
##' @param data [data.frame] dataset (in the long format) containing the observations.
##' @param method.fit [character] Should Restricted Maximum Likelihoood (\code{"REML"}) or Maximum Likelihoood (\code{"ML"}) be used to estimate the model parameters?
##' @param type.information [character] Should the expected information be computed  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param df [logical] Should the degree of freedom be computed using a Satterthwaite approximation?
##' @param weights [formula or character] variable in the dataset used to weight the log-likelihood and its derivative. Should be constant within cluster.
##' @param scale.Omega [formula or character] variable in the dataset used to rescale the residual variance-covariance matrix. Should be constant within cluster.
##' @param trace [interger, >0] Show the progress of the execution of the function.
##' @param control [list] Control values for the optimization method. The element \code{optimizer} indicates which optimizer to use and additional argument will be pass to the optimizer.
##'
##' @details \bold{Computation time} the \code{lmm} has not been developped to be a fast function as, by default, it uses REML estimation with the observed information matrix and uses a Satterthwaite approximation to compute degrees of freedom (this require to compute the third derivative of the log-likelihood which is done by numerical differentiation). The computation time can be substantially reduced by using ML estimation with the expected information matrix and no calculation of degrees of freedom: arguments \code{method.fit="ML"}, \code{type.information="expected"}, \code{df=FALSE}. This will, however, lead to less accurate p-values and confidence intervals in small samples.
##'
##' By default, the estimation of the model parameters will be made using the \code{nlme::gls} function.
##' See argument optimizer in \code{\link{LMMstar.options}}
##'
##' @seealso
##' \code{\link{summary.lmm}} for a summary of the model fit. \cr
##' \code{\link{model.tables.lmm}} for a data.frame containing estimates with their uncertainty. \cr
##' \code{\link{plot.lmm}} for a graphical display of the model fit or diagnostic plots. \cr
##' \code{\link{levels.lmm}} to display the reference level. \cr
##' \code{\link{anova.lmm}} for testing linear combinations of coefficients (F-test, multiple Wald tests) \cr
##' \code{\link{getVarCov.lmm}} for extracting estimated residual variance-covariance matrices. \cr
##' \code{\link{residuals.lmm}} for extracting residuals or creating residual plots (e.g. qqplots). \cr
##' \code{\link{predict.lmm}} for evaluating mean and variance of the outcome conditional on covariates or other outcome values.

##' @return an object of class \code{lmm} containing the estimated parameter values, the residuals, and relevant derivatives of the likelihood.

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
##' if(require(emmeans)){
##'   dummy.coef(eCS.lmm)
##' }
##' 
##' ## all parameters
##' coef(eCS.lmm, effects = "all")
##' model.tables(eCS.lmm, effects = "all")
##' confint(eCS.lmm, effects = "all")
##'
##' ## variance-covariance structure
##' getVarCov(eCS.lmm)
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
##' predict(eCS.lmm, newdata = newd, keep.newdata = TRUE)
##' ## conditional on covariates and outcome
##' newd <- dL[1:3,]
##' newd$Y[3] <- NA
##' predict(eCS.lmm, newdata = newd, type = "dynamic", keep.newdata = TRUE)
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
##' e.lmm0 <- lmm(Y~1, data = dfW)
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

    out <- list(call = match.call(), data.original = data)
    options <- LMMstar.options()

    ## ** check and normalize user imput
    if(is.null(trace)){
        trace <- options$trace
    }
    if(trace>=1){cat("1. Check and normalize user imput \n")}

    ## *** data
    if(!inherits(data,"data.frame")){
        stop("Argument \'data\' should be a data.frame or inherit from it. \n")
    }
    if("XXindexXX" %in% names(data)){
        stop("Argument \'data\' should not contain a column named \"XXindexXX\" as this name is used internally by the lmm function. \n")
    }
    if("XXtimeXX" %in% names(data)){
        stop("Argument \'data\' should not contain a column named \"XXtimeXX\" as this name is used internally by the lmm function. \n")
    }
    if("XXtime.indexXX" %in% names(data)){
        stop("Argument \'data\' should not contain a column named \"XXtime.indexXX\" as this name is used internally by the lmm function. \n")
    }
    if("XXclusterXX" %in% names(data)){
        stop("Argument \'data\' should not contain a column named \"XXclusterXX\" as this name is used internally by the lmm function. \n")
    }
    if("XXstrataXX" %in% names(data)){
        stop("Argument \'data\' should not contain a column named \"XXstrataXX\" as this name is used internally by the lmm function. \n")
    }
    if("XXstrata.indexXX" %in% names(data)){
        stop("Argument \'data\' should not contain a column named \"XXstrata.indexXX\" as this name is used internally by the lmm function. \n")
    }
    
    ## *** objective function
    if(is.null(method.fit)){
        method.fit <- options$method.fit
    }else{
        method.fit <- match.arg(method.fit, choices = c("ML","REML"))
    }
    out$method.fit <- method.fit

    ## *** degrees of freedom
    if(is.null(df)){
        df <- options$df
    }else if(!is.logical(df)){
        stop("Argument \'df\' should be TRUE or FALSE. \n")
    }
    
    ## *** type of information
    if(is.null(type.information)){
        type.information <- options$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }
    
    ## *** optimizer
    if(is.null(control$optimizer)){
        optimizer <- options$optimizer
    }else{
        optimx.method <- c("BFGS", "CG", "Nelder-Mead", "nlminb", "bobyqa")
        optimizer <- match.arg(control$optimizer, c("gls","FS",optimx.method)) ## FS = fisher scoring
        control$optimizer <- NULL
    }

    ## *** repetition 
    if(trace>=2){cat("- repetition ")}
    if(missing(repetition)){
        if(missing(structure) || identical(structure,"ID")){
            var.cluster  <- NA
            var.time  <- NA
            var.strata  <- NA
            structure <- "ID"
        }else if(inherits(structure,"structure")){
            var.cluster <- structure$name$cluster
            var.time <- structure$name$time
            var.strata <- structure$name$strata
            if(is.na(var.cluster)){
                stop("Missing argument \'repetition\': should be a formula like ~ time|cluster or strata ~ time|cluster.")
            }else if(var.cluster %in% names(data) == FALSE){
                stop("Argument \'structure\' is inconsistent with argument \'data\'. \n",
                     "Could not find column \"",var.cluster,"\" indicating the cluster in argument \'data\'. \n")
            }
            if(is.na(var.time)){
                stop("Missing argument \'repetition\': should be a formula like ~ time|cluster or strata ~ time|cluster.")
            }else if(var.time %in% names(data) == FALSE){
                stop("Argument \'structure\' is inconsistent with argument \'data\'. \n",
                     "Could not find column \"",var.time,"\" indicating the cluster in argument \'data\'. \n")
            }
            if(!is.na(var.strata) && var.strata %in% names(data) == FALSE){
                stop("Argument \'structure\' is inconsistent with argument \'data\'. \n",
                     "Could not find column \"",var.strata,"\" indicating the strata in argument \'data\'. \n")
            }
        }else{
            stop("Missing argument \'repetition\': should be a formula like ~ time|cluster or strata ~ time|cluster.")
        }
    }else{
        if(!inherits(repetition,"formula")){
            stop("Argument \'repetition\' must be of class formula, something like: ~ time|cluster or strata ~ time|cluster. \n")
        }
        res.split <- strsplit(deparse(repetition),"|", fixed = TRUE)[[1]]
        if(length(res.split)>2){
            stop("Incorrect specification of argument \'repetition\'. \n",
                 "The symbol | should only exacly once, something like: ~ time|cluster or strata ~ time|cluster. \n")
        }
        if(length(res.split)==1){
            var.cluster <- NA
        }else{
            var.cluster <- trimws(res.split[2], which = "both")
            if(length(var.cluster)==0){var.cluster <- NA}
        }
        var.time <- all.vars(stats::update(stats::as.formula(res.split[1]),0~.))
        if(length(var.time)==0){var.time <- NA}
        var.strata <- setdiff(all.vars(stats::as.formula(res.split[1])), var.time)
        if(length(var.strata)==0){var.strata <- NA}

        if(length(var.cluster)>1){
            stop("Incorrect specification of argument \'repetition\': too many cluster variables. \n",
                 "There should be exactly one variable after the grouping symbol (|), something like: ~ time|cluster or strata ~ time|cluster. \n")
        }else if(length(var.cluster)==1 && !is.na(var.cluster) && var.cluster %in% names(data) == FALSE){
            stop("Argument \'repetition\' is inconsistent with argument \'data\'. \n",
                 "Could not find column \"",var.cluster,"\" indicating the cluster in argument \'data\' \n")
        }
        if(length(var.time)>1){
            stop("Incorrect specification of argument \'repetition\': too many time variables. \n",
                 "There should be exactly one variable before the grouping symbol (|), something like: ~ time|cluster or strata ~ time|cluster. \n")
        }else if(length(var.time)==1 && !is.na(var.time) && var.time %in% names(data) == FALSE){
            stop("Argument \'repetition\' is inconsistent with argument \'data\'. \n",
                 "Could not find column \"",var.time,"\" indicating the cluster in argument \'data\' \n")
        }
        if(length(var.strata)>1){
            stop("Incorrect specification of argument \'repetition\': too many strata variables. \n",
                 "There should be at most one variable on the left hand side, something like: ~ time|cluster or strata ~ time|cluster. \n")
        }else if(length(var.strata)==1 && !is.na(var.strata) && var.strata %in% names(data) == FALSE){
            stop("Argument \'repetition\' is inconsistent with argument \'data\'. \n",
                 "Could not find column \"",var.strata,"\" indicating the strata in argument \'data\' \n")
        }
        if(length(var.time)==1 && is.na(var.time)){ ## add time when only one obs per cluster
            if(length(var.cluster)==1 && is.na(var.cluster)){
                stop("Incorrect specification of argument \'repetition\': missing time and cluster variable. \n",
                     "Should have exactly one variable after the grouping symbol (|), something like: ~ time|cluster or strata ~ time|cluster. \n")
            }else if(any(duplicated(data[[var.cluster]]))){
                stop("Incorrect specification of argument \'repetition\': missing time variable. \n",
                     "There should be exactly one variable before the grouping symbol (|), something like: ~ time|cluster or strata ~ time|cluster. \n")
            }
        }
    }

    if(trace>=2){cat("\n")}

    ## *** residual variance-covariance structure
    if(trace>=2){cat("- residual variance-covariance structure  ")}

    if(missing(structure)){
        if(!is.na(var.cluster) && any(duplicated(data[[var.cluster]]))){
            structure <- "UN"
        }else if(!is.na(var.time) && sum(!duplicated(data[[var.time]]))>1){
            structure <- "IND"
        }else{
            structure <- "ID"
        }
    }
    
    if(inherits(structure,"structure")){
        if(optimizer=="gls"){
            stop("Argument \'structure\' must be a character when using gls optimizer. \n",
                 "Otherwise add argument control = list(optimizer = \"FS\"). \n")
        }
        if(!missing(repetition)){
            if(!is.na(structure$name$time) && structure$name$time != var.time){
                stop("Argument \'structure\' is inconsistent with argument \'repetition\'. \n",
                     "Not the same time variable: ",structure$name$time," vs. ",var.time,".\n")
            }
            structure$name$time <- "XXtimeXX"

            if(!is.na(structure$name$cluster) && structure$name$cluster != var.cluster){
                stop("Argument \'structure\' is inconsistent with argument \'repetition\'. \n",
                     "Not the same cluster variable: ",structure$name$cluster," vs. ",var.cluster,".\n")
            }
            structure$name$cluster <- "XXclusterXX"
        }
        type.structure <- structure$type
    }else if(inherits(structure,"character")){
        type.structure <- match.arg(structure, c("ID","IND","CS","UN"))
        if(!is.na(var.cluster) && any(duplicated(data[[var.cluster]]))){
            add.time <- TRUE
        }else if(!is.na(var.time) && sum(!duplicated(data[[var.time]]))>1){
            add.time <- TRUE
        }else{
            add.time <- FALSE
        }
        
        n.time <- data[[var.time]]
        if(is.na(var.strata)){
            if(structure %in% c("IND","UN") && add.time){
                structure <- do.call(type.structure, list(stats::as.formula(paste("~",var.time)), var.cluster = "XXclusterXX", var.time = "XXtimeXX"))
            }else{
                structure <- do.call(type.structure, list(formula = ~1, var.cluster = "XXclusterXX", var.time = "XXtimeXX"))
            }
        }else{
            if(structure %in% c("IND","UN") && add.time){
                structure <- do.call(type.structure, list(stats::as.formula(paste(var.strata,"~",var.time)), var.cluster = "XXclusterXX", var.time = "XXtimeXX"))
            }else{
                structure <- do.call(type.structure, list(stats::as.formula(paste(var.strata,"~1")), var.cluster = "XXclusterXX", var.time = "XXtimeXX"))
            }
        }
    }else{
        stop("Argument \'structure\' must either be a character or a structure object. \n")
    }

    ## *** mean structure
    if(trace>=2){cat("- mean structure")}
    if(!inherits(formula,"formula")){
        stop("Argument \'formula\' must be of class formula \n",
             "Something like: outcome ~ fixedEffect1 + fixedEffect2 \n")
    }
    name.mean <- all.vars(formula)
    if(any(name.mean %in% names(data) == FALSE)){
        invalid <- name.mean[name.mean %in% names(data) == FALSE]
        stop("Argument \'formula\' is inconsistent with argument \'data\'. \n",
             "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
             sep = "")
    }
    
    var.outcome <- lhs.vars(formula)
    if(length(var.outcome)!=1){
        stop("Argument \'formula\' must be contain exactly one outcome variable \n")
    }
    var.X <- rhs.vars(formula)
    if(any(grepl("Intercept",var.X)) ||any(grepl("(Intercept)",var.X))){
        stop("Argument \'formula\' should not contain a variable called \"Intercept\" or \"(Intercept)\". \n")
    }
    if(any(grepl(":",var.X,fixed=TRUE))){
        stop("Argument \'formula\' should not contain a variable whose name contain \":\". \n")
    }
    formula.terms <- stats::terms(formula)
    if(any(attr(formula.terms,"order")>2)){
        stop("Cannot handle interaction involving more than two variables. \n")
    }

    if(!is.na(var.strata)){
        if(optimizer == "gls"){
            var.X.withinStrata <- setdiff(var.X, var.strata)
            if(length(var.X.withinStrata)==0){
                formula.design <- stats::as.formula(paste0(var.strata,"~1")) ## no covariate within strata
            }else{
                terms.mean <- stats::terms(formula)
                newterm.labels <- gsub(paste0("\\:",var.strata,"$"),"",gsub(paste0("^",var.strata,"\\:"),"",setdiff(attr(terms.mean,"term.labels"),var.strata)))  ## remove interaction or main effect with the strata variable
                if(attr(terms.mean,"intercept")>0){
                    formula.design <- stats::update(formula, paste0(var.strata,"~",paste0(unique(newterm.labels),collapse=" + ")))
                }else{
                    formula.design <- stats::update(formula, paste0(var.strata,"~-1+",paste0(unique(newterm.labels),collapse=" + ")))
                }
            }
        }else{
            formula.design <- stats::as.formula(stats::delete.response(stats::terms(formula)))
        }
    }else{
        formula.design <- stats::as.formula(stats::delete.response(stats::terms(formula)))
    }
    out$formula <- list(mean = formula,
                        mean.design = formula.design,
                        var.design = structure$formula$var,
                        cor.design = structure$formula$cor)
    var.Z <- c(all.vars(out$formula$var.design),all.vars(out$formula$var.design))
    if(trace>=2){cat("\n")}

    ## *** data
    if(trace>=2){cat("- data")}
    data <- .prepareData(data, var.cluster = var.cluster, var.time = var.time, var.strata = var.strata)

    ## cluster
    U.cluster <- levels(data$XXclusterXX)
    n.cluster <- length(U.cluster)
    
    ## time
    U.time <- levels(data$XXtimeXX)
    n.time <- length(U.time)

    test.duplicated <- tapply(data$XXtimeXX, data$XXclusterXX, function(iT){any(duplicated(iT))})
    if(any(test.duplicated)){
        stop("Incorrect specification of argument \'repetition\'. \n",
             "The time variable (first variable before |) should contain unique values within clusters \n")
    }
    
    ## strata
    U.strata <- levels(data$XXstrataXX)
    n.strata <- length(U.strata)
        
    if(n.strata > 1){
        if(any(rowSums(table(data$XXclusterXX,data$XXstrata.indexXX)>0)!=1)){
            stop("When a variable is used to stratify the variance structure, all observations belonging to each cluster must belong to a single strata. \n")
        }

        if(optimizer=="gls"){
            tocheck <- setdiff(var.X,var.strata)
            if(var.strata %in% var.X == FALSE){
                stop("When a variable is used to stratify the variance structure, it should also be use to stratify the mean structure. \n",
                     "Consider adding all interactions with \"",var.strata,"\" in the argument \'formula\'. \n",
                     "Or using \"FS\" optimizer (see control argument). \n")
            }
            
            sapply(tocheck, function(iX){ ## iX <- "age"
                iCoef <- which(attr(formula.terms,"factors")[iX,]>=1)
                iInteraction <- attr(formula.terms,"factors")[var.strata,iCoef,drop=FALSE]
                if(all(iInteraction==0)){
                    stop("When a variable is used to stratify the variance structure, it should also be used to stratify the mean structure. \n",
                         "Consider adding an interaction between \"",iX,"\" and \"",var.strata,"\" in the argument \'formula\'. \n",
                         "Or using \"FS\" optimizer (see control argument). \n")
                }
            })
        }

        test.length <- tapply(data$XXtimeXX, data$XXstrata.indexXX, function(iT){list(unique(iT))})
            
        if(length(unique(sapply(test.length,length)))>1){
            stop("The time variable should contain the same number of unique values in each strata \n")
        }
        test.unique <- do.call(rbind,lapply(test.length, sort))
            
        if(any(apply(test.unique, MARGIN = 2, FUN = function(x){length(unique(x))})!=1)){
            stop("The time variable should contain the same unique values in each strata \n")
        }

        test.order <- do.call(rbind,test.length)
            
        if(any(apply(test.order, MARGIN = 2, FUN = function(x){length(unique(x))})!=1)){
            stop("The order of the unique values in the time variable should be the same in each strata \n")
        }
    }

    ## weights
    if(!is.null(weights)){
        if(optimizer=="gls"){
            stop("Cannot use argument \'weights\' with optimizer \"gls\".\n",
                 "Consider using \"FS\" optimizer instead (see control argument). \n")
        }
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

    ## scale.Omega
    if(!is.null(scale.Omega)){
        if(optimizer=="gls"){
            stop("Cannot use argument \'scale.Omega\' with optimizer \"gls\".\n",
                 "Consider using \"FS\" optimizer instead (see control argument). \n")
        }
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
        options$precompute.moments <- FALSE
    }else{
        var.scale.Omega <- NA
    }

    ## store
    out$strata <- list(n = n.strata, levels = U.strata, var = var.strata)
    out$time <- list(n = n.time, levels = U.time, var = var.time)
    out$cluster <- list(n = n.cluster, levels = U.cluster, var = var.cluster)
    out$weights <- list(var = c("logLik" = var.weights, "Omega" = var.scale.Omega))
    out$outcome <- list(var = var.outcome)
    out$data <- data

    ## *** missing values
    var.all <- unname(unique(stats::na.omit(c(var.strata,var.outcome,var.X,var.time,var.cluster,var.Z))))
    index.na <- which(rowSums(is.na(data[,var.all,drop=FALSE]))>0)
    if(length(index.na)>0){
        attr(index.na, "cluster") <- data$XXclusterXX[index.na]
        attr(index.na, "time") <- data$XXtimeXX[index.na]
        data <- data[-index.na,, drop=FALSE]
    }else{
        index.na <- NULL
    }

    out$index.na <- index.na

    if(trace>=2){cat("\n")}

    ## *** design matrices
    if(trace>=1){cat("- extract design matrices")}
    out$design <- .model.matrix.lmm(formula.mean = out$formula$mean.design,
                                    structure = structure,
                                    data = data, var.outcome = var.outcome, var.weights = out$weights$var,
                                    U.strata = U.strata,
                                    U.time = U.time,
                                    stratify.mean = optimizer=="gls",
                                    precompute.moments = options$precompute.moments)

    ## note use model.frame to handline splines in the formula
    out$xfactor <- c(list(mean = stats::.getXlevels(stats::terms(out$formula$mean.design),stats::model.frame(out$formula$mean.design,data))),
                     out$design$vcov$xfactor)
    out$design$vcov$xfactor <- NULL

    if(optimizer=="gls"){
        ## move from design matrix to dataset (useful when doing baseline adjustment)
        data.fit <- as.data.frame(out$design$mean)
        ## colnames(data.fit) <- gsub(" ","_",gsub("(Intercept)","Intercept",gsub(":","_",colnames(data.fit), fixed = TRUE), fixed = TRUE))
        colnames(data.fit) <- paste0("p",1:NCOL(data.fit))
    
        ## add outcome,strata,time,id to the dataset
        add.col <- unique(c(var.outcome, "XXclusterXX", "XXtimeXX", "XXtime.indexXX", "XXstrataXX", "XXstrata.indexXX", all.vars(out$formula$var.design),all.vars(out$formula$cor.design)))
        if(any(add.col %in% names(data.fit) == FALSE)){
            data.fit <- cbind(data.fit, data[,add.col[add.col %in% names(data.fit) == FALSE],drop=FALSE])
        }
        ## order by strata, time, and cluster (strata, cluster, and time does not provide satisfactory results, mixing k-parameters)
        col.pattern <- factor(out$design$vcov$X$pattern.cluster, out$design$vcov$X$Upattern$name)[as.character(data.fit[["XXclusterXX"]])]
        data.fit <- data.fit[order(data.fit[["XXstrata.indexXX"]],col.pattern,data.fit[["XXclusterXX"]],data.fit[["XXtime.indexXX"]]),,drop=FALSE]
        if(n.strata==1){
            txt.data <- "data.fit"
        }else{
            txt.data <- paste0("data.fit[data.fit$XXstrataXX==\"",U.strata,"\",,drop=FALSE]")
        }
        ## change cluster name so that they appear in order (important when there are missing values as this influences the reference level for the weights)
        data.fit[["XXclusterXX"]] <- as.factor(as.numeric(factor(data.fit[["XXclusterXX"]],levels=unique(data.fit[["XXclusterXX"]]))))
        ##  update formula
        txt.formula <- tapply(paste0("p",1:NCOL(out$design$mean)),out$design$param$strata.mu, function(iStrata){
            paste0(var.outcome, "~ 0 + ",  paste(iStrata, collapse = " + "))
        })
    }

    if(trace>=2){cat("\n")}


    ## ** Estimate model parameters
    if(trace>=1){cat("2. Estimate model parameters")}

    if(optimizer=="gls"){
        name.var <- stats::na.omit(setdiff(structure$name$var,structure$name$strata))
        if(length(name.var)>0){
            form.var <- stats::as.formula(paste0("~1|", paste(name.var[[1]], collapse = "*")))
            txt.var <- paste0("weights = nlme::varIdent(form = ", deparse(form.var), "), ")
        }else{
            txt.var <- NULL
        }
        if(max(out$design$cluster$nobs)==1 || type.structure %in% c("IND","ID")){
            txt.cor <- NULL
        }else if(type.structure == "CS"){
            form.cor <- ~1|XXclusterXX
            txt.cor <- paste0("correlation = nlme::corCompSymm(form = ",deparse(form.cor),"), ")
        }else if(type.structure == "UN"){
            form.cor <- ~XXtime.indexXX|XXclusterXX
            txt.cor <- paste0("correlation = nlme::corSymm(form = ",deparse(form.cor),"), ")
        }

        txt.gls <- paste0("nlme::gls(",txt.formula,", ",
                          txt.var, ## optional weights argument
                          txt.cor, ## optional correlation argument
                          "method = ",deparse(method.fit),", data = ",txt.data,", control = control)")

        out$gls <- stats::setNames(lapply(txt.gls, function(iTxt){eval(parse(text = iTxt))}),
                                   U.strata)
        out$gls.call <- lapply(out$gls, function(iM){
            paste0(gsub(",",",\n    ",gsub(" ","",paste(deparse(iM$call), collapse = ""))),"\n")
        })
        param.mu <- lapply(U.strata, function(iS){ coef(out$gls[[iS]]) })
        param.sigma <- lapply(U.strata, function(iS){ stats::sigma(out$gls[[iS]]) })
        param.value <- c(stats::setNames(unlist(param.mu),out$design$param$mu),
                         stats::setNames(unlist(param.sigma), out$design$param$sigma))
        if(length(out$design$param$k)>0){
            param.k <- lapply(U.strata, function(iS){ ## iS <- "1"
                iCoef <- coef(out$gls[[iS]]$modelStruct$varStruct, unconstrained = FALSE)
            })            
            param.value <- c(param.value,stats::setNames(unlist(param.k), out$design$param$k))
        }
        if(length(out$design$param$rho)>0){
            param.rho <- lapply(U.strata, function(iS){ coef(out$gls[[iS]]$modelStruct$corStruct, unconstrained = FALSE) })
            param.value <- c(param.value,stats::setNames(unlist(param.rho), out$design$param$rho))
        }        
        out$opt <- list(name = "gls")

    }else{
        outEstimate <- .estimate(design = out$design, time = out$time, method.fit = method.fit, type.information = type.information,
                                 transform.sigma = options$transform.sigma, transform.k = options$transform.k, transform.rho = options$transform.rho,
                                 precompute.moments = options$precompute.moments, 
                                 optimizer = optimizer, init = control$init, n.iter = control$n.iter, tol.score = control$tol.score, tol.param = control$tol.param, trace = control$trace)
        param.value <- outEstimate$estimate
        out$opt <- c(name = optimizer, outEstimate[c("cv","n.iter","score","previous.estimate")])
        
        if(out$opt$cv==FALSE){
            warning("Convergence issue: no stable solution has been found. \n")
        }
        
    }
    out$param <- list(value = param.value,
                      type = out$design$param$type,
                      strata = stats::setNames(out$strata$levels[out$design$param$strata],names(out$design$param$strata)))

    if(trace>=1){cat("\n")}

    ## ** Compute likelihood derivatives
    if(trace>=1){cat("3. Compute likelihood derivatives \n")}

    outMoments <- .moments.lmm(value = out$param$value, design = out$design, time = out$time, method.fit = method.fit, type.information = type.information,
                               transform.sigma = options$transform.sigma, transform.k = options$transform.k, transform.rho = options$transform.rho,
                               logLik = TRUE, score = TRUE, information = TRUE, vcov = TRUE, df = df, indiv = FALSE, effects = c("mean","variance","correlation"), robust = FALSE,
                               trace = trace>=2, precompute.moments = options$precompute.moments, method.numDeriv = options$method.numDeriv, transform.names = FALSE)
    
    out[names(outMoments)] <- outMoments

    
    

    if(trace>=1){cat("\n")}

    ## ** convert to lmm and export
    if(optimizer=="gls" && any(abs(out$score)>0.1)){
        warning("Large score value - incorrect model convergence or interface with nlme::gls. \n",
                "Consider switching to internal optimizer using control = list(optimizer = \"FS\") when calling LMM. \n")
    }
    class(out) <- "lmm"
    return(out)
}

## * .prepareData
## convert to data.frame
## generate factor time, cluster, strata and related indexes
.prepareData <- function(data, var.cluster, var.time, var.strata){

    ## ** convert to data.frame
    data <- as.data.frame(data)

    ## ** index
    data$XXindexXX <- 1:NROW(data)
    
    ## ** cluster
    if(is.na(var.cluster)){
        data$XXclusterXX <- as.factor(1:NROW(data))
    }else if(var.cluster %in% names(data)){
        if(is.factor(data[[var.cluster]])){
            data$XXclusterXX <- droplevels(data[[var.cluster]])
        }else{
            data$XXclusterXX <- factor(data[[var.cluster]], levels = sort(unique(data[[var.cluster]])))
        }
    }
    
    ## ** time
    if(is.na(var.time)){
        iTime <- tapply(data$XXclusterXX, data$XXclusterXX, function(iC){paste0("t",1:length(iC))})
        iIndex <- tapply(1:NROW(data), data$XXclusterXX, function(iC){iC})
        data[iIndex,"XXtimeXX"] <- as.factor(iTime)
        data$XXtime.indexXX <- as.numeric(data$XXtimeXX)
    }else if(var.time %in% names(data)){
        if(is.factor(data[[var.time]])){
            data$XXtimeXX <- droplevels(data[[var.time]])
        }else{
            data$XXtimeXX <- factor(data[[var.time]], levels = sort(unique(data[[var.time]])))
        }
        data$XXtime.indexXX <- as.numeric(data$XXtimeXX)
    }
    
    
    ## ** strata
    if(is.na(var.strata)){
        var.strata <- "XXstrata.indexXX"
        data$XXstrataXX <- factor(1)
        data$XXstrata.indexXX <- 1
    }else if(var.strata %in% names(data)){
        if(is.factor(data[[var.strata]])){
            data[["XXstrataXX"]] <- droplevels(data[[var.strata]])
        }else{
            data[["XXstrataXX"]] <- factor(data[[var.strata]], levels = sort(unique(data[[var.strata]])))
        }
        data[["XXstrata.indexXX"]] <- as.numeric(data[["XXstrataXX"]])
    }

    ## ** export
    return(data)

}

######################################################################
### lmm.R ends here
