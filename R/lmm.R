### lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: maj 12 2023 (09:43) 
##           By: Brice Ozenne
##     Update #: 2210
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

    ## ** 0. extract call and default from package
    out <- list(call = match.call(), data.original = data)
    options <- LMMstar.options()
    precompute.moments <- options$precompute.moments
    if(is.null(trace)){
        trace <- options$trace
    }
    
    ## ** 1. check and normalize user input
    if(trace>=1){cat("1. Check and normalize user input")}
    missing.repetition <- missing(repetition) ## argument repetition will be updated later which can be confusing
    missing.structure <- missing(structure) ## argument structure will be updated later which can be confusing   

    ## *** formula
    detail.formula <- formula2var(formula)
    formula <- detail.formula$formula$regressor ## remove possible random effects

    var.all_meanformula <- detail.formula$vars$all 
    var.X <- detail.formula$vars$regressor
    var.outcome <- detail.formula$vars$response

    if(detail.formula$special == "repetition"){
        stop("Random effects in argument \'formula\' should be wrapped into parentheses. \n",
             "Something like Y ~ X1 + (1|id). Otherwise consider using argument \'repetition'. \n",
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
    out$outcome <- list(var = var.outcome)

    if(any(grepl("Intercept",var.X)) ||any(grepl("(Intercept)",var.X))){
        stop("Argument \'formula\' should not contain a variable called \"Intercept\" or \"(Intercept)\". \n")
    }
    if(any(grepl(":",var.X,fixed=TRUE))){
        stop("Argument \'formula\' should not contain a variable whose name contain \":\". \n")
    }
       
    ## *** structure
    if(!missing.structure){
        if(detail.formula$special=="ranef"){
            stop("Argument \'formula\' and \'structure\' are inconsistent. \n",
                 "Either specify random effect (argument \'formula\') or a covariance structure (argument \'structure\') but not both. \n")
        }else if(inherits(structure,"character") || inherits(structure,"function")){
            structure <- do.call(structure, args = list(~1))
        }else if(inherits(structure,"structure")){
            ## nothing
        }else{
            stop("Argument \'structure\' must either be a character or a structure object. \n")
        }
        var.strata <- structure$name$strata
    }else if(detail.formula$special=="ranef"){
        if(length(detail.formula$vars$time)>0){
            stop("Incorrect argument \'formula\', \n",
                 "Current version can only handle random intercepts (i.e. no covariates in random effects). \n")
        }
        if(attr(detail.formula$special,"crossed") && attr(detail.formula$special,"nested")){
            stop("Incorrect argument \'formula\', \n",
                 "Current version cannot handle crossed and nested random effects. \n")
        }
        if(attr(detail.formula$special,"crossed")){            
            if(missing.repetition){
                repetition <- as.formula(paste0("~",paste(detail.formula$vars$cluster,collapse="+")))
            }
            ff.structure <- paste0("~",paste0(detail.formula$vars$ranef, collapse = "+"))
            structure <- CS(list(~1,as.formula(ff.structure)),
                            heterogeneous = -1,
                            ranef = detail.formula$special)
        }else{
            if(missing.repetition){
                repetition <- as.formula(paste0("~1|",detail.formula$vars$cluster))
            }
            structure.ranef <- list(type = data.frame(crossed = attr(detail.formula$special,"crossed"),
                                                      nested = attr(detail.formula$special,"crossed")),
                                    formula = detail.formula$formula$ranef,
                                    vars = detail.formula$vars$ranef,
                                    terms = detail.formula$terms$ranef,
                                    hierarchy = detail.formula$vars$hierarchy)
                                            
            if(attr(detail.formula$special,"nested")){
                ff.structure <- paste0("~",paste0(detail.formula$vars$ranef[-1], collapse = "+"))
                structure <- CS(as.formula(ff.structure),
                                heterogeneous = FALSE,
                                ranef = structure.ranef)
            }else{
                structure <- CS(~1,
                                ranef = structure.ranef)
            }
        }
        var.strata <- NA
    }else{
        var.strata <- NA
    }
    
    ## *** repetition
    if(missing(repetition)){
        update.strataStructure <- FALSE
        if(missing(structure)){
            var.cluster  <- NA
            var.time  <- NA
            var.strata  <- NA
            structure <- ID(~1)
        }else if(inherits(structure,"structure")){
            var.cluster <- structure$name$cluster
            var.time <- structure$name$time
            var.strata <- structure$name$strata
        }
    }else{
        if(inherits(try(repetition,silent=TRUE),"try-error")){
            stop("Could not evaluate argument \'repetition\'. Maybe the symbol \'~\' is missing in the formula. \n")
        }
        if(!inherits(repetition,"formula")){
            stop("Argument \'repetition\' must be of class formula, something like: ~ time or ~ time|cluster. \n")
        }
        ## extract variables from formula
        detail.repetition <- formula2var(repetition, name.argument = "repetition")
        if(detail.repetition$special=="ranef"){
            stop("Incorrect specification of argument \'repetition\'. \n",
                 "Should be something like: ~ time or ~ time|cluster. \n")
        }
        if(detail.repetition$special == "repetition"){
            var.cluster <- detail.repetition$var$cluster
            var.time <- detail.repetition$var$time
        }else if(detail.repetition$special=="none"){            
            var.cluster <- NA
            var.time <- detail.repetition$var$regressor
        }
        var.strata2 <- detail.repetition$var$response
        if(length(var.strata2)>0){ ## catch strata variable 
            if(any(!is.na(var.strata)) && !identical(var.strata,var.strata2)){
                stop("Inconsistency between the strata defined via the \'repetition\' and the \'structure\' argument. \n",
                     "\"",paste(var.strata2, collapse = "\" \""),"\" vs. \"",paste(var.strata, collapse = "\" \""),"\" \n")
            }else{
                update.strataStructure <- TRUE
                var.strata <- var.strata2
            }
        }else{
            update.strataStructure <- FALSE
        }
    }

    ## compatibility structure/repetition
    if(!missing.structure){
        
        if(all(is.na(var.cluster)) && structure$type %in% c("CS","UN")){
            stop("Incorrect specification of argument \'repetition\': missing cluster variable. \n",
                 "Should have exactly one variable after the grouping symbol (|), something like: ~ time|cluster or strata ~ time|cluster. \n")
        }
        if(all(is.na(var.time)) && structure$type %in% c("IND","UN") && all(is.na(structure$name$var[[1]]))){
            stop("Incorrect specification of argument \'repetition\': missing time variable. \n",
                 "Should have exactly one variable before the grouping symbol (|), something like: ~ time|cluster or strata ~ time|cluster. \n")
        }
    }

    ## sanity check: variability within cluster (for clusters with more than a single value)
    if(any(!is.na(var.cluster))){

        value.cluster <- as.character(interaction(data[var.cluster]))
        index.clusterMY <- which(value.cluster %in% unique(value.cluster[duplicated(value.cluster)]))

        test.rep <- tapply(data[[var.outcome]][index.clusterMY],index.clusterMY[index.clusterMY],function(iValue){
            ## identify clusters with constant value (including all NA)
            if(length(iValue)==1){
                return(FALSE)
            }else if(all(is.na(iValue))){
                return(TRUE)
            }else{
                return(abs(max(iValue,na.rm=TRUE)-min(iValue,na.rm=TRUE))<1e-12)
            }
        })
        if(length(test.rep)>0 & all(test.rep)){
            warning("Constant outcome value within cluster. \n")
        }
    }

    ## *** objective function
    if(is.null(method.fit)){
        if(length(var.X)==0 && attr(stats::terms(formula), "intercept") == 0){
            method.fit <- "ML"
        }else{
            method.fit <- options$method.fit
        }
    }else{
        method.fit <- match.arg(method.fit, choices = c("ML","REML"))
        if(length(var.X)==0 && method.fit == "REML"){
            message("Revert back to ML estimation as there is no mean structure. \n")
            method.fit <- "ML"
        }
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
    
    ## *** data
    data <- .prepareData(data,
                         var.cluster = var.cluster,
                         var.time = var.time,
                         var.strata = var.strata,
                         missing.repetition = missing.repetition,
                         droplevels = TRUE)

    ## cluster
    if(is.na(var.cluster)){
        var.clusterOLD <- var.cluster
        var.cluster <- "XXclusterXX"
        attr(var.cluster, "original") <- var.clusterOLD
    }else{
        attr(var.cluster, "original") <- var.cluster
    }
    U.cluster <- levels(data$XXclusterXX)
    n.cluster <- length(U.cluster)
    out$cluster <- list(n = n.cluster, levels = U.cluster, var = var.cluster)

    ## time
    if(all(is.na(var.time))){
        var.timeOLD <- var.time
        var.time <- "XXtimeXX"
        attr(var.time, "original") <- var.timeOLD
    }else{
        attr(var.time, "original") <- var.time
    }
    U.time <- levels(data$XXtimeXX)
    if(length(var.time)>0){
        ls.timeOriginal <- by(data[,var.time,drop=FALSE],data$XXtimeXX,function(iDF){iDF[1,,drop=FALSE]}, simplify = FALSE)
        attr(U.time,"original") <- do.call(rbind,ls.timeOriginal)
    }
    n.time <- length(U.time)
    out$time <- list(n = n.time, levels = U.time, var = var.time)

    ## strata
    if(is.na(var.strata)){
        var.strataOLD <- var.strata
        var.strata <- "XXstrataXX"
        attr(var.strata, "original") <- var.strataOLD
    }else{
        attr(var.strata, "original") <- var.strata
    }
    U.strata <- levels(data$XXstrataXX)
    n.strata <- length(U.strata)
    out$strata <- list(n = n.strata, levels = U.strata, var = var.strata)
     
    ## *** optimizer
    if(is.null(control$optimizer)){
        optimizer <- options$optimizer
    }else{
        optimx.method <- c("BFGS", "CG", "Nelder-Mead", "nlminb", "bobyqa")
        optimizer <- match.arg(control$optimizer, c("FS",optimx.method)) ## FS = fisher scoring
        control$optimizer <- NULL
    }
    if(is.null(control$trace)){
        trace.control <- trace-1
    }else{
        trace.control <- control$trace
    }

    ## *** weights
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

    ## *** scale.Omega
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
        precompute.moments <- FALSE
    }else{
        var.scale.Omega <- NA
    }

    out$weights <- list(var = c("logLik" = var.weights, "Omega" = var.scale.Omega))

    if(trace>=1){cat("\n")}

    ## ** 2. Design matrix and precomputation
    if(trace>=1){cat("2. Design matrix and precomputations \n")}

    ## *** mean structure
    if(trace>=2){cat("- formula for the mean structure")}

    formula.design <- stats::as.formula(stats::delete.response(stats::terms(formula)))
    out$formula <- list(mean = formula,
                        mean.design = formula.design)

    if(trace>=2){cat("\n")}

    ## *** residual variance-covariance structure
    if(trace>=2){cat("- residual variance-covariance structure")}

    if(missing(structure)){
        if(any(duplicated(data[["XXclusterXX"]]))){
            if(all(is.na(attr(var.time,"original")))){
                structure <- CS(~1)
            }else{
                structure <- UN(~1)
            }
        }else if(length(levels(data[["XXtimeXX"]]))>1){
            structure <- IND(~1)
        }else{
            structure <- ID(~1)
        }
    }

    if(!missing.repetition){
        if(any(!is.na(structure$name$time)) && any(sort(structure$name$time) != sort(var.time))){
            warning("Argument \'structure\' is inconsistent with argument \'repetition\'. \n",
                    "Not the same time variable: ",paste(structure$name$time, collapse = " ")," vs. ",paste(var.time, collapse =  " "),".\n")
        }
        if(!is.na(structure$name$cluster) && structure$name$cluster != var.cluster){
            warning("Argument \'structure\' is inconsistent with argument \'repetition\'. \n",
                    "Not the same cluster variable: ",structure$name$cluster," vs. ",var.cluster,".\n")
        }
    }
    type.structure <- structure$type
    call.structure <- as.list(structure$call)
    args.structure <- call.structure[-1]

    if(("add.time" %in% names(args.structure) == FALSE || identical(args.structure$add.time,TRUE)) && type.structure %in% c("IND","UN","EXP","TOEPLITZ") && n.time>1){
        if(type.structure == "TOEPLITZ" && ("heterogeneous" %in% names(args.structure) && is.null(args.structure$heterogeneous))){
            args.structure$add.time <- "XXtimeXX"
        }else{
            args.structure$add.time <- var.time
        }
    }
    if(is.na(structure$name$cluster)){
        args.structure$var.cluster <- "XXcluster.indexXX"
    }
    if(is.na(structure$name$cluster)){
        args.structure$var.time <- "XXtime.indexXX"
    }

    ## add strata to the call
    if(update.strataStructure){
        if(is.list(args.structure$formula)){
            args.structure$formula <- list(stats::update(stats::as.formula(args.structure$formula[[1]], stats::as.formula(paste0(var.strata2,"~.")))),
                                           stats::update(stats::as.formula(args.structure$formula[[2]], stats::as.formula(paste0(var.strata2,"~.")))))
        }else{
            args.structure$formula <- stats::update(stats::as.formula(args.structure$formula), stats::as.formula(paste0(var.strata2,"~.")))
        }
    }
    if(inherits(call.structure[[1]], "function")){
        structure <- do.call(call.structure[[1]], args = args.structure)
    }else{
        structure <- do.call(deparse(call.structure[[1]]), args = args.structure)
    }
    if(structure$type=="CUSTOM"){precompute.moments <- FALSE}
    
    out$formula$var.design <- structure$formula$var
    out$formula$cor.design <- structure$formula$cor
    var.Z <- c(all.vars(out$formula$var.design),all.vars(out$formula$cor.design))

    if(!identical(structure$name$var[[1]], NA) && any(structure$name$var[[1]] %in% names(data) == FALSE)){
        invalid <- structure$name$var[[1]][structure$name$var[[1]] %in% names(data) == FALSE]
        stop("Variance structure inconsistent with argument \'data\'. \n",
             "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
             sep = "")
    }
    if(!identical(structure$name$cor[[1]], NA) && any(structure$name$cor[[1]] %in% names(data) == FALSE)){
        invalid <- structure$name$cor[[1]][structure$name$cor[[1]] %in% names(data) == FALSE]
        stop("Correlation structure inconsistent with argument \'data\'. \n",
             "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
             sep = "")
    }

    ## update transformation
    if(structure$type=="CUSTOM" && (is.null(structure$d2FCT.sigma) || is.null(structure$d2FCT.rho)) && (df || method.fit=="REML" || type.information=="observed")){
        ## need second derivative but transformation based on dJacobian not implemented!
        options$transform.sigma <- "none"
        options$transform.k <- "none"
        options$transform.rho <- "none"
    }
    
    if(trace>=2){cat("\n")}

    ## *** missing values
    if(trace>=2){cat("- remove missing values")}
    var.all <- unname(unique(stats::na.omit(c(var.strata,var.outcome,var.X,var.time,var.cluster,var.Z))))
    index.na <- which(rowSums(is.na(data[,var.all,drop=FALSE]))>0)
    data.save <- data

    if(length(index.na) == NROW(data)){
        var.na <- var.all[colSums(!is.na(data[,var.all,drop=FALSE]))==0]
        if(length(var.na)==0){
            stop("All observations have at least one missing data. \n")
        }else if(length(var.na)==1){
            stop("Variable \"",var.na,"\" contains only missing data. \n")
        }else{
            stop("Variables \"",paste(var.na, collapse="\" \""),"\" contain only missing data. \n")
        }
    }

    test.naOutcome <- is.na(data[[var.outcome]])
    test.naOther <- rowSums(is.na(data[setdiff(var.all,var.outcome)]))>0
    if( any(test.naOther > test.naOutcome) ){
        index.row <- which(test.naOther > test.naOutcome)
        test.naOther2 <- colSums(is.na(data[index.row,setdiff(var.all,var.outcome),drop=FALSE]))
        name.naOther2 <- names(test.naOther2)[test.naOther2>0]

        if(length(name.naOther2)==1){
            warning("Can only handle missing values in the outcome variable. \n",
                    "Observation(s) with missing values in \"",name.naOther2,"\" will be removed. \n")
        }else{
            warning("Can only handle missing values in the outcome variable. \n",
                    "Observation(s) with missing values in \"",paste(name.naOther2, collapse="\" \""),"\" will be removed. \n")
        }
    }

    if(length(index.na)>0){        
        attr(index.na, "cluster") <- data[index.na,var.cluster]
        attr(index.na, "cluster.index") <- data[index.na,"XXcluster.indexXX"]
        attr(index.na, "cluster.index")[attr(index.na, "cluster.index") %in% unique(data[-index.na,"XXcluster.indexXX"]) == FALSE] <- NA
        
        attr(index.na, "time") <- data[index.na,var.time]
        attr(index.na, "time.index") <- data[index.na,"XXtime.indexXX"]
        attr(index.na, "time.index")[attr(index.na, "time.index") %in% unique(data[-index.na,"XXtime.indexXX"]) == FALSE] <- NA

        data$XXindexXX <- NULL
        data$XXclusterXX <- NULL
        data$XXcluster.indexXX <- NULL
        data$XXtimeXX <- NULL
        data$XXtime.indexXX <- NULL
        data$XXstrataXX <- NULL
        data$XXstrata.indexXX <- NULL

        data <- .prepareData(data[-index.na,, drop=FALSE],
                             var.cluster = attr(var.cluster,"original"),
                             var.time = attr(var.time,"original"),
                             var.strata = attr(var.strata,"original"),
                             missing.repetition = missing.repetition,
                             droplevels = TRUE)
    }else{
        index.na <- NULL
    }

    out$data <- data.save
    out$index.na <- index.na

    if(trace>=2){cat("\n")}

    ## *** design matrices
    if(trace>=2){cat("- design matrices")}
    out$design <- .model.matrix.lmm(formula.mean = out$formula$mean.design,
                                    structure = structure,
                                    data = data, var.outcome = var.outcome, var.weights = out$weights$var,
                                    precompute.moments = precompute.moments,
                                    drop.X = options$drop.X)

    if(!is.na(attr(var.cluster,"original"))){
        if(is.factor(data[[attr(var.cluster,"original")]])){
            out$design$cluster$levels.original <- levels(data.save[[attr(var.cluster,"original")]])
        }else{
            out$design$cluster$levels.original <- sort(unique(data.save[[attr(var.cluster,"original")]]))
        }
    }
    if(any(!is.na(attr(var.time,"original")))){
        if(length(attr(var.time,"original"))==1){
            if(is.factor(data[[attr(var.time,"original")]])){
                out$design$time$levels.original <- levels(data.save[[attr(var.time,"original")]])
            }else{
                out$design$time$levels.original <- sort(unique(data.save[[attr(var.time,"original")]]))
            }
        }else{
            out$design$time$levels.original <- sort(unique(as.character(interaction(data.save[,attr(var.time,"original")], drop = TRUE))))
        }
    }
    ## note use model.frame to handline splines in the formula
    out$xfactor <- c(list(mean = stats::.getXlevels(stats::terms(out$formula$mean.design),stats::model.frame(out$formula$mean.design,data))),
                     out$design$vcov$xfactor)
    out$design$vcov$xfactor <- NULL
    
    if(trace>=2){cat("\n")}


    ## ** 3. Estimate model parameters
    if(trace>=1){cat("3. Estimate model parameters")}

    valid.control <- c("init","n.iter","tol.score","tol.param","trace")
    if(any(names(control) %in% valid.control  == FALSE)){
        stop("Incorrect elements in argument \'control\': \"",paste(names(control)[names(control) %in% valid.control  == FALSE], collapse = "\" \""),"\". \n",
             "Valid elements: \"",paste(valid.control, collapse = "\" \""),"\".\n")
    }
    if(identical(control$init,"lmer")){
        ## check feasibility
        requireNamespace("lme4")
        if(inherits(out$design$vcov,"RE")){
            stop("Initializer \"lmer\" only available for random effect structures.")
        }
        if(n.strata>1){
            stop("Initializer \"lmer\" cannot handle multiple strata.")
        }
        if(attr(detail.formula$special,"crossed") & attr(detail.formula$special,"nested")){
            stop("Initializer \"lmer\" cannot handle both crossed and random effects. \n")
        }

        ## fit lmer model
        e.lmer <- lme4::lmer(detail.formula$formula$all, data = data, REML = method.fit=="REML")

        ## extract lmer estimates
        lmer.beta <- lme4::fixef(e.lmer)
        lmer.sigma <- stats::setNames(stats::sigma(e.lmer)^2,"sigma")
        lmer.tau <- stats::setNames((sapply(lme4::VarCorr(e.lmer),attr,"stddev"))^2, sapply(strsplit(names(lme4::VarCorr(e.lmer)),split=":"),"[",1))

        ## convert to LMMstar estimates
        init.sigma <- sqrt(lmer.sigma+sum(lmer.tau))
        browser()
        if(attr(detail.formula$special,"crossed")){
            init.tau <- cumsum(rev(lmer.tau))/init.sigma^2
        }else if(attr(detail.formula$special,"crossed")){
            init.tau <- lmer.tau/lmer.sigma^2
        }else{
            init.tau <- lmer.tau/lmer.sigma^2
        }
        names(init.tau) <- table.nesting$param[match(gsub("(Intercept)",out$cluster$var,table.nesting$variable, fixed = TRUE),names(init.tau))]
        if(!identical(sort(names(c(lmer.beta,init.sigma,init.tau))), sort(out$design$param$name))){
            stop("Could not identify all coefficients from the lmer model. \n")
        }
            
        control$init <- c(lmer.beta,init.sigma,init.tau)[out$design$param$name]
    }
    if(trace>0){
        if(trace.control>0){cat("\n")}
        if(trace.control>1){cat("\n")}
    }
    outEstimate <- .estimate(design = out$design, time = out$time, method.fit = method.fit, type.information = type.information,
                             transform.sigma = options$transform.sigma, transform.k = options$transform.k, transform.rho = options$transform.rho,
                             precompute.moments = precompute.moments, 
                             optimizer = optimizer, init = control$init, n.iter = control$n.iter, tol.score = control$tol.score, tol.param = control$tol.param, trace = trace.control)
    param.value <- outEstimate$estimate
    out$opt <- c(name = optimizer, outEstimate[c("cv","n.iter","score","previous.estimate","previous.logLik","control")])
    if((trace==0 && trace.control>0)){
        cat("\n")
    }
    if(out$opt$cv<=0){
        warning("Convergence issue: no stable solution has been found. \n")
    }
        
    out$param <- param.value

    if(trace>=1){cat("\n")}

    ## ** 4. Compute likelihood derivatives
    if(trace>=1){cat("4. Compute likelihood derivatives \n")}
    outMoments <- .moments.lmm(value = out$param, design = out$design, time = out$time, method.fit = method.fit, type.information = type.information,
                               transform.sigma = options$transform.sigma, transform.k = options$transform.k, transform.rho = options$transform.rho,
                               logLik = TRUE, score = TRUE, information = TRUE, vcov = TRUE, df = df, indiv = FALSE, effects = c("mean","variance","correlation"), robust = FALSE,
                               trace = trace>=2, precompute.moments = precompute.moments, method.numDeriv = options$method.numDeriv, transform.names = FALSE)
    
    out[names(outMoments)] <- outMoments
    

    if(trace>=1){cat("\n")}

    ## ** 5. convert to lmm and export
    class(out) <- "lmm"
    return(out)
}

## * .prepareData
## convert to data.frame
## generate factor time, cluster, strata and related indexes
.prepareData <- function(data, var.cluster, var.time, var.strata, missing.repetition, droplevels){

    ## ** convert to data.frame
    if(!inherits(data,"data.frame")){
        stop("Argument \'data\' should be a data.frame or inherit from it. \n")
    }
    if(is.null(var.cluster)){var.cluster <- NA}
    if(is.null(var.time)){var.time <- NA}
    if(is.null(var.strata)){var.strata <- NA}
    data <- as.data.frame(data)
    
    ## ** tests
    if(is.null(missing.repetition)){
        name.arg <- NULL
    }else if(missing.repetition){
        name.arg <- "structure"
    }else{
        name.arg <- "repetition"
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
    if("XXcluster.indexXX" %in% names(data)){
        stop("Argument \'data\' should not contain a column named \"XXcluster.indexXX\" as this name is used internally by the lmm function. \n")
    }
    if("XXstrataXX" %in% names(data)){
        stop("Argument \'data\' should not contain a column named \"XXstrataXX\" as this name is used internally by the lmm function. \n")
    }
    if("XXstrata.indexXX" %in% names(data)){
        stop("Argument \'data\' should not contain a column named \"XXstrata.indexXX\" as this name is used internally by the lmm function. \n")
    }
    if(any(!is.na(var.cluster)) && any(var.cluster %in% names(data) == FALSE)){

        if(!is.null(name.arg)){
            stop("Argument \'",name.arg,"\' is inconsistent with argument \'data\'. \n",
                 "Could not find column \"",paste(var.cluster, collapse = "\" \""),"\" indicating the cluster in argument \'data\'. \n", sep="")
        }else{
            stop("Could not find column \"",paste(var.cluster, collapse = "\" \""),"\" indicating the cluster in argument \'data\'. \n", sep="")
        }
    }
    if(any(!is.na(var.time)) && any(var.time %in% names(data) == FALSE)){
        if(!is.null(name.arg)){
            stop("Argument \'",name.arg,"\' is inconsistent with argument \'data\'. \n",
                 "Could not find column \"",paste(var.time, collapse = "\" \""),"\" indicating the time in argument \'data\'. \n", sep="")
        }else{
            stop("Could not find column \"",paste(var.time, collapse = "\" \""),"\" indicating the time in argument \'data\'. \n", sep="")
        }
    }
    if(any(!is.na(var.strata)) && any(var.strata %in% names(data) == FALSE)){
        if(!is.null(name.arg)){
            stop("Argument \'",name.arg,"\' is inconsistent with argument \'data\'. \n",
                 "Could not find column \"",paste(var.strata, collapse = "\" \""),"\" indicating the strata in argument \'data\'. \n",sep="")
        }else{
            stop("Could not find column \"",paste(var.strata, collapse = "\" \""),"\" indicating the strata in argument \'data\'. \n",sep="")
        }
    }
    if(length(var.cluster)>1){
        if(!is.null(missing.repetition) && missing.repetition){
            stop("Incorrect specification of argument \'",name.arg,"\': too many cluster variables. \n")
        }else{
            stop("Incorrect specification of argument \'",name.arg,"\': too many cluster variables. \n",
                 "There should be exactly one variable after the grouping symbol (|), something like: ~ time|cluster or strata ~ time|cluster. \n", sep = "")
        }
    }
    
    ## ** index
    data$XXindexXX <- 1:NROW(data)
       
    ## ** cluster
    if(is.na(var.cluster)){
        if(identical(attr(var.cluster,"crossed"),TRUE)){
            if(any(duplicated(data[,var.time]))){
                stop("Incorrect specification of argument \'repetition\': missing cluster variable. \n",
                     "Cannot guess the cluster variable with non-unique levels for the cross random effects. \n")
            }
            data$XXclusterXX <- as.factor("1")
        }else{
            data$XXclusterXX <- as.factor(sprintf(paste0("%0",ceiling(log10(NROW(data)))+0.1,"d"), 1:NROW(data)))
        }
    }else if("XXclusterXX" %in% names(data) == FALSE){
        if(is.factor(data[[var.cluster]])){
            if(droplevels){
                data$XXclusterXX <- droplevels(data[[var.cluster]])
            }else{
                data$XXclusterXX <- data[[var.cluster]]
            }
        }else if(is.numeric(data[[var.cluster]]) && all(data[[var.cluster]] %% 1 == 0)){
            data$XXclusterXX <- as.factor(sprintf(paste0("%0",ceiling(log10(max(abs(data[[var.cluster]])))+0.1),"d"), data[[var.cluster]]))
        }else{
            data$XXclusterXX <- factor(data[[var.cluster]], levels = sort(unique(data[[var.cluster]])))
        }
    }
    data$XXcluster.indexXX <- as.numeric(droplevels(data$XXclusterXX))

    ## ** time
    if(length(var.time)>1){
        ## create new variable summarizing all variables
        data[["XXtimeXX"]] <- interaction(lapply(var.time, function(iX){as.factor(data[[iX]])}), drop = TRUE)
        
        ## ## try to handle the case where time in cluster 1 is (C1,S1) (C1,S2) (C1,S3) while it is (C2,S1) (C2,S2) (C2,S3)
        ## ## i.e. unify time across clusters
        ## interaction.tempo <-  <- interaction(lapply(var.time, function(iX){as.factor(data[[iX]])}), drop = TRUE)
        ## data[["XXtime.indexXX"]] <- NA
        ## for(iCluster in unique(data[[var.cluster]])){ ## iCluster <- 1
        ##     data[["XXtime.indexXX"]][data[[var.cluster]]==iCluster] <- as.numeric(droplevels(interaction.tempo[data[[var.cluster]]==iCluster]))
        ## }
        ## data[["XXtimeXX"]] <- factor(data[["XXtime.indexXX"]], labels = levels(interaction.tempo)[1:max(data[["XXtime.indexXX"]])])

    }else if(is.na(var.time) || (identical(var.time,"XXtimeXX") && "XXtimeXX" %in% names(data) == FALSE)){
        iTime <- tapply(data$XXclusterXX, data$XXclusterXX, function(iC){1:length(iC)})
        iIndex <- split(1:NROW(data), data$XXclusterXX)
        if(is.list(iIndex)){
            iIndex <- unlist(iIndex)
            iTime <- unlist(iTime)
        }
        data[iIndex,"XXtimeXX"] <- as.factor(sprintf(paste0("%0",ceiling(log10(max(iTime))+0.1),"d"), iTime))
        data$XXtime.indexXX <- as.numeric(droplevels(data$XXtimeXX))
    }else if(is.null(missing.repetition)){
        ## from model.matrix: prediction for new data where the time variable has already be nicely reformated to factor
        data$XXtimeXX <- data[[var.time]]
    }else{
        if(is.factor(data[[var.time]])){
            if(droplevels){
                data$XXtimeXX <- droplevels(data[[var.time]])
            }else{
                data$XXtimeXX <- data[[var.time]]
            }
        }else if(is.numeric(data[[var.time]]) && all(data[[var.time]] %% 1 == 0)){
            data$XXtimeXX <- as.factor(sprintf(paste0("%0",ceiling(log10(max(abs(data[[var.time]])))+0.1),"d"), data[[var.time]]))
        }else{
            data$XXtimeXX <- factor(data[[var.time]], levels = sort(unique(data[[var.time]])))
        }
    }
    data$XXtime.indexXX <- as.numeric(data$XXtimeXX)

    ## ** strata
    if(length(var.strata)>1){
        ## create new variable summarizing all variables
        data[["XXstrataXX"]] <- interaction(lapply(var.strata, function(iX){as.factor(data[[iX]])}), drop = TRUE)
    }else if(is.na(var.strata) || (identical(var.strata,"XXindexXX") && "XXindexXX" %in% names(data) == FALSE)){
        var.strata <- "XXstrata.indexXX"
        data$XXstrataXX <- factor(1)
    }else{
        if(is.factor(data[[var.strata]]) && !is.null(missing.repetition)){
            if(droplevels){
                data$XXstrataXX <- droplevels(data[[var.strata]])
            }else{
                data$XXstrataXX <- data[[var.strata]]
            }
        }else if(is.numeric(data[[var.strata]]) && all(data[[var.strata]] %% 1 == 0)){
            data$XXstrataXX <- as.factor(sprintf(paste0("%0",ceiling(log10(max(abs(data[[var.strata]])))+0.1),"d"), data[[var.strata]]))
        }else{
            data$XXstrataXX <- factor(data[[var.strata]], levels = sort(unique(data[[var.strata]])))
        }
    }
    data$XXstrata.indexXX <- as.numeric(data$XXstrataXX)

    ## ** final checks
    n.cluster <- max(data$XXcluster.indexXX)
    n.time <- max(data$XXtime.indexXX)
    n.strata <- max(data$XXstrata.indexXX)
    if(n.cluster >= n.time){
        test.duplicated <- tapply(data$XXcluster.indexXX, data$XXtime.indexXX, function(iT){any(duplicated(iT))})
    }else{
        test.duplicated <- tapply(data$XXtime.indexXX, data$XXcluster.indexXX, function(iT){any(duplicated(iT))})
    }
    if(any(test.duplicated)){
        stop("Argument \'",name.arg,"\' is inconsistent with argument \'data\'. \n",
             "The time variable should contain unique values within clusters \n", sep="")
    }

    if(n.strata>1 & n.time >1){
        test.sameStrata <- tapply(data$XXstrata.indexXX, data$XXcluster.indexXX, function(iT){any(iT[1]!=iT[-1])})
        if(any(test.sameStrata)){
            stop("Argument \'",name.arg,"\' is inconsistent with argument \'data\'. \n",
                 "When a variable is used to stratify the variance structure, all observations within a cluster must belong to same strata. \n")
        
        }
    }
    ## library(microbenchmark)
    ## microbenchmark(a = tapply(data$XXcluster.indexXX, data$XXtime.indexXX, function(iT){any(duplicated(iT))}),
    ##                b = tapply(data$XXtime.indexXX, data$XXcluster.indexXX, function(iT){any(duplicated(iT))}))
    ## microbenchmark(a = any(tapply(data$XXstrata.indexXX, data$XXcluster.indexXX, function(iT){any(iT[1]!=iT[-1])})),
    ##                b = any(rowSums(table(data$XXcluster.indexXX,data$XXstrata.indexXX)>0)!=1))
    
    ## ** export
    return(data)
}

######################################################################
### lmm.R ends here
