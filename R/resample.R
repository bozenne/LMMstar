### resample.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 31 2022 (10:09) 
## Version: 
## Last-Updated: jul 18 2024 (16:25) 
##           By: Brice Ozenne
##     Update #: 784
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * resample (documentation)
##' @title Inference via Resampling for Linear Mixed Model
##' @description Non-parametric bootstrap or permutation test for Linear Mixed Models.
##' @name resample
##'
##' @param object a \code{lmm} object.
##' @param type [character] should permutation test (\code{"perm-var"} or \code{"perm-res"}) or non-parametric bootstrap (\code{"boot"}) be used?
##' @param effects [character vector] the variable(s) to be permuted or the effect(s) to be tested via non-parametric bootstrap.
##' Can also be a function of the model parameters when performing non-parametric bootstrap.
##' @param n.sample [integer] the number of samples used.
##' @param studentized [logical] should a studentized boostrap or permutation test be used?
##' @param trace [logical] should the execution of the resampling be traced?
##' @param seed [integer, >0] Random number generator (RNG) state used when starting resampling.
##' @param cpus [integer, >0] number of child-processes for parallel evaluation.
##' If \code{NULL} no state is set.
##' @param export.cpus [character vector] name of the variables to export to each cluster.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details All approach are carried at the cluster level: \itemize{
##' \item Bootstrap: sampling with replacement clusters. If a cluster is picked twice then different cluster id is used for each pick.
##' \item Permutation: permuting covariate values between clusters (this only lead to the null hypothesis when the covariate values are constant within clusters)
##' or assigning new outcome values by, under the null, permuting the normalized residuals, rescaling to residuals, and adding back the permuted fixed effect (any mean effect under H1 would be 0 because of the permutation if the variance-covariance structure is correct). The later procedure originates from Oliver et al (2012).
##' }
##'
##' The studentized bootstrap keep the original estimate and standard error but uses the samples to evaluates the quantiles of the distribution used to form the confidence intervals.
##' The studentized permutation test approximate the distribution of the test statistic under the null (instead of the distribution of the estimate under the null.).
##'
##' P-values for the bootstrap are computed by centering the bootstrap distribution of the estimate or test statistic around 0 and evaluating the frequency at which it takes values more extreme than the observed estimate or test statistics.
##' 
##' @references Oliver E. Lee and Thomas M. Braun (2012), \bold{Permutation Tests for Random Effects in Linear Mixed Models}. \emph{Biometrics}, Journal 68(2).
##'
##' @keywords htest
##' 
##' @examples
##' \dontrun{
##'
##' #### univariate linear regression ####
##' data(gastricbypassW, package = "LMMstar")
##' ## rescale to ease optimization
##' gastricbypassW$weight1 <- scale(gastricbypassW$weight1)
##' gastricbypassW$weight2 <- scale(gastricbypassW$weight2)
##' gastricbypassW$glucagonAUC1 <- scale(gastricbypassW$glucagonAUC1)
##'
##' e.lm <- lmm(weight2~weight1+glucagonAUC1, data = gastricbypassW)
##' model.tables(e.lm)
##'
##' ## non-parametric bootstrap
##' resample(e.lm, type = "boot", effects = c("weight1","glucagonAUC1"), seed = 10)
##' ## permutation test
##' resample(e.lm, type = "perm-var", effects = "weight1", seed = 10) 
##' resample(e.lm, type = "perm-var", effects = "glucagonAUC1", seed = 10)
##' ## using multiple cores
##' resample(e.lm, type = "boot", effects = c("weight1","glucagonAUC1"), cpus = 4)
##' 
##' #### random intercept model ####
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$weight <- scale(gastricbypassL$weight)
##' gastricbypassL$glucagonAUC <- scale(gastricbypassL$glucagonAUC)
##' gastricbypassL$gender <- as.numeric(gastricbypassL$id) %% 2
##' gastricbypassLR <- na.omit(gastricbypassL)
##' 
##' eCS.lmm <- lmm(weight~glucagonAUC+gender, data = gastricbypassLR,
##'                repetition = ~visit|id, structure = "CS")
##' model.tables(eCS.lmm)
##'
##' ## non-parametric bootstrap
##' resample(eCS.lmm, type = "boot", effects = c("glucagonAUC","gender"), seed = 10, trace = FALSE)
##' ## permutation test
##' resample(eCS.lmm, type = "perm-var", effects = "gender", seed = 10)
##' resample(eCS.lmm, type = "perm-res", effects = "glucagonAUC", seed = 10) 
##' }
##' 
#' @export
`resample` <-
  function(object, type, ...) UseMethod("resample")


## * resample.lmm (code)
##' @export
##' @rdname resample
resample.lmm <- function(object, type, effects, n.sample = 1e3, studentized = TRUE, 
                         trace = TRUE, seed = NULL, cpus = 1, export.cpus = NULL,
                         ...){

    options <- LMMstar.options()

    ## ** check user input
    ## dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## type
    if(type=="bootstrap"){
        type <- "boot"
    }
    type <- match.arg(type, c("perm-var","perm-res","boot"))
    ## boot: non-parametric bootstrap
    ## perm-var: exchange covariate value between clusters (need to be constant within cluster)
    ## perm-res: outcome becomes permuted normalized residuals under the null to which the (permuted) fixed effect under the null are added.

    name.meanvar <- attr(object$design$mean,"variable")
    value.meancoef <- coef(object, effects = "mean")
    sd.meancoef <- sqrt(diag(vcov(object, effects = "mean")))

    name.meancoef <- names(value.meancoef)
    name.coef <- names(coef(object, effects = "all", transform.sigma = FALSE, transform.k = FALSE, transform.rho = FALSE))
    if(type == "perm-var"){
        if(is.function(effects)){
            stop("Argument \'effects\' cannot be a function when using permutation. \n",
                 "Consider setting the argument \'type\' to \"boot\". \n")
        }
        if(any(is.character(effects)==FALSE)){
            stop("Argument \'effects\' should contain character strings refering to a variable for the mean structure. \n")
        }
        if(length(effects)==0){
            stop("Argument \'effects\' should have length at least 1. \n")
        }
        if(any(effects %in% name.meanvar == FALSE)){
            stop("Argument \'effects\' should be one of \"",paste(name.meanvar, collapse = "\", \""),"\". \n")
        }

        M.factor <- attr(stats::terms(object$formula$mean.design),"factor")
        colnames(M.factor) <- setdiff(name.meancoef,"(Intercept)")

        index.keepcoef <- which(colSums(M.factor[effects,,drop=FALSE]!=0)>0)
        name.keepcoef <- names(index.keepcoef)

        if(any(name.keepcoef %in% c("seed","sample","convergence"))){
            stop("The exported coefficient(s) should not be named \"",paste(c("seed","sample","convergence")[c("seed","sample","convergence") %in% name.keepcoef], collapse = "\" \""),"\". \n",
                 "This name is reserved for internally use and storage of the results. \n")
        }
        if(any(name.keepcoef %in% paste0("se.",name.keepcoef))){
            stop("The exported coefficient(s) should not be named with the prefix \".se\" at it will be used internally to store the standard errors. \n")
        }

        ls.variableName <- list(variance = stats::variable.names(object, effects = "variance"),
                                correlation = stats::variable.names(object, effects = "correlation"))
        if(!is.null(ls.variableName$variance) & !is.null(ls.variableName$correlation)){
            effects.vcov <- any(effects %in% union(ls.variableName$variance, ls.variableName$correlation))
        }else{
            effects.vcov <- FALSE
        }
        effects.fct <- FALSE

        ## test whether the covariate is constant within cluster
        test.Wvar <- unlist(by(object$data[effects], object$data[[object$cluster$var]], function(iData){sum(!duplicated(iData))>1}, simplify = FALSE))
        if(any(test.Wvar)){
            stop("The covariate(s) indicated by the argument \"effects\" vary within cluster. \n",
                 "This is not support when using permutations. Consider using type=\"boot\" instead of type=\"perm\". \n")
        }
        Uvar <- by(object$data[effects], object$data[[object$cluster$var]], function(iData){iData[1,,drop=FALSE]}, simplify = FALSE)

    }else if(type == "perm-res"){
        if(is.function(effects)){
            stop("Argument \'effects\' cannot be a function when using permutation. \n",
                 "Consider setting the argument \'type\' to \"boot\". \n")
        }
        if(any(is.character(effects)==FALSE)){
            stop("Argument \'effects\' should contain character strings refering to a parameter of the mean structure. \n")
        }
        if(any(effects %in% name.meancoef == FALSE)){
            stop("Argument \'effects\' should be one of \"",paste(name.meancoef, collapse = "\", \""),"\". \n")
        }
        effects.vcov <- FALSE
        effects.fct <- FALSE
        name.keepcoef <- effects
    }else if(type == "boot"){

        if(is.function(effects)){
            effects.formals <- names(formals(effects))            
            if(studentized){
                if(length(effects.formals) != 2){
                    stop("When a function, argument \'effects\' should have two arguments (studentized bootstrap) \n",
                         "The first argument refers to the estimates and the second to thevariance-covariance matrix of the estimates \n")
                }
                effects.estimate <- effects(coef(object, effects = "all", transform.sigma = FALSE, transform.k = FALSE, transform.rho = FALSE),
                                            vcov(object, effects = "all", transform.sigma = FALSE, transform.k = FALSE, transform.rho = FALSE))
                if(!is.matrix(effects.estimate) || !is.numeric(effects.estimate) || NROW(effects.estimate)!=2){
                    stop("The output of the function defined in the argument \'effects\' must be a numeric matrix with two rows (studentized bootstrap). \n",
                         "The first row should contain estimates and the second row corresponding standard errors.")
                }
                if(is.null(colnames(effects.estimate))){
                    stop("The columns of the matrix output by the function defined in the argument \'effects\' must be named. \n")
                }
                name.keepcoef <- colnames(effects.estimate)
                effects.fct <- TRUE
            }else if(!studentized){
                if(length(effects.formals) != 1){
                    stop("When a function, argument \'effects\' should have one or two arguments (non-studentized bootstrap) \n",
                         "The argument refers to the estimates. \n")
                }
                effects.estimate <- effects(coef(object, effects = "all", transform.sigma = FALSE, transform.k = FALSE, transform.rho = FALSE))
                if(!is.vector(effects.estimate) || !is.numeric(effects.estimate)){
                    stop("The output of the function defined in the argument \'effects\' must be a numeric vector (non-studentized bootstrap). \n")
                }
                if(is.null(names(effects.estimate))){
                    stop("The output of the function defined in the argument \'effects\' must be named. \n")
                }
                name.keepcoef <- names(effects.estimate)
                effects.fct <- TRUE
            }
            
        }else{
            if(any(is.character(effects)==FALSE)){
                stop("Argument \'effects\' should be a function or a character strings refering to model parameters. \n")
            }
            if(any(effects %in% name.coef == FALSE)){
                if(grepl("=",effects)){
                    stop("Argument \'effects\': \"",paste(effects[effects %in% name.coef == FALSE], collapse = "\", \""),"\" refer to model parameter(s) not to an equation. \n",
                         "Consider removing: =",strsplit(effects,split = "=")[[1]][2],"\n")
                }else{
                    stop("Incorrect argument \'effects\': \"",paste(effects[effects %in% name.coef == FALSE], collapse = "\", \""),"\" does not match any model parameter. \n")
                }
            }
            name.keepcoef <- effects
            effects.fct <- FALSE
        }
        effects.vcov <- TRUE
        
    }

    ## cpus
    max.cpus <- parallel::detectCores()
    if(length(cpus)!=1){
        stop("Argument \'cpus\' should have length 1.\n ")
    }else if(identical(cpus,"all")){
        cpus <- max.cpus
    }else if(!is.numeric(cpus) || cpus <=0 || cpus %% 1 != 0){
        stop("Argument \'cpus\' should be an integer between 1 and ",max.cpus," or \'all\'.\n ")
    }else if(cpus>1 && cpus>parallel::detectCores()){
        stop("Argument \'cpus\' exceeds the number of available CPUs.\n ",
             "It should be an integer between 1 and ",max.cpus," or \'all\'.\n ")
    }

    ## n.sample
    if(length(n.sample)!=1){
        stop("Argument \'n.sample\' should have lenght 1. \n")
    }
    if(!is.numeric(n.sample) || n.sample<0 || n.sample %% 1 >0 ){
        stop("Argument \'n.sample\' should be a positive integer. \n")
    }

    ## null
    if(type == "boot"){
        all.null <- model.tables(object, effects = "all", columns = "null", transform.k = "none", transform.rho = "none")
        if(is.function(effects)){
            null <- try(effects(all.null$null), silent = TRUE)
            if(inherits(null, "try-error")){
                null <- stats::setNames(rep(NA, length(name.keepcoef)),name.keepcoef)
            }
        }else{
            null <- stats::setNames(all.null[name.keepcoef,],name.keepcoef)
        }
    }else if(type %in% c("perm-var","perm-res")){
        null <- NULL
    }
    
    
    ## ** initialize
    ## *** data
    data <- object$data
    var.cluster <- object$cluster$var
    vec.Uid <- object$cluster$level
    n.cluster <- object$cluster$n
    n.obs <- NROW(data)
    index.cluster <- object$design$index.cluster
    Vindex.cluster <- unlist(index.cluster)
    for(iCluster in names(index.cluster)){
        names(index.cluster[[iCluster]]) <- rep(iCluster, length(index.cluster[[iCluster]]))
    }
    precompute.moments <- !is.null(object$design$precompute.XX)
  
    var.all <- stats::variable.names(object, effects = "all")
    var.mean <- stats::variable.names(object, effects = "mean")
    XX.all <- c("XXindexXX", "XXclusterXX", "XXcluster.indexXX", "XXtimeXX", "XXtime.indexXX", "XXstrataXX", "XXstrata.indexXX")

    var.weights <- object$weights$var
    var.outcome <- object$outcome$var

    ## *** parameters
    if(inherits(object$design$vcov,"ID") || inherits(object$design$vcov,"IND")){
        param.init <- NULL ## OLS solver + empirical residual standard deviation can directly find the solution
    }else{
        param.init <- object$param
        if(type == "perm-var"){
            param.init[name.keepcoef] <- 0
        }else if(type == "perm-res"){
            param.init[effects] <- 0
        }
    }
    ## ## *** missing patterns
    ## if(type == "perm-var"){
    ##     Upattern <- object$design$vcov$Upattern
    ##     ## find patterns sharing the same times 
    ##     ls.time.pattern <- unique(Upattern$time)  
    ##     time.pattern <- lapply(ls.time.pattern, paste, collapse="|")
    ##     n.patterns <- length(time.pattern)

    ##     ## gather clusters with patterns sharing the same times 
    ##     match.tp <- match(lapply(Upattern$time,paste,collapse="|"), time.pattern)
    ##     ls.indexcluster.pattern <- as.list(by(Upattern, match.tp, function(iData){unlist(iData$index.cluster)}, simplify = FALSE))
    ##     n.idpatterns <- lapply(ls.indexcluster.pattern,length)
    ## }

    ## *** residuals
    if(type == "perm-res"){
        ## re-estimate the model under the null hypothesis
        call0 <- object$call
        call0$formula <- stats::update(eval(call0$formula), paste0("~.-",effects))
        if(!is.null(object$index.na)){
            call0$data <- object$data.original[-object$index.na,,drop=FALSE]
        }
        object0 <- eval(call0)
        Xbeta0 <- stats::predict(object0, newdata = data[var.mean], se = FALSE, df = FALSE, simplify = TRUE)
        epsilon0.norm <- stats::residuals(object0, type = "normalized")
        OmegaChol0 <- lapply(stats::sigma(object0, chol = TRUE, cluster = as.character(vec.Uid)), FUN = base::t)
    }

    ## *** seed
    if(!is.null(seed)){

        if(length(seed)!=1 && length(seed) != n.sample){
         
            stop("Incorrect length for argument \'seed\': should either have length 1 or the number of simulations (here ",n.sample,"). \n",
                 "Current length: ",length(seed),". \n")

        }else if(length(seed)==n.sample){

            test.seed <- TRUE
            seqSeed <- seed

        }else if(length(seed)==1){

            tol.seed <- 10^(floor(log10(.Machine$integer.max))-1)
            if(n.sample>tol.seed){
                stop("Cannot set a seed per simulation when considering more than ",tol.seed," similations. \n")
            }
            set.seed(seed)
            test.seed <- TRUE
            seqSeed <- sample.int(tol.seed, n.sample,  replace = FALSE)

        }

        if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
            old <- .Random.seed # to save the current seed
            on.exit(try(.Random.seed <<- old, silent = TRUE)) # restore the current seed (before the call to the function)
        }else{
            on.exit(rm(.Random.seed, envir=.GlobalEnv))
        }
            
    }else{
        test.seed <- FALSE
        seqSeed <- NULL
    }


    ## ** function
    warperResample <- function(iSample){

        ## *** resample
        if(type == "perm-var"){
            ## permute X-values (constant within cluster)
            iPerm <- sample(n.cluster, replace = FALSE)
            iData <- data
            iData[unlist(index.cluster),effects] <- unlist(mapply(x = index.cluster, y = Uvar[iPerm], FUN = function(x,y){rep(y,length(x))}))
                        
            ## permute X-values between individuals with same missing data pattern
            ## iData <- data
            ## for(iPattern in 1:n.patterns){ ## iPattern <- 1
            ##     iPerm <- sample(n.idpatterns[[iPattern]], replace = FALSE)
            ##     iPatternindex.cluster <-  index.cluster[ls.indexcluster.pattern[[iPattern]]]
            ##     iData[[effects]][unlist(iPatternindex.cluster)] <- data[[effects]][unlist(iPatternindex.cluster[iPerm])]
            ## }
        }else if(type == "perm-res"){            
            iData <- data
            ## permute residuals (including the fixed effect to be tested)
            iPerm <- sample(1:n.obs, replace = FALSE)
            iEpsilon0.norm <- epsilon0.norm[iPerm]
            ## rescale residuals
            iData[[var.outcome]][unlist(index.cluster)] <- unlist(lapply(1:n.cluster, function(iC){OmegaChol0[[iC]] %*% iEpsilon0.norm[index.cluster[[iC]]]}))
            ## add fixed effects
            iData[[var.outcome]] <- iData[[var.outcome]] + Xbeta0[iPerm]
            ## range(iData[[var.outcome]] - data[,var.outcome])
        }else if(type == "boot"){
            iIndex.bootcluster <- index.cluster[sample(n.cluster, replace = TRUE)]
            iData <- data[unlist(iIndex.bootcluster),,drop=FALSE]
            iData[[var.cluster]] <- unlist(mapply(x=1:n.cluster, times=lengths(iIndex.bootcluster), FUN = rep, SIMPLIFY = FALSE))
            rownames(iData) <- NULL
            ## lmm(Y~X1*X2+X5,data = iData[,c("Y","X1","X2","X5")])
        }

        ## *** update design
        if(effects.vcov){
            ## update design according to the permutation (mean and variance)
            iStructure <- object$design$vcov
            iStructure$var <- NULL
            iStructure$cor <- NULL
            iStructure$param <- NULL
            iStructure$pattern <- NULL
            iStructure$Upattern <- NULL
            iStructure$pair.vcov <- NULL
            iStructure$pair.meanvcov <- NULL

            iData2 <- .lmmNormalizeData(iData[,var.all,drop=FALSE],
                                        var.outcome = var.outcome,
                                        var.cluster = attr(var.cluster, "original"),
                                        var.time = attr(object$time$var, "original"),
                                        var.strata = attr(object$strata$var, "original"),
                                        droplevels = TRUE,
                                        initialize.cluster = iStructure$ranef$crossed,
                                        initialize.time = setdiff(iStructure$ranef$vars, iStructure$var.cluster),
                                        na.rm = TRUE)$data

            iDesign <- .model.matrix.lmm(formula.mean = object$formula$mean.design,
                                         structure = iStructure,
                                         data = iData2,
                                         var.outcome = object$outcome$var,
                                         var.weights = object$weights$var,
                                         precompute.moments = precompute.moments,
                                         drop.X = object$design$drop.X,
                                         options = options)

        }else if(type == "perm-res"){ ## change in the Y values

            iDesign <- object$design
            iDesign$Y <- iData[[var.outcome]]

            ## update pre-computation              
            if(precompute.moments && NCOL(iDesign$mean)>0){
                if(is.na(var.weights[1])){
                    iwY <- cbind(iData[[var.outcome]])
                }else{
                    iwY <- cbind(iData[[var.outcome]]*sqrt(iData[[var.weights[1]]]))
                }
                iDesign$precompute.XY <- .precomputeXR(X = iDesign$precompute.XX$Xpattern, residuals = iwY, pattern = iDesign$vcov$Upattern$name,
                                                       pattern.ntime = stats::setNames(iDesign$vcov$Upattern$n.time, iDesign$vcov$Upattern$name),
                                                       pattern.cluster = attr(iDesign$vcov$pattern,"list"), index.cluster = iDesign$index.cluster)
            }

        }else{ ## change in the X values
            ## update design matrix according to the permutation (only mean)
            iDesign <- object$design
            iDesign$mean <- stats::model.matrix(object, newdata = iData[,var.mean,drop=FALSE])
            attr(iDesign$mean, "assign") <- attr(object$design$mean, "assign")
            attr(iDesign$mean, "contrasts") <- attr(object$design$mean, "contrasts")
            attr(iDesign$mean, "variable") <- attr(object$design$mean, "variable")
            attr(iDesign$mean, "terms") <- attr(object$design$mean, "terms")

            ## update pre-computation                
            if(precompute.moments && NCOL(iDesign$mean)>0){
                if(is.na(var.weights[1])){
                    iwX.mean <- iDesign$mean
                    iwY <- cbind(iData[[var.outcome]])
                }else{
                    iwX.mean <- sweep(iDesign$mean, FUN = "*", MARGIN = 1, STATS = sqrt(iData[[var.weights[1]]]))
                    iwY <- cbind(iData[[var.outcome]]*sqrt(iData[[var.weights[1]]]))
                }
                iIndex.cluster <- .extractIndexData(data = iData, structure = iDesign$vcov)$index.cluster
                iDesign$precompute.XX <- .precomputeXX(X = iwX.mean, pattern = iDesign$vcov$Upattern$name, 
                                                       pattern.ntime = stats::setNames(iDesign$vcov$Upattern$n.time, iDesign$vcov$Upattern$name),
                                                       pattern.cluster = attr(iDesign$vcov$pattern,"list"), index.cluster = iIndex.cluster)
                iDesign$precompute.XY <- .precomputeXR(X = iDesign$precompute.XX$Xpattern, residuals = iwY, pattern = iDesign$vcov$Upattern$name,
                                                       pattern.ntime = stats::setNames(iDesign$vcov$Upattern$n.time, iDesign$vcov$Upattern$name),
                                                       pattern.cluster = attr(iDesign$vcov$pattern,"list"), index.cluster = iIndex.cluster)
            }
        }

        ## *** re-estimate
        iEstimate <- try(.estimate(design = iDesign,
                                   time = object$time,
                                   method.fit = object$args$method.fit,
                                   type.information = object$args$type.information,
                                   transform.sigma = object$reparametrize$transform.sigma,
                                   transform.k = object$reparametrize$transform.k,
                                   transform.rho = object$reparametrize$transform.rho,
                                   precompute.moments = precompute.moments, 
                                   optimizer = object$args$control$optimizer,
                                   init = param.init,
                                   n.iter = object$opt$control["n.iter"],
                                   tol.score = object$opt$control["tol.score"],
                                   tol.param = object$opt$control["tol.param"],
                                   n.backtracking = object$opt$control["n.backtracking"],
                                   trace = FALSE), silent = TRUE)

        if(!inherits(iEstimate,"try-error") && studentized){
            iVcov <- .moments.lmm(value = iEstimate$estimate,
                                  design = iDesign,
                                  time = object$time,
                                  method.fit = object$args$method.fit,
                                  type.information = object$args$type.information,
                                  transform.sigma = "none", ## object$reparametrize$transform.sigma,
                                  transform.k = "none", ## object$reparametrize$transform.k,
                                  transform.rho = "none", ## object$reparametrize$transform.rho,
                                  logLik = FALSE, score = FALSE, information = FALSE, vcov = TRUE, df = FALSE, indiv = FALSE,
                                  effects = list("mean",c("mean","variance"))[[effects.vcov+1]], robust = FALSE,
                                  trace = FALSE, precompute.moments = precompute.moments, method.numDeriv = "simple", transform.names = FALSE)$vcov
            
        }

        ## *** export
        if(inherits(iEstimate,"try-error")){
            return(iEstimate)
        }else if(!studentized){
            if(effects.fct){
                return(c(convergence = iEstimate$cv, effects(iEstimate$estimate)))
            }else{
                return(c(convergence = iEstimate$cv, iEstimate$estimate[name.keepcoef]))
            }
        }else if(studentized){
            if(effects.fct){
                iRes <- effects(iEstimate$estimate, iVcov)
                return(c(convergence = iEstimate$cv, iRes[1,], se = iRes[2,]))
            }else{
                return(c(convergence = iEstimate$cv, iEstimate$estimate[name.keepcoef], se = sqrt(diag(iVcov)[name.keepcoef])))
            }
        }

    }
    ## res <- warperResample(1)
    ## ** run
    if(trace){
        if(type == "perm-var"){
            if(studentized){
                cat("\tStudentized permutation test with ",n.sample," replicates \n",
                    "\t(permutation of the regressor values between clusters)\n\n", sep = "")
            }else{
                cat("\tPermutation test with ",n.sample," replicates \n",
                    "\t(permutation of the regressor values between clusters)\n\n", sep = "")
            }
        }else if(type == "perm-res"){
            if(studentized){
                cat("\tStudentized permutation test with ",n.sample," replicates \n",
                    "\t(permutation of the normalized residuals under the null)\n\n", sep = "")
            }else{
                cat("\tPermutation test with ",n.sample," replicates \n",
                    "\t(permutation of the normalized residuals under the null)\n\n", sep = "")
            }
        }else if(type == "boot"){
            if(studentized){
                cat("\tNon-parametric studentized bootstrap with ",n.sample," replicates \n\n", sep = "")
            }else{
                cat("\tNon-parametric bootstrap with ",n.sample," replicates \n\n", sep = "")
            }
        }
    }
    
    
    if(cpus==1){
        if (trace > 0 & requireNamespace("pbapply")) {
            method.loop <- pbapply::pblapply
        }else{
            method.loop <- lapply
        }
        ls.sample <- do.call(method.loop,
                             args = list(X = 1:n.sample,
                                         FUN = function(iX){
                                             if(test.seed){set.seed(seqSeed[iX])}
                                             iOut <- warperResample(iX)
                                             if(!is.null(seed) && !inherits(iOut,"try-error")){
                                                 return(c(sample = iX, seed = seqSeed[iX],iOut))
                                             }else{
                                                 return(c(sample = iX, iOut))
                                             }
                                         })
                             )
    }else if(cpus>1){
        ## split into a 100 jobs
        split.resampling <- parallel::splitIndices(nx = n.sample, ncl = min(max(100,10*cpus), n.sample))
        nsplit.resampling <- length(split.resampling)

        ## define cluster
        if(trace>1){
            ## display all output from each core
            cl <- parallel::makeCluster(cpus, outfile = "")
        }else{
            cl <- parallel::makeCluster(cpus)
        }

        ## progress bar
        if(trace>0){
            pb <- utils::txtProgressBar(max = nsplit.resampling, style = 3)          
            progress <- function(n){utils::setTxtProgressBar(pb, n)}
            opts <- list(progress = progress)
        }else{
            opts <- list()
        }

        ## link to foreach
        doSNOW::registerDoSNOW(cl)
        ## to export from current environment and private functions from the package
        LMMstar.fct <- ls(asNamespace("LMMstar"), all.names = TRUE)

        toExport <- c("addLeading0",".augmodel.matrix",
                      grep("^collapse",LMMstar.fct, value = TRUE),
                      ".estimate",".extractIndexData",
                      grep("^\\.findUpatterns",LMMstar.fct, value = TRUE),
                      ".model.matrix.lmm",".model.matrix_regularize",".lmmNormalizeData",".moments.lmm",
                      ".precomputeXR",".precomputeXX",
                      grep("^\\.skeleton",LMMstar.fct, value = TRUE),
                      "unorderedPairs",
                      ".vcov.matrix.lmm")
        if(!is.null(export.cpus)){
            toExport <- c(toExport, "export.cpus")
        } 

        ## run
        iBlock <- NULL ## [:forCRANcheck:] foreach        a
        ls2.sample <- foreach::`%dopar%`(
                                   foreach::foreach(iBlock=1:nsplit.resampling,
                                                    .export = toExport,
                                                    .packages = c("LMMstar","nlme"),
                                                    .options.snow = opts), {

                                                       iOut <- lapply(split.resampling[[iBlock]], function(iSplit){
                                                           if(!is.null(seed)){set.seed(seqSeed[iSplit])}
                                                           iOut <- warperResample(iSplit)                                                       
                                                           if(!is.null(seed) && !inherits(iOut,"try-error")){
                                                               return(c(sample = iSplit, seed = seqSeed[iSplit],iOut))
                                                           }else{
                                                               return(c(sample = iSplit, iOut))
                                                           }
                                                       })
                                                       return(iOut)

                                                   })
        if(trace>0){close(pb)}
        parallel::stopCluster(cl)

        ## collect
        ls.sample <- do.call("c",ls2.sample)
    }


    ## ** post-process
    index.noerror <- which(sapply(ls.sample, inherits, "try-error")==FALSE)
    M.sample <- do.call(rbind,ls.sample[index.noerror])
    
    out <- list(call = match.call(),
                args = list(type = type, effects = effects, n.sample = n.sample, studentized = studentized, seed = seed))

    if(effects.fct){
        if(studentized){
            out$estimate <- effects.estimate[1,]
            out$se <- effects.estimate[2,]
        }else{
            out$estimate <- effects.estimate
        }
    }else{
        out$estimate <- value.meancoef[name.keepcoef]
        out$se <- sd.meancoef[name.keepcoef]
    }

    out$sample.estimate <- matrix(NA, nrow = n.sample, ncol = length(name.keepcoef),
                                  dimnames = list(NULL,name.keepcoef))
    out$sample.estimate[M.sample[,"sample"],] <- M.sample[,name.keepcoef,drop=FALSE]

    out$null <- null

    if(studentized){
        out$sample.se <- matrix(NA, nrow = n.sample, ncol = length(name.keepcoef),
                                dimnames = list(NULL,name.keepcoef))
        out$sample.se[M.sample[,"sample"],] <- M.sample[,paste0("se.",name.keepcoef),drop=FALSE]
    }

    out$cv <- rep(FALSE, length = n.sample)
    out$cv[M.sample[,"sample"]] <- M.sample[,"convergence"]

    if(!is.null(seqSeed)){
        out$seed <- seqSeed
    }

    ## ** export
    class(out) <- append("resample",class(out))
    return(out)
}

## * resample.mlmm
##' @export
##' @rdname resample
resample.mlmm <- function(object, type, method = NULL, cluster = NULL, n.sample = 1e3, studentized = TRUE,
                          trace = TRUE, seed = NULL, cpus = 1, export.cpus = NULL,
                          ...){

    options <- LMMstar.options()
    pool.method <- options$pool.method

    ## ** normalize arguments
    ## alternative
    if(object$glht[[1]][[1]]$alternative!="two.sided"){
        stop("Can only perform two sided tests. \n")
    }
    
    ## n.sample
    if(length(n.sample)!=1){
        stop("Argument \'n.sample\' should have lenght 1. \n")
    }
    if(!is.numeric(n.sample) || n.sample<0 || n.sample %% 1 >0 ){
        stop("Argument \'n.sample\' should be a positive integer. \n")
    }
    
    ## cluster
    if(is.null(cluster)){
        ls.nameCluster <- lapply(object$model, function(iModel){attr(variable.names(iModel),"cluster")})
        if(length(unique(unlist(ls.nameCluster)))==1 && all(lengths(ls.nameCluster)==1)){
            cluster <- unname(ls.nameCluster[[1]])
        }else{
            stop("Argument \'cluster\' cannot be guessed from the object. \n")
        }
        
    }
    if(length(cluster)!=1){
        stop("Argument \'cluster\' should have lenght 1. \n")
    }
    if(!is.character(cluster)){
        stop("Argument \'cluster\' should be a character. \n")
    }

    ## data    
    object.call <- attr(object, "call")
    by <- object.call$by
    level.by <- names(object$model)
    data <- try(eval(object.call$data), silent = TRUE)

    if(inherits(data, "try-error")){ ## reconstruct dataset based on what was stored in the object
        object.model <- object$model
        n.model <- length(object.model)
        object.manifest <- lapply(object.model, manifest)
        object.data <- lapply(object.model, "[[", "data")

        if(any(sapply(object.data, function(iData){cluster %in% names(iData)})==FALSE)){
            stop("Argument \'cluster\' could not be found in the dataset(s) stored in the object. \n")
        }
        object.Udata <- mapply(FUN = function(iData,iName,iBy){
            iOut <- cbind(iBy,iData[c(cluster,iName)])
            names(iOut)[1] <- by
            return(iOut)
        }, iData = object.data, iName = object.manifest, iBy = level.by, SIMPLIFY = FALSE)
        data <- object.Udata[[1]]

        if(n.model>1){
            for(iModel in 2:n.model){ ## iModel <- 3
                data <- merge(data, object.Udata[[iModel]], by = intersect(names(data),c(object.manifest[[iModel]],by,cluster)), all=TRUE)
            }
        }
    
    }else{
        if(cluster %in% names(data)==FALSE){
            stop("Argument \'cluster\' could not be found in the dataset. \n")
        }
    }
    data <- as.data.frame(data)

    ## type
    if(tolower(type) == "permutation"){
        type <- "perm"
    }else if(tolower(type) == "bootstrap"){
        type <- "boot"
    }
    type <- match.arg(type, c("perm","boot"))
    if(type=="perm"){

        effect.param <- object$univariate$parameter
        if(any(object$univariate$type!="mu")){
            stop("Can only test mean parameters when using permutations. \n")
        }
        ls.variable.perm <- mapply(FUN = function(iParam,iObject){ ## iParam <- "X2b" ## iObject <- object$model[[1]]
            iX <- iObject$design$mean 
            attr(iX,"variable")[attr(iX,"assign")[colnames(iX)==iParam]]
        }, iParam = effect.param, iObject = object$model, SIMPLIFY = FALSE)
        variable.perm <- unname(unique(unlist(ls.variable.perm)))
        if(length(variable.perm)>1){
            stop("Permutation test not available for effects depending on more than one variable. \n",
                 "Can only permute a single variable. \n")
        }

        ## test whether the covariate is constant within cluster
        test.Wvar <- unlist(by(data[[variable.perm]], data[[cluster]], function(iData){sum(!duplicated(iData))>1}, simplify = FALSE))
        if(any(test.Wvar)){
            stop("The covariate(s) indicated by the argument \"effects\" vary within cluster. \n",
                 "This is not support when using permutations. Consider using type=\"boot\" instead of type=\"perm\". \n")
        }
        Uvar <- by(data[[variable.perm]], data[[cluster]], function(iData){iData[1]}, simplify = FALSE)
    }

    ## cpus
    max.cpus <- parallel::detectCores()
    if(length(cpus)!=1){
        stop("Argument \'cpus\' should have length 1.\n ")
    }else if(identical(cpus,"all")){
        cpus <- max.cpus
    }else if(!is.numeric(cpus) || cpus <=0 || cpus %% 1 != 0){
        stop("Argument \'cpus\' should be an integer between 1 and ",max.cpus," or \'all\'.\n ")
    }else if(cpus>1 && cpus>parallel::detectCores()){
        stop("Argument \'cpus\' exceeds the number of available CPUs.\n ",
             "It should be an integer between 1 and ",max.cpus," or \'all\'.\n ")
    }

    ## dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** prepare
    ## *** seed
    if(!is.null(seed)){

        if(length(seed)!=1 && length(seed) != n.sample){
         
            stop("Incorrect length for argument \'seed\': should either have length 1 or the number of simulations (here ",n.sample,"). \n",
                 "Current length: ",length(seed),". \n")

        }else if(length(seed)==n.sample){

            test.seed <- TRUE
            seqSeed <- seed

        }else if(length(seed)==1){

            tol.seed <- 10^(floor(log10(.Machine$integer.max))-1)
            if(n.sample>tol.seed){
                stop("Cannot set a seed per simulation when considering more than ",tol.seed," similations. \n")
            }
            set.seed(seed)
            test.seed <- TRUE
            seqSeed <- sample.int(tol.seed, n.sample,  replace = FALSE)

        }

        if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
            old <- .Random.seed # to save the current seed
            on.exit(try(.Random.seed <<- old, silent = TRUE)) # restore the current seed (before the call to the function)
        }else{
            on.exit(rm(.Random.seed, envir=.GlobalEnv))
        }
            
    }else{
        test.seed <- FALSE
        seqSeed <- NULL
    }

    ## *** recover dataset
    keep.cols <- c("by","parameter","estimate","se","df","statistic","lower","upper","null","p.value")
        
    Ucluster <- sort(unique(unlist(data[[cluster]])))
    n.cluster <- length(Ucluster)
    index.cluster <- tapply(1:NROW(data), data[[cluster]], identity)

    ## ** function
    warperResample <- function(iSample){

        ## *** resample
        if(type == "perm"){
            ## permute X-values (constant within cluster)
            iPerm <- sample(n.cluster, replace = FALSE)
            iData <- data
            for(iCluster in 1:n.cluster){ ## iCluster <- 1
                iData[index.cluster[[iCluster]],variable.perm] <- Uvar[[iPerm[iCluster]]] ## seems to properly expand Uvar over multiple timepoints
            }
        }else if(type == "boot"){
            ils.Boot <- index.cluster[sample(n.cluster, replace = TRUE)]
            iBoot <- unlist(ils.Boot)
            iBoot.cl <- unlist(mapply(x = 1:n.cluster, y = lengths(ils.Boot), function(x,y){rep(x,y)}, SIMPLIFY = FALSE))

            iData <- data[iBoot,,drop=FALSE]
            iData[[cluster]] <- iBoot.cl
            rownames(iData) <- NULL            
        }

        ## *** re-estimate
        iCall <- object.call
        iCall$trace <- 0
        iCall$data <- iData
        if(all(method %in% "p.rejection" == FALSE)){
            iCall$df <- FALSE
        }
        iEstimate <- try(eval(iCall), silent = TRUE)

        ## *** export
        ## use merge in case effects for some categories are not estimated
        if(inherits(iEstimate,"try-error")){
            iOut <- iEstimate
        }else{
            mean.cv <- mean(sapply(iEstimate$model,function(iM){iM$opt$cv}))
            if(studentized){
                iTable <- model.tables(iEstimate, columns = keep.cols, method = method, df = FALSE)
                iOut <- c(convergence =  mean.cv, iTable$estimate, iTable$se)
                names(iOut)[-1] <- paste(rownames(iTable),paste0("se.",rownames(iTable)))
            }else{
                iCoef <- coef(iEstimate, method = method)
                iOut <- c(convergence =  mean.cv, iCoef)
            }
        }
        return(iOut)
    }

    ## ** run
    if(trace){
        if(type == "permutation"){
            cat("\tPermutation test with ",n.sample," replicates \n",
                "\t(permutation of ",variable.perm," values between clusters)\n\n", sep = "")
        }else if(type == "boot"){
            cat("\tNon-parametric bootstrap with ",n.sample," replicates \n\n", sep = "")
        }
    }
    


    if(cpus==1){
        if (trace > 0 & requireNamespace("pbapply")) {
            method.loop <- pbapply::pblapply
        }else{
            method.loop <- lapply
        }

        ls.sample <- do.call(method.loop,
                             args = list(X = 1:n.sample,
                                         FUN = function(iX){ ## iX <- 1
                                             if(test.seed){set.seed(seqSeed[iX])}
                                             iOut <- warperResample(iX)
                                             if(!is.null(seed) && !inherits(iOut,"try-error")){
                                                 return(c(sample = iX, seed = seqSeed[iX],iOut))
                                             }else{
                                                 return(c(sample = iX, iOut))
                                             }
                                         })
                             )
    }else if(cpus>1){
        ## split into a 100 jobs
        split.resampling <- parallel::splitIndices(nx = n.sample, ncl = min(max(100,10*cpus), n.sample))
        nsplit.resampling <- length(split.resampling)

        ## define cluster
        if(trace>1){
            ## display all output from each core
            cl <- parallel::makeCluster(cpus, outfile = "")
        }else{
            cl <- parallel::makeCluster(cpus)
        }

        ## progress bar
        if(trace>0){
            pb <- utils::txtProgressBar(max = nsplit.resampling, style = 3)          
            progress <- function(n){utils::setTxtProgressBar(pb, n)}
            opts <- list(progress = progress)
        }else{
            opts <- list()
        }

        ## link to foreach
        doSNOW::registerDoSNOW(cl)
        ## to export from current environment and private functions from the package
        toExport <- NULL
        if(!is.null(export.cpus)){
            toExport <- c(toExport, "export.cpus")
        } 

        ## run
        iBlock <- NULL ## [:forCRANcheck:] foreach        a
        ls.sample <- foreach::`%dopar%`(
                                  foreach::foreach(iBlock=1:n.split.resampling,
                                                   .export = toExport,
                                                   .packages = c("LMMstar","nlme"),
                                                   .options.snow = opts), {

                                                       iOut <- lapply(split.resampling[[iBlock]], function(iSplit){
                                                           if(!is.null(seed)){set.seed(seqSeed[iSplit])}
                                                           iOut <- warperResample(iSplit)
                                                           if(!is.null(seed) && !inherits(iOut,"try-error")){
                                                               return(c(sample = iSplit, seed = seqSeed[iSplit],iOut))
                                                           }else{
                                                               return(c(sample = iSplit, iOut))
                                                           }
                                                       })
                                                       return(iOut)

                                                   })
        if(trace>0){close(pb)}
        parallel::stopCluster(cl)

        ## collect
        ls.sample <- do.call("c",ls2.sample)
    }
    
    ## ** post process
    index.noerror <- which(sapply(ls.sample, inherits, "try-error")==FALSE)
    M.sample <- do.call(rbind,ls.sample[index.noerror])

    out <- list(call = match.call(),
                args = list(type = type, method = method, n.sample = n.sample, studentized = studentized, seed = seed))

    if(studentized){
        object.table <- model.tables(object, columns = c("estimate","se","null"), method = method, df = FALSE)
        out$estimate <- stats::setNames(object.table$estimate, object.table$parameter)
        out$se <- stats::setNames(object.table$se, object.table$parameter)
        out$null <- stats::setNames(object.table$null, object.table$parameter)
    }else{
        object.table <- model.tables(object, columns = c("estimate","null"), method = method, df = FALSE)
        out$estimate <- stats::setNames(object.table$estimate, object.table$parameter)
        out$null <- stats::setNames(object.table$null, object.table$parameter)
    }

    name.keepcoef <- names(out$estimate)
    out$sample.estimate <- matrix(NA, nrow = n.sample, ncol = length(name.keepcoef),
                                  dimnames = list(NULL,name.keepcoef))
    out$sample.estimate[M.sample[,"sample"],] <- M.sample[,name.keepcoef,drop=FALSE]

    if(studentized){
        out$sample.se <- matrix(NA, nrow = n.sample, ncol = length(name.keepcoef),
                                dimnames = list(NULL,name.keepcoef))
        out$sample.se[M.sample[,"sample"],] <- M.sample[,paste0("se.",name.keepcoef),drop=FALSE]
    }

    out$cv <- rep(FALSE, length = n.sample)
    out$cv[M.sample[,"sample"]] <- M.sample[,"convergence"]

    if(!is.null(seqSeed)){
        out$seed <- seqSeed
    }

    ## ** export
    class(out) <- append("resample",class(out))
    return(out)
}


##----------------------------------------------------------------------
### resample.R ends here
