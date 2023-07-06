### resample.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 31 2022 (10:09) 
## Version: 
## Last-Updated: jul  4 2023 (18:58) 
##           By: Brice Ozenne
##     Update #: 438
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
##' @param n.sample [integer] the number of samples used.
##' @param studentized [logical] should a studentized boostrap or permutation test be used?
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
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
##' resample(eCS.lmm, type = "boot", effects = c("glucagonAUC","gender"), seed = 10)
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
                         level = 0.95, correction = TRUE,
                         trace = TRUE, seed = NULL, cpus = 1, export.cpus = NULL,
                         ...){

    ## ** check user input
    alpha <- 1-level
    type <- match.arg(type, c("perm-var","perm-res","boot"))
    ## boot: non-parametric bootstrap
    ## perm-var: exchange covariate value between clusters (need to be constant within cluster)
    ## perm-res: outcome becomes permuted normalized residuals under the null to which the (permuted) fixed effect under the null are added.
    if(!is.logical(correction) || correction %in% 0:1 == FALSE){
        stop("Argument \'correction\' must be binary or logical. \n")
    }

    name.meanvar <- attr(object$design$mean,"variable")
    value.meancoef <- coef(object, effects = "mean")
    sd.meancoef <- sqrt(diag(vcov(object, effects = "mean")))

    name.meancoef <- names(value.meancoef)
    if(type == "perm-var"){
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
        effects.vcov <- any(effects %in% manifest(object, effects = c("variance","correlation")))
        
        ## test whether the covariate is constant within cluster
        test.Wvar <- unlist(by(object$data[effects], object$data[[object$cluster$var]], function(iData){sum(!duplicated(iData))>1}, simplify = FALSE))
        if(any(test.Wvar)){
            stop("The covariate(s) indicated by the argument \"effects\" vary within cluster. \n",
                 "This is not support when using permutations. Consider using type=\"boot\" instead of type=\"perm\". \n")
        }
        Uvar <- by(object$data[effects], object$data[[object$cluster$var]], function(iData){iData[1,,drop=FALSE]}, simplify = FALSE)

    }else if(type == "perm-res"){
        if(any(is.character(effects)==FALSE)){
            stop("Argument \'effects\' should contain character strings refering to a parameter of the mean structure. \n")
        }
        if(any(effects %in% name.meancoef == FALSE)){
            stop("Argument \'effects\' should be one of \"",paste(name.meancoef, collapse = "\", \""),"\". \n")
        }
        effects.vcov <- FALSE
        name.keepcoef <- effects
    }else if(type == "boot"){
        if(any(is.character(effects)==FALSE)){
            stop("Argument \'effects\' should contain character strings refering to a parameter of the mean structure. \n")
        }
        if(any(effects %in% name.meancoef == FALSE)){
            stop("Argument \'effects\' should be one of \"",paste(name.meancoef, collapse = "\", \""),"\". \n")
        }
        effects.vcov <- TRUE
        name.keepcoef <- effects
    }

    if(length(cpus)!=1){
        stop("Argument \'cpus\' should have length 1.\n ")
    }else if(identical(cpus,"all")){
        cpus <- parallel::detectCores()
    }else if(!is.numeric(cpus) || cpus <=0 || cpus %% 1 != 0){
        stop("Argument \'cpus\' should be a positive integer or \'all\'.\n ")
    }else if(cpus>1 && cpus>parallel::detectCores()){
        stop("Only ",parallel::detectCores()," CPU cores are available. \n")
    }

    ## dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** initialize
    ## *** data
    data <- object$data
    if(length(object$index.na)>0){
        data <- data[-object$index.na,]
    }
    var.cluster <- object$cluster$var
    vec.Uid <- object$cluster$level
    n.cluster <- object$cluster$n
    n.obs <- NROW(data)
    index.cluster <- object$design$index.cluster
    Vindex.cluster <- unlist(index.cluster)
    for(iCluster in names(index.cluster)){
        names(index.cluster[[iCluster]]) <- rep(iCluster, length(index.cluster[[iCluster]]))
    }
    nTime.cluster <- lapply(index.cluster,length)
    precompute.moments <- !is.null(object$design$precompute.XX)
  
    var.all <- manifest(object, effects = "all")
    var.mean <- manifest(object, effects = "mean")
    XX.all <- c("XXindexXX", "XXclusterXX", "XXcluster.indexXX", "XXtimeXX", "XXtime.indexXX", "XXstrataXX", "XXstrata.indexXX")

    var.weights <- object$weights$var
    var.outcome <- object$outcome$var

    ## *** parameters
    param.init <- object$param
    if(type == "perm-var"){
        param.init[name.keepcoef] <- 0
    }else if(type == "perm-res"){
        param.init[effects] <- 0
    }

    ## ## *** missing patterns
    ## if(type == "perm-var"){
    ##     Upattern <- object$design$vcov$X$Upattern
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
        Xbeta0 <- stats::predict(object0, newdata = data[,var.mean,drop=FALSE], se = FALSE, df = FALSE)[,1]
        epsilon0.norm <- stats::residuals(object0, type = "normalized")
        OmegaChol0 <- lapply(stats::sigma(object0, chol = TRUE, cluster = as.character(vec.Uid)), FUN = base::t)
    }

    ## ** function
    warperResample <- function(iSample){

        ## *** resample
        if(type == "perm-var"){
            ## permute X-values (constant within cluster)
            iPerm <- sample(n.cluster, replace = FALSE)
            iData <- data
            for(iCluster in 1:n.cluster){ ## iCluster <- 1
                iData[index.cluster[[iCluster]],effects] <- Uvar[[iPerm[iCluster]]] ## seems to properly expand Uvar over multiple timepoints
            }
            
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
            ils.Boot <- index.cluster[sample(n.cluster, replace = TRUE)]
            names(ils.Boot) <- 1:n.cluster
            iBoot <- unlist(ils.Boot)
            
            iData <- data[iBoot,,drop=FALSE]
            iData[[var.cluster]] <- as.numeric(factor(names(iBoot), levels = unique(names(iBoot))))
            rownames(iData) <- NULL
            ## lmm(Y~X1*X2+X5,data = iData[,c("Y","X1","X2","X5")])
        }

        ## *** update design
        if(effects.vcov){
            ## update design according to the permutation (mean and variance)
            iStructure <- object$design$vcov
            iStructure$X <- NULL
            iStructure$param <- NULL

            iData2 <- .prepareData(iData[,var.all,drop=FALSE],
                                   var.cluster = attr(var.cluster, "original"),
                                   var.time = attr(object$time$var, "original"),
                                   var.strata = attr(object$strata$var, "original"),
                                   missing.repetition = object$time$n>1,
                                   droplevels = TRUE)
            
            iDesign <- .model.matrix.lmm(formula.mean = object$formula$mean.design,
                                         structure = iStructure,
                                         data = iData2,
                                         var.outcome = object$outcome$var,
                                         var.weights = object$weights$var,
                                         precompute.moments = precompute.moments,
                                         drop.X = object$design$drop.X)
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
                iDesign$precompute.XY <- .precomputeXR(X = iDesign$precompute.XX$Xpattern, residuals = iwY, pattern = iDesign$vcov$X$Upattern$name,
                                                       pattern.ntime = stats::setNames(iDesign$vcov$X$Upattern$n.time, iDesign$vcov$X$Upattern$name),
                                                       pattern.cluster = iDesign$vcov$X$Upattern$index.cluster, index.cluster = iDesign$index.cluster)
            }

        }else{ ## change in the X values
            ## update design matrix according to the permutation (only mean)
            iDesign <- object$design
            iDesign$mean <- model.matrix(object, data = iData[,var.mean,drop=FALSE])
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
                iDesign$precompute.XX <- .precomputeXX(X = iwX.mean, pattern = iDesign$vcov$X$Upattern$name, 
                                                       pattern.ntime = stats::setNames(iDesign$vcov$X$Upattern$n.time, iDesign$vcov$X$Upattern$name),
                                                       pattern.cluster = iDesign$vcov$X$Upattern$index.cluster, index.cluster = iIndex.cluster)

                iDesign$precompute.XY <- .precomputeXR(X = iDesign$precompute.XX$Xpattern, residuals = iwY, pattern = iDesign$vcov$X$Upattern$name,
                                                       pattern.ntime = stats::setNames(iDesign$vcov$X$Upattern$n.time, iDesign$vcov$X$Upattern$name),
                                                       pattern.cluster = iDesign$vcov$X$Upattern$index.cluster, index.cluster = iIndex.cluster)
            }
        }

        ## *** re-estimate
        iEstimate <- try(.estimate(design = iDesign,
                                   time = object$time,
                                   method.fit = object$method.fit,
                                   type.information = attr(object$information,"type.information"),
                                   transform.sigma = object$reparametrize$transform.sigma,
                                   transform.k = object$reparametrize$transform.k,
                                   transform.rho = object$reparametrize$transform.rho,
                                   precompute.moments = precompute.moments, 
                                   optimizer = object$opt$name,
                                   init = param.init,
                                   n.iter = object$opt$control["n.iter"],
                                   tol.score = object$opt$control["tol.score"],
                                   tol.param = object$opt$control["tol.param"],
                                   n.backtracking = object$opt$control["n.backtracking"],
                                   trace = FALSE))

        if(!inherits(iEstimate,"try-error") && studentized){
            iVcov <- .moments.lmm(value = iEstimate$estimate,
                                  design = iDesign,
                                  time = object$time,
                                  method.fit = object$method.fit,
                                  type.information = attr(object$information,"type.information"),
                                  transform.sigma = object$reparametrize$transform.sigma,
                                  transform.k = object$reparametrize$transform.k,
                                  transform.rho = object$reparametrize$transform.rho,
                                  logLik = FALSE, score = FALSE, information = FALSE, vcov = TRUE, df = FALSE, indiv = FALSE, effects = c("mean"), robust = FALSE,
                                  trace = FALSE, precompute.moments = precompute.moments, method.numDeriv = "simple", transform.names = FALSE)$vcov
            
        }

        ## *** export
        if(inherits(iEstimate,"try-error")){
            return(iEstimate)
        }else if(!studentized){
            return(c(iEstimate$cv, iEstimate$estimate[name.keepcoef]))
        }else if(studentized){
            c(iEstimate$cv, iEstimate$estimate[name.keepcoef], se = sqrt(diag(iVcov)[name.keepcoef]))
        }

    }
    ## res <- warperResample(1)
    ## ** run
    if(trace){
        if(type == "perm-var"){
            if(studentized){
                cat("\tStudentized permutation test with ",n.sample," replicates \n",
                    "\tn(permutation of the regressor values between clusters)\n\n", sep = "")
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
    
    if(!is.null(seed)){
        tol.seed <- 10^(floor(log10(.Machine$integer.max))-1)
        if(n.sample>tol.seed){
            stop("Cannot set a seed per simulation when considering more than ",tol.seed," similations. \n")
        }
        if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
            old <- .Random.seed # to save the current seed
            on.exit(.Random.seed <<- old) # restore the current seed (before the call to the function)
        }else{
            on.exit(rm(.Random.seed, envir=.GlobalEnv))
        }
        set.seed(seed)
        seqSeed <- sample.int(tol.seed, n.sample,  replace = FALSE)        
    }else{
        seqSeed <- NULL
    }

    if(cpus==1){
        if (trace > 0) {
            requireNamespace("pbapply")
            method.loop <- pbapply::pblapply
        }else{
            method.loop <- lapply
        }
    }else if(cpus>1){
        cl <- parallel::makeCluster(cpus)
        ## link to foreach
        doSNOW::registerDoSNOW(cl)
        ## export from user
        if(!is.null(export.cpus)){
            parallel::clusterExport(cl, export.cpus)
        } 
        ## export seed 
        if (!is.null(seed)) {
            parallel::clusterExport(cl, varlist = "seqSeed", envir = environment())
        }
        ## export BuyseTest 
        fct2export <- c(".estimate",".precomputeXR",".precomputeXX",".extractIndexData",".model.matrix.lmm",".prepareData")
        parallel::clusterExport(cl = cl, varlist = fct2export, envir = as.environment(asNamespace("LMMstar")))
        
    }

    if (cpus == 1) {
        ls.sample <- do.call(method.loop,
                             args = list(X = 1:n.sample,
                                         FUN = function(X){
                                             if(!is.null(seed)){set.seed(seqSeed[X])}
                                             iOut <- warperResample(X)
                                             return(iOut)
                                         })
                             )
    }else if(cpus > 1){
        if(trace>0){
            pb <- utils::txtProgressBar(max = n.sample, style = 3)          
            progress <- function(n){utils::setTxtProgressBar(pb, n)}
            opts <- list(progress = progress)
        }else{
            opts <- list()
        }

        ls.sample <- foreach::`%dopar%`(
                                  foreach::foreach(i=1:n.sample, .options.snow = opts), {
                                      if(!is.null(seed)){set.seed(seqSeed[i])}
                                      iOut <- warperResample(X)
                                      return(iOut)
                                  })
        if(trace>0){close(pb)}
        parallel::stopCluster(cl)
    }


    ## ** post-process
    M.sample <- do.call(rbind,ls.sample[sapply(ls.sample, inherits, "try-error")==FALSE])
    Mcv.sample <- M.sample[M.sample[,1]==1,-1,drop=FALSE]
    ncv.sample <- NROW(Mcv.sample)

    if(studentized){
        Mcv.sample.se <- Mcv.sample[,(length(name.keepcoef)+1):(2*length(name.keepcoef)),drop=FALSE]
        Mcv.sample <- Mcv.sample[,1:length(name.keepcoef),drop=FALSE]
    }
    
    out <- model.tables(object)[name.keepcoef,,drop=FALSE]*NA
    out$estimate <- as.double(value.meancoef[name.keepcoef])

    if(type %in% c("perm-var","perm-res") ){

        if(studentized){
            Mcv.estimate <- matrix(value.meancoef[name.keepcoef]/sd.meancoef[name.keepcoef], nrow = ncv.sample, ncol = length(name.keepcoef), byrow = TRUE,
                                   dimnames = list(NULL,name.keepcoef))
            Mcv.sample.Wald <- Mcv.sample/Mcv.sample.se

            out$se <- as.double(sd.meancoef[name.keepcoef])
            out$p.value <- (colSums(abs(Mcv.sample.Wald) > abs(Mcv.estimate), na.rm = TRUE)+correction)/(colSums(!is.na(Mcv.sample.Wald))+correction)
        }else{
            Mcv.estimate <- matrix(value.meancoef[name.keepcoef], nrow = ncv.sample, ncol = length(name.keepcoef), byrow = TRUE,
                                   dimnames = list(NULL,name.keepcoef))

            out$p.value <- (colSums(abs(Mcv.sample) > abs(Mcv.estimate), na.rm = TRUE)+correction)/(colSums(!is.na(Mcv.sample))+correction)
        }
    }else if(type == "boot"){
        if(studentized){
            Mcv.sample.Wald0 <- apply(Mcv.sample, MARGIN = 2, FUN = scale, scale = FALSE, center = TRUE)/Mcv.sample.se  ## center around the null
            
            out$se <- as.double(sd.meancoef[name.keepcoef])
            out$lower <- out$estimate + out$se * apply(Mcv.sample.Wald0, MARGIN = 2, FUN = stats::quantile, probs = alpha/2, na.rm = TRUE)
            out$upper <- out$estimate + out$se * apply(Mcv.sample.Wald0, MARGIN = 2, FUN = stats::quantile, probs = 1-alpha/2, na.rm = TRUE)
            out$p.value <- sapply(name.keepcoef, function(iName){
                boot2pvalue(stats::na.omit(out[iName,"estimate"] + out[iName,"se"] * Mcv.sample.Wald0[,iName]),
                            null = 0,
                            estimate = out[iName,"estimate"],
                            alternative = "two.sided",
                            add.1 = correction)
            })
            
        }else{
            Mcv.estimate <- matrix(value.meancoef[name.keepcoef], nrow = ncv.sample, ncol = length(name.keepcoef), byrow = TRUE,
                                   dimnames = list(NULL,name.keepcoef))

            out$se <- apply(Mcv.sample, MARGIN = 2, FUN = stats::sd, na.rm = TRUE)
            out$lower <- apply(Mcv.sample, MARGIN = 2, FUN = stats::quantile, probs = alpha/2, na.rm = TRUE)
            out$upper <- apply(Mcv.sample, MARGIN = 2, FUN = stats::quantile, probs = 1-alpha/2, na.rm = TRUE)
            out$p.value <- sapply(name.keepcoef, function(iName){
                boot2pvalue(stats::na.omit(Mcv.sample[,iName]),
                            null = 0,
                            estimate = out[iName,"estimate"],
                            alternative = "two.sided",
                            add.1 = correction)
            })
        }
    }

    ## ** export
    attr(out,"call") <- match.call()
    attr(out,"args") <- list(type = type, effects = effects, n.sample = n.sample, studentized = studentized, seed = seed)
    attr(out,"M.sample") <- M.sample
    attr(out,"n.sample") <- n.sample
    class(out) <- append("resample",class(out))
    return(out)

}

## * resample.mlmm
##' @export
##' @rdname resample
resample.mlmm <- function(object, type, method = NULL, cluster = NULL, n.sample = 1e3,
                          level = 0.95, correction = TRUE,
                          trace = TRUE, seed = NULL, cpus = 1,
                          ...){

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
        ls.nameCluster <- lapply(object$model, function(iModel){manifest(iModel)["cluster"]})
        if(length(unique(unlist(ls.nameCluster)))==1 && all(sapply(ls.nameCluster,length)==1)){
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

    ## type
    if(tolower(type) == "perm"){
        type <- "permutation"
    }else if(tolower(type) == "boot"){
        type <- "bootstrap"
    }
    type <- match.arg(type, c("permutation","bootstrap"))
    if(type=="permutation"){

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
    }else if(type=="bootstrap"){
        calc.pvalue <- try(requireNamespace("BuyseTest")) ## to get the p-values
    }
    
    ## cpus
    if(length(cpus)!=1){
        stop("Argument \'cpus\' should have length 1.\n ")
    }else if(identical(cpus,"all")){
        cpus <- parallel::detectCores()
    }else if(!is.numeric(cpus) || cpus <=0 || cpus %% 1 != 0){
        stop("Argument \'cpus\' should be a positive integer or \'all\'.\n ")
    }else if(cpus>1 && cpus>parallel::detectCores()){
        stop("Only ",parallel::detectCores()," CPU cores are available. \n")
    }

    ## dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** recover dataset
    method.merge <- c("average", "pool.fixse", "pool.se", "pool.gls", "pool.gls1", "pool.rubin")
    keep.cols <- c("by","parameter","estimate","se","df","statistic","lower","upper","null","p.value")
        
    Ucluster <- sort(unique(unlist(data[[cluster]])))
    n.cluster <- length(Ucluster)
    index.cluster <- split(1:NROW(data), data[[cluster]], identity)
    template <- list(multivariate = object$multivariate,
                     univariate = object$univariate)
    template$multivariate[,c("statistic","df.num","df.denom","p.value")] <- NULL
    template$univariate[,c("estimate","se","df","statistic","lower","upper","p.value")] <- NULL
                     
    warperResample <- function(iSample){

        ## *** resample
        iNewCluster <- stats::setNames(sample(n.cluster, n.cluster, replace = TRUE),
                                       1:n.cluster)
        iNewData <- data[]
        if(type == "permutation"){
            ## permute X-values (constant within cluster)
            iPerm <- sample(n.cluster, replace = FALSE)
            iData <- data
            for(iCluster in 1:n.cluster){ ## iCluster <- 1
                iData[index.cluster[[iCluster]],variable.perm] <- Uvar[[iPerm[iCluster]]] ## seems to properly expand Uvar over multiple timepoints
            }
        }else if(type == "bootstrap"){
            ils.Boot <- index.cluster[sample(n.cluster, replace = TRUE)]
            names(ils.Boot) <- 1:n.cluster
            iBoot <- unlist(ils.Boot)
            
            iData <- data[iBoot,,drop=FALSE]
            iData[[cluster]] <- as.numeric(factor(names(iBoot), levels = unique(names(iBoot))))
            rownames(iData) <- NULL            
        }

        ## *** re-estimate
        iCall <- object.call
        iCall$trace <- 0
        iCall$data <- iData
        iEstimate <- eval(iCall)
        
        ## *** export
        ## use merge in case effects for some categories are not estimated
        iEstimate$univariate2 <- model.tables(iEstimate, columns = keep.cols, method = method)

        iOut <- template
        if(is.null(method) || (method %in% method.merge == FALSE)){
            iOut$multivariate <- merge(iOut$multivariate, iEstimate$multivariate, by = c("type","test","null"), all = TRUE)
            iOut$univariate <- merge(iOut$univariate, iEstimate$univariate2, by = c("by","parameter"), all = TRUE)
        }else{
            iOut$univariate <- iEstimate$univariate2
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
    
    if(!is.null(seed)){set.seed(seed)}

    if(cpus>1){
        cl <- parallel::makeCluster(cpus)
        fct2export <- NULL
        parallel::clusterExport(cl = cl, varlist = fct2export, envir = as.environment(asNamespace("LMMstar")))
    }else{
        cl <- NULL
    }

    if(!is.null(cl) || trace){
        ls.sample <- pbapply::pblapply(1:n.sample, warperResample, cl = cl)
    }else{
        ls.sample <- lapply(1:n.sample, warperResample)
    }

    if(cpus>1){
        parallel::stopCluster(cl)
    }

    ## ** post process
    keep.cols2 <- c("estimate", "se", "df", "statistic", "lower", "upper", "null", "p.value")
    
    out <- list()
    sample.multivariate <- do.call(cbind,lapply(ls.sample, function(iSample){iSample$multivariate$statistic}))
    sample.univariate <- do.call(cbind,lapply(ls.sample, function(iSample){iSample$univariate$statistic}))

    if(type == "permutation"){
        obs.multivariate <- matrix(object$multivariate[,"statistic"], nrow = NROW(object$multivariate), ncol = NCOL(sample.multivariate))
        obs.univariate <- matrix(object$univariate[,"statistic"], nrow = NROW(object$univariate), ncol = NCOL(sample.univariate), byrow = FALSE)
        
        if(is.null(method) || (method %in% method.merge == FALSE)){
            out$multivariate <- object$multivariate
            out$multivariate$df.num <- NA
            out$multivariate$df.denom <- NA
            out$multivariate$p.value <- (rowSums(abs(sample.multivariate) > abs(obs.multivariate), na.rm = TRUE)+1)/(rowSums(!is.na(sample.multivariate))+1)
            attr(out$univariate,"sample") <- sample.multivariate
        }

        if(is.null(method) || (method %in% method.merge == FALSE)){
            out$univariate <- object$univariate
        }else{
            out$univariate <- model.tables(object, method = method, columns = keep.cols2)
        }
        out$univariate$se <- NA
        out$univariate$df <- NA
        out$univariate$p.value <- (rowSums(abs(sample.univariate) > abs(obs.univariate), na.rm = TRUE)+1)/(rowSums(!is.na(sample.univariate))+1)
        attr(out$univariate,"sample") <- sample.univariate

    }else if(type == "bootstrap"){
        alpha <- 1-level

        if(is.null(method) || (method %in% method.merge == FALSE)){
            out$univariate <- object$univariate
        }else{
            out$univariate <- model.tables(object, method = method, columns = keep.cols2)
        }
        out$univariate$se <- apply(sample.univariate, MARGIN = 1, FUN = stats::sd, na.rm = TRUE)
        out$univariate$statistic <- out$univariate$estimate/out$univariate$se
        out$univariate$df <- NA
        out$univariate$lower <- apply(sample.univariate, MARGIN = 1, FUN = stats::quantile, probs = alpha/2, na.rm = TRUE)
        out$univariate$upper <- apply(sample.univariate, MARGIN = 1, FUN = stats::quantile, probs = 1-alpha/2, na.rm = TRUE)

        if(!inherits(calc.pvalue,"try-error")){
            out$univariate$p.value <- sapply(name.keepcoef, function(iName){
                BuyseTest::boot2pvalue(stats::na.omit(sample.univariate[,iName]),
                                       null = 0,
                                       estimate = out[iName,"estimate"],
                                       alternative = "two.sided",
                                       add.1 = correction)
            })
        }
        attr(out$univariate,"sample") <- sample.univariate
    }
    
    ## ** export
    return(out)
}

##----------------------------------------------------------------------
### resample.R ends here
