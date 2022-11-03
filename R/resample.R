### resample.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 31 2022 (10:09) 
## Version: 
## Last-Updated: nov  3 2022 (11:27) 
##           By: Brice Ozenne
##     Update #: 164
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Inference via Resampling for Linear Mixed Model
##' @description Non-parametric bootstrap or permutation test for Linear Mixed Models.
##'
##' @param object a \code{lmm} object.
##' @param type [character] should permutation test (\code{"perm"}) or non-parametric bootstrap (\code{"boot"}) be used?
##' @param effects [character vector] the variable(s) to be permuted or the effect(s) to be tested via non-parametric bootstrap.
##' @param n.sample [integer] the number of samples used.
##' @param trace [logical] should the execution of the resampling be traced?
##' @param seed [integer] Random number generator (RNG) state used when starting resampling.
##' @param cl a cluster object passed to \code{pbapply::pblapply}. On plateform other than Windows can also be an integer indicating the number of child-processes for parallel evaluation.
##'
##' @details Both approach are carried at the cluster level: \itemize{
##' \item Bootstrap: sampling with replacement clusters. If a cluster is picked twice then different cluster id is used for each pick.
##' \item Permutation: permuting covariate values between clusters. This only lead to the null hypothesis when the covariate values are constant within clusters.
##' }
##' 
##' @export
resample <- function(object, type, effects, n.sample = 1e3, trace = TRUE, seed = NULL, cl = NULL){

    ## ** check user input
    type <- match.arg(type, c("perm","perm-var","boot"))
    if(type=="perm"){
        type <- "perm-var"
    }

    name.meanvar <- attr(object$design$mean,"variable")
    value.meancoef <- coef(object, effects = "mean")
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
                 "Consider using type=\"perm-res\" or type=\"boot\" instead of type=\"perm-var\". \n")
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
    if(object$opt$name=="gls"){
        stop("Function resample not compatible with \"gls\" optimizer\n.")
    }

    ## ** initialize
    ## *** data
    data <- object$data
    var.cluster <- object$cluster$var
    vec.id <- data[[var.cluster]]
    vec.Uid <- unique(vec.id)
    n.id <- length(vec.Uid)
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
        object0 <- .constrain.lmm(object, effects = stats::setNames(rep(0,length(effects)),effects))
        Xbeta0 <- stats::predict(object0, newdata = data[,var.mean,drop=FALSE], se = FALSE, df = FALSE)[,1]
        epsilon.norm <- stats::residuals(object0, type = "normalized")
        OmegaChol <- stats::sigma(object0, chol = TRUE, cluster = as.character(vec.Uid))
    }

    ## ** function
    warperResample <- function(iSample){

        ## *** resample
        if(type == "perm-var"){
            ## permute X-values
            iPerm <- sample(n.id, replace = FALSE)
            iData <- data
            for(iCluster in 1:n.id){ ## iCluster <- 1
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
            iEpsilon.norm <- epsilon.norm[iPerm]
            ## rescale residuals
            for(iCluster in 1:n.id){ ## iCluster <- 1
                iData[[var.outcome]][index.cluster[[iCluster]]] <- OmegaChol[[iCluster]] %*% iEpsilon.norm[index.cluster[[iCluster]]]
            }
            ## add fixed effects
            iData[[var.outcome]] <- iData[[var.outcome]] + Xbeta0
            ## range(iData[[var.outcome]] - data[iPerm,var.outcome])

        }else if(type == "boot"){
            ils.Boot <- index.cluster[sample(n.id, replace = TRUE)]
            names(ils.Boot) <- 1:n.id
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
                                         stratify.mean = FALSE,
                                         precompute.moments = precompute.moments)
        }else{
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
            }else{
                iPrecompute.XX <- NULL
                iPrecompute.XY <- NULL
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
                                   trace = FALSE))


        ## *** export
        if(inherits(iEstimate,"try-error")){
            return(NULL)
        }else{
            return(iEstimate$estimate[name.keepcoef])
        }

    }
    ## res <- warperResample(1)

    ## ** run
    if(!is.null(seed)){set.seed(seed)}

    if((!is.null(cl) || trace) && requireNamespace("pbapply")){
        ls.sample <- pbapply::pblapply(1:n.sample, warperResample, cl = cl)
    }else{
        ls.sample <- lapply(1:n.sample, warperResample)
    }

    ## ** post-process
    M.sample <- do.call(rbind,ls.sample)
    n.sample <- NROW(M.sample)

    out <- model.tables(object)[name.keepcoef,,drop=FALSE]*NA
    out$estimate <- as.double(value.meancoef[name.keepcoef])

    M.estimate <- matrix(value.meancoef[name.keepcoef], nrow = n.sample, ncol = length(name.keepcoef), byrow = TRUE,
                         dimnames = list(NULL,name.keepcoef))

    if(type %in% c("perm-var","perm-res") ){
        
        out$p.value <- (colSums(abs(M.sample) > abs(M.estimate), na.rm = TRUE)+1)/(colSums(!is.na(M.sample))+1)
        
    }else if(type == "boot"){

        out$se <- apply(M.sample, MARGIN = 2, FUN = stats::sd, na.rm = TRUE)
        out$lower <- apply(M.sample, MARGIN = 2, FUN = stats::quantile, probs = 0.025, na.rm = TRUE)
        out$upper <- apply(M.sample, MARGIN = 2, FUN = stats::quantile, probs = 0.975, na.rm = TRUE)

        M.sampleH0 <- scale(M.sample, center = TRUE, scale = FALSE)
        out$p.value <- (colSums(abs(M.sampleH0) > abs(M.estimate), na.rm = TRUE)+1)/(colSums(!is.na(M.sample))+1)
        ## out$p.value <- unlist(lapply(name.keepcoef, function(iName){ BuyseTest::boot2pvalue(x = M.sample[,iName], null = 0, estimate = value.meancoef[iName])} ))
    }
    attr(out,"M.sample") <- M.sample
    attr(out,"n.sample") <- n.sample

    return(out)

}

##----------------------------------------------------------------------
### resample.R ends here
