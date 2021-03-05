### lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: mar  5 2021 (23:08) 
##           By: Brice Ozenne
##     Update #: 297
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * lmm (documentation)
##' @title Linear mixed model
##' @description Fit a linear mixed model using either a compound symmetry structure or an unstructured covariance matrix.
##' This is essentially an interface to the \code{nlme::gls} function.
##'
##' @param formula [formula] Specify the model for the mean.
##' On the left hand side the outcome and on the right hand side the covariates affecting the mean value.
##' E.g. Y ~ Gender + Gene.
##' @param covariance [formula] Specify the model for the covariance.
##' On the right hand side the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' On the left hand side, a possible stratification variable, e.g. group ~ time|id. In that case the mean structure should only be stratified on this variable using interactions.
##' @param structure [character] type of covariance structure, either \code{"CS"} (compound symmetry) or \code{"UN"} (unstructured).
##' @param data [data.frame] dataset (in the long format) containing the observations.
##' @param df [logical] Should the degree of freedom be computed using a Satterthwaite approximation?
##' @param ... passed to \code{nlme::gls}.

## * lmm (examples)
##' @examples
##' ## simulate data in the wide format
##' library(lava)
##' m <- lvm(c(Y1,Y2,Y3,Y4) ~ age + gender)
##' categorical(m, labels = c("male","female")) <- ~gender
##' transform(m, id~gender) <- function(x){1:NROW(x)}
##' distribution(m, ~age) <- gaussian.lvm(mean = 50, sd = 10)
##'
##' set.seed(10)
##' dW <- lava::sim(m, 1e2)
##'
##' ## move to the long format
##' name.varying <- paste0("Y",1:4)
##' dL <- reshape(dW, direction  = "long",
##'               idvar = c("id","age","gender"),
##'               varying = name.varying,
##'               v.names = "Y",
##'               timevar = "visit")
##' rownames(dL) <- NULL
##' dL$visit <- factor(dL$visit,
##'                    levels = 1:length(name.varying),
##'                    labels = name.varying)
##' head(dL)
##' 
##' ## fit mixed model
##' eCS.lmm <- lmm(Y ~ visit + age + gender, covariance = ~visit|id, structure = "CS", data = dL, debug = 2)
##' eCS.lmm <- lmm(Y ~ visit, covariance = ~visit|id, structure = "CS", data = dL, debug = 2)
##' vcov(eCS.lmm, type = "lmm")
##' vcov(eCS.lmm, type = "gls")
##' 
##' 
##' eUN.lmm <- lmm(Y ~ visit*gender + age* gender, covariance = ~visit|id, structure = "UN", data = dL, debug = 2)
##' 
##' eCSs.lmm <- lmm(Y ~ visit*gender + age* gender, covariance = gender~visit|id, structure = "CS", data = dL, debug = 2)
##' summary(e0.lmm)
##' 
##'
##' e.lmm <- lmm(Y ~ visit + age + gender, covariance = ~visit|id, data = dL)
##' cat(attr(e.lmm,"code")) ## code used to fit the model
##' head(attr(e.lmm,"data")) ## data used to fit the model
##' summary(e.lmm)

## * lmm (code)
##' @export
lmm <- function(formula, covariance, structure, data, df = FALSE, debug = FALSE, ...){

    ## ** check and normalize user imput
    if(debug>=1){cat("1. check and normalize user imput \n")}
    
    ## *** formula
    if(debug>=2){cat("- formula")}
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
    var.outcome <- name.mean[1]
    var.X <- name.mean[-1]
    formula.terms <- terms(formula)
    if(any(attr(formula.terms,"order")>2)){
        stop("Cannot handle interaction involving more than two variables. \n")
    }

    if(debug>=2){cat("\n")}
    ## *** covariance 
    if(debug>=2){cat("- covariance ")}
    if(!inherits(covariance,"formula")){
        stop("Argument \'covariance\' must be of class formula, something like: ~ time|id or group ~ time|id. \n")
    }

    name.vcov <- all.vars(covariance)
    if(any(name.vcov %in% names(data) == FALSE)){
        invalid <- name.vcov[name.vcov %in% names(data) == FALSE]
        stop("Argument \'formula\' is inconsistent with argument \'data\'. \n",
             "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
             sep = "")
    }

    ## - right hand side
    if(debug>=2){cat(" (rhs ")}
    if(length(rhs.vars(covariance))!=2){
        stop("Incorrect specification of argument \'covariance\'. \n",
             "There should be exactly two variables on the right hand side, something like: ~ time|id or group ~ time|id. \n")
    }

    if(!grepl("|",deparse(covariance),fixed = TRUE)){
        stop("Incorrect specification of argument \'covariance\'. \n",
             "No | symbol found so no grouping variable could be defined. \n",
             "Shoud be something like: ~ time|id or group ~ time|id. \n")
    }

    if(length(grepl("|",deparse(covariance),fixed = TRUE))>1){
        stop("Incorrect specification of argument \'covariance\'. \n",
             "The symbol | should only appear once, something like: ~ time|id or group ~ time|id. \n")
    }
    res.split <- strsplit(deparse(covariance),"|", fixed = TRUE)[[1]]
    var.cluster <- trimws(res.split[2], which = "both")
    var.time <- all.vars(update(stats::as.formula(res.split[1]),0~.))
    if(length(var.time)!=1){
        stop("Incorrect specification of argument \'covariance\'. \n",
             "Should have exactly one variable before the grouping symbol (|), something like: ~ time|id or group ~ time|id. \n")
    }
    if(length(var.cluster)!=1){
        stop("Incorrect specification of argument \'covariance\'. \n",
             "Should have exactly one variable after the grouping symbol (|), something like: ~ time|id or group ~ time|id. \n")
    }
    test.duplicated <- tapply(data[[var.time]], data[[var.cluster]], function(iT){any(duplicated(iT))})
    if(any(test.duplicated)){
        stop("Incorrect time variable, it should contain unique values within clusters \n")
    }

    U.cluster <- sort(unique(data[[var.cluster]]))
    n.cluster <- length(U.cluster)
    U.time <- sort(unique(data[[var.time]]))
    n.time <- length(U.time)
    n.obspercluster <- tapply(data[[var.time]],data[[var.cluster]], length)

    

    ## *** left hand side
    if(debug>=2){cat(" lhs) ")}
    if(length(lhs.vars(covariance))==0){
        var.strata <- NULL
    }else{
        if(length(lhs.vars(covariance))!=1){
            stop("Incorrect specification of argument \'covariance\'. \n",
                 "There should be at most one variable on the left hand side, something like: ~ time|id or group ~ time|id. \n")
        }else{
            var.strata <- name.vcov[1]
            U.strata <- sort(unique(data[[var.strata]]))
            n.strata <- length(U.strata)
            
            tocheck <- setdiff(var.X,var.strata)

            if(var.strata %in% var.X == FALSE){
                stop("When a variable is used to stratify the covariance structure, it should also be use to stratify the mean structure. \n",
                     "Consider adding all interactions with \"",var.strata,"\" in the argument \'formula\'. \n")
            }
            
            if(any(rowSums(table(data[[var.cluster]],data[[var.strata]])>0)!=1)){
                stop("When a variable is used to stratify the covariance structure, all observations belonging to each cluster must belong to a single strata. \n")
            }

            sapply(tocheck, function(iX){ ## iX <- "age"
                iCoef <- which(attr(formula.terms,"factors")[iX,]==1)
                iInteraction <- attr(formula.terms,"factors")[var.strata,iCoef,drop=FALSE]
                if(length(unique(iInteraction))!=n.strata){
                    stop("When a variable is used to stratify the covariance structure, it should also be used to stratify the mean structure. \n",
                         "Consider adding an interaction between \"",iX,"\" and \"",var.strata,"\" in the argument \'formula\'. \n")
                }
            })

            if(n.strata>max(5,n.cluster/10)){
                warning("Large number of strata - ",n.strata," - may lead to optimization issue. \n")
            }


        }
    }


    if(debug>=2){cat("\n")}

    ## *** data
    if(debug>=2){cat("- data")}
    data <- as.data.frame(data)
    if("XXindexXX" %in% names(data)){
        stop("Incorrect specification of argument \'data\'. \n",
             "The variable \"XXindexXX\" is used internally but already exists in \'data\' \n")
    }
    data$XXindexXX <- 1:NROW(data)
    
    if(is.null(var.strata)){
        var.strata <- "XXstrata.indexXX"
        U.strata <- 1
        n.strata <- 1
        if("XXstrata.indexXX" %in% names(data)){
            stop("Incorrect specification of argument \'data\'. \n",
                 "The variable \"XXstrata.indexXX\" is used internally but already exists in \'data\' \n")
        }
        data$XXstrata.indexXX <- 1
    }else{
        data[[var.strata]] <- as.character(data[[var.strata]])
    }
    
    var.time.index <- paste0("XX",var.time,".indexXX")
    if(var.time.index %in% names(data)){
        stop("Incorrect specification of argument \'data\'. \n",
             "The variable ",var.time.index," is used internally but already exists in \'data\' \n")
    }
    data[[var.time.index]] <- as.numeric(as.factor(data[[var.time]]))

    if(debug>=2){cat("\n")}

    ## *** structure
    if(debug>=2){cat("- structure")}
    structure <- match.arg(toupper(structure), c("CS","UN"))
    if(n.strata==1){
        txt.data <- "data"
    }else{
        txt.data <- paste0("data[data$",var.strata,"==\"",U.strata,"\",,drop=FALSE]")
    }

    
    if(n.strata>1 && var.strata %in% var.X){
        term.exclude <- c(var.strata,paste0(setdiff(var.X, var.strata),":",var.strata))
        formula.mean <- update(formula, paste0(".~.-",paste0(term.exclude,collapse="-")))
    }else{
        formula.mean <- formula
    }


    if(debug>=2){cat("\n")}

    ## ** fit mixed model
    if(debug>=1){cat("2. fit mixed model")}

    if(all(n.obspercluster==1)){
        if(structure == "CS"){
            txt.gls <- paste0("nlme::gls(",deparse(formula.mean),",
                                         data = ",txt.data,", ..., )")
        }else if(structure == "UN"){
            form.var <- stats::as.formula(paste0("~1|",var.time))
            txt.gls <- paste0("nlme::gls(",deparse(formula.mean),",
                                         weights = nlme::varIdent(form = ",deparse(form.var),"),
                                         data = ",txt.data,", ...)")
        }
    }else{
        if(structure == "CS"){
            form.cor <- stats::as.formula(paste0("~1|",var.cluster))
            txt.gls <- paste0("nlme::gls(",deparse(formula.mean),",
                                         correlation = nlme::corCompSymm(form = ",deparse(form.cor),"),
                                         data = ",txt.data,", ...)")
        }else if(structure == "UN"){
            form.var <- stats::as.formula(paste0("~1|",var.time))
            form.cor <- stats::as.formula(paste0("~",var.time.index,"|",var.cluster))
            txt.gls <- paste0("nlme::gls(",deparse(formula.mean),",
                                          correlation = nlme::corSymm(form = ",deparse(form.cor),"),
                                          weights = nlme::varIdent(form = ",deparse(form.var),"),
                                          data = ",txt.data,", ...)")
        }
    }
    e.gls <- setNames(lapply(txt.gls, function(iTxt){eval(parse(text = iTxt))}),
                      U.strata)

    if(debug>=1){cat("\n")}

    ## ** convert to lmm object
    if(debug>=1){cat("3. convert to lmm object \n")}

    ## *** extract information from gls
    if(debug>=2){cat("- extract information from gls")}
    callGLS <- lapply(e.gls, function(iM){
        paste0(gsub(",",",\n    ",gsub(" ","",paste(deparse(iM$call), collapse = ""))),"\n")
    })

    if(n.strata==1){
        info <- .getInfoGLS(object = e.gls[[1]], data = data, strata = NULL,
                            name.time = U.time, name.cluster = U.cluster, name.outcome = var.outcome)
        M.beta <- rbind(names(info$coef$beta))
        M.sigma <- rbind(names(info$coef$sigma))
        M.k <- rbind(names(info$coef$k))
        M.rho <- rbind(names(info$coef$rho))

        X <- info$design
        info.cluster <- info$cluster
        info.time <- info$time

        info.coef <- info$coef$all
        info.type <- info$coef$type
        info.df <- info$df
        info.betavcov <- info$betavcov
        info.score <- info$core
        info.logLik <- info$logLik
        
        info.resvcov <- info$resvcov
    }else{
        browser()
        info <- lapply(U.strata, function(iStrata){
            .getInfoGLS(object = e.gls[[iStrata]], data = data, strata = iStrata,
                        name.time = U.time, name.outcome = var.outcome)
        })
    }
    
    if(debug>=2){cat("\n")}
    
    ## *** merge into lmm
    if(debug>=2){cat("- merge into lmm")}

    e.lmm <- list(call = match.call(),
                  gls = e.gls,
                  data = data,
                  mean.structure = list(formula = formula, coef.beta = M.beta),
                  variance.structure = list(formula = covariance, type = structure, coef.sigma = M.sigma, coef.k = M.k, coef.rho = M.rho),
                  X = X,
                  cluster = info.cluster,
                  time = info.time,
                  variable = list(outcome = var.outcome,
                                  covariates = var.X,
                                  strata = if(n.strata>1){var.strata}else{NULL},
                                  time = var.time,
                                  cluster = var.cluster),
                  coef = info.coef,
                  coeftype = info.type,
                  df = info.df,
                  betavcov = info.betavcov,
                  logLik = info.logLik,
                  score = info.score,
                  resvcov = info.resvcov,
                  code = callGLS
                  )
                
    
    if(debug>=2){cat("\n")}

    ## ** export
    class(e.lmm) <- append("lmm",class(e.lmm))
    return(e.lmm)
}

## * .getInfoGLS
.getInfoGLS <- function(object, data, strata, name.time, name.cluster, name.outcome){
    data <- eval(object$call$data)

    ## ** extract correlation/variance structure according to data
    if(!is.null(object$modelStruct$corStruct)){
        corgroups <- object$groups
    }
    if(!is.null(object$modelStruct$varStruct)){
        timegroups <- getGroups(object$modelStruct$varStruct)
    }else if(!is.null(object$modelStruct$corStruct)){
        timegroups <- corgroups
        for(iCluster in unique(corgroups)){
            timegroups[corgroups == iCluster] <- attr(object$modelStruct$corStruct, "covariate")[[iCluster]]
        }
        timegroups <- droplevels(timegroups)
    }else{
        timegroups <- rep(name.outcome, NROW(data))
    }

    ## ** number and position of the clusters among the observations
    if(is.null(object$modelStruct$corStruct)){ ## no correlation structure
        index.cluster <- 1:NROW(data)
        n.cluster <- NROW(data)

    }else{ ## correlation structure
        index.cluster <- setNames(as.numeric(corgroups), corgroups)
        n.cluster <- attr(object$modelStruct$corStruct,"Dim")$M
    }

    ## ** extract number and position of unique residual variance-covariance structures
    if(is.null(object$modelStruct$varStruct)){ ## no variance structure
        index.vargroup <- setNames(rep(1,n.cluster), sort(unique(names(index.cluster))))
        n.vargroup <-  1
        Sigma.pattern <- list(rep("sigma",length(name.time)))
        
    }else if(is.null(object$modelStruct$corStruct)){ ## no correlation structure
        index.vargroup <- as.numeric(as.factor(timegroups))        
        n.vargroup <-  max(index.vargroup)
        Sigma.pattern <- as.list(levels(as.factor(timegroups)))
        
    }else{ ## variance and correlation structure
        ## variance parameter within cluster
        variance.per.cluster <- tapply(X = timegroups, INDEX = corgroups, FUN = function(iVec){list(iVec)}) ## WARNING MAY MESS UP THE ORDER
        variance.per.cluster <- variance.per.cluster[levels(corgroups)] ## PROPERLY REORDER
            
        ## unique variance patterns
        Sigma.pattern <- unique(variance.per.cluster)
        names(Sigma.pattern) <- sapply(Sigma.pattern, paste, collapse = "|")
        n.vargroup <- length(Sigma.pattern)
    
        ## associate each cluster to a variance structure
        index.vargroup <- sapply(variance.per.cluster, function(x){
            as.double(which(unlist(lapply(Sigma.pattern, identical, x))))
        })

    }
    
    Sigma.pattern.full <- Sigma.pattern[sapply(Sigma.pattern, length) == length(name.time)]
    names(Sigma.pattern.full) <- sapply(Sigma.pattern.full, paste, collapse = "|")

    ## ** coefficients
    theta <- getCoef(object, effects = c("mean","variance"), add.type = TRUE)
    beta <- setNames(theta[theta$type=="mean","estimate"],rownames(theta)[theta$type=="mean"])
    sigma <- setNames(theta[theta$type=="std.residual","estimate"],rownames(theta)[theta$type=="std.residual"])
    k <- setNames(theta[theta$type=="factor.std.residual","estimate"],rownames(theta)[theta$type=="factor.std.residual"])
    rho <- setNames(theta[theta$type=="correlation","estimate"],rownames(theta)[theta$type=="correlation"])
    
    ## ** residual variance-covariance structure
    Omega <- .getVarCov(object,
                        sigma = sigma,
                        k = k,
                        rho = rho,
                        name.time = name.time,
                        Sigma.pattern = Sigma.pattern,
                        Sigma.pattern.full = Sigma.pattern.full
                        )
    precision.vargroup <- lapply(Omega$variance.vargroup, solve)

    ## ** Likelihood and derivatives
    Y <- data[[name.outcome]]
    X <- model.matrix(stats::formula(object), data) 
    logLik <- .logLik(Y = Y, X = X, beta = beta, sigma = sigma, k = k, rho = rho, precision = precision.vargroup,
                      index.variance = index.vargroup, index.cluster = index.cluster, indiv = FALSE, REML = object$method=="REML")
    score <- .score(Y = Y, X = X, beta = beta, sigma = sigma, k = k, rho = rho, precision = precision.vargroup,
                    index.variance = index.vargroup, index.cluster = index.cluster, indiv = TRUE, REML = object$method=="REML")
    hessian <- .hessian(X = X, beta = beta, sigma = sigma, k = k, rho = rho, precision = precision.vargroup,
                        index.variance = index.vargroup, index.cluster = index.cluster, indiv = TRUE, REML = object$method=="REML")
    information <- -apply(hessian, FUN = sum, MARGIN = 2:3)

    ## ** merge
    if(!is.null(strata)){
        names(sigma) <- paste0("sigma:",strata)
        names(k) <- paste0(names(k),strata)
        names(rho) <- paste0(names(rho),strata)
        names(beta) <- paste0(names(beta),strata)
    }else{
        names(sigma) <- "sigma"
    }
    out <- list(coef = list(beta = beta, sigma = sigma, k = k, rho = rho,
                            all = c(beta,sigma,k,rho),
                            type = c(rep("mean",length(beta)),"sigma",rep("k",length(k)),rep("rho",length(rho)))),
                design = X,
                resvcov = list(n = length(Omega$variance.vargroup), full = Omega$variance.full, all = Omega$variance.vargroup, all.inverse = precision.vargroup, index = index.vargroup),
                cluster = list(n = n.cluster, name = name.cluster, index = data[index.cluster,"XXindexXX"]),
                time = list(n = length(name.time), name = name.time, index = timegroups),
                logLik = logLik,
                score = score,
                hessian = hessian,
                information = information,
                betavcov = solve(information),
                df = object$dim$N - object$dim$p * (object$method=="REML")
                )

    ## ** check
    if(object$method=="ML"){
        test.logLik <- abs(out$logLik - as.double(logLik(object)))
        if(test.logLik>1e-10){
            warning("Mismatch between the gls and lmm log likelihood (largest difference ",test.logLik,"). \n",
                    "Consider contacting the package manager. \n")
        }

        test.vcov <- max(abs(vcov(object) * (object$dim$N-object$dim$p) / object$dim$N - out$betavcov))
        if(test.vcov>1e-10){
            warning("Mismatch between the gls and lmm variance covariance matrix (largest difference ",test.vcov,"). \n",
                    "Consider contacting the package manager. \n")
        }

        
        
    }else{
        test.logLik <- abs(out$logLik - as.double(logLik(object)))
        if(test.logLik>1e-10){
            warning("Mismatch between the gls and lmm log likelihood (largest difference ",test.logLik,"). \n",
                    "Consider contacting the package manager. \n")
        }
        
        logLikML - logLikREML
        test.vcov <- max(abs(vcov(object) - out$betavcov))
        if(test.vcov>1e-10){
            warning("Mismatch between the gls and lmm variance covariance matrix (largest difference ",test.vcov,"). \n",
                    "Consider contacting the package manager. \n")
        }
    }

    ## ** export
    return(out)
}


######################################################################
### lmm.R ends here
