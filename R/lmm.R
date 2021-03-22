### lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: mar 22 2021 (23:27) 
##           By: Brice Ozenne
##     Update #: 460
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
##' @param variance [formula] Specify the model for the covariance.
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
##' eCS.lmm <- lmm(Y ~ visit + age + gender, variance = ~visit|id, structure = "CS", data = dL, debug = 2)
##' eCS.lmm <- lmm(Y ~ visit, variance = ~visit|id, structure = "CS", data = dL, debug = 2)
##' vcov(eCS.lmm, type = "lmm")
##' vcov(eCS.lmm, type = "gls")
##' 
##' 
##' eUN.lmm <- lmm(Y ~ visit*gender + age* gender, variance = ~visit|id, structure = "UN", data = dL, debug = 2)
##' 
##' eCSs.lmm <- lmm(Y ~ visit*gender + age* gender, variance = gender~visit|id, structure = "CS", data = dL, debug = 2)
##' summary(e0.lmm)
##' 
##'
##' e.lmm <- lmm(Y ~ visit + age + gender, variance = ~visit|id, data = dL)
##' cat(attr(e.lmm,"code")) ## code used to fit the model
##' head(attr(e.lmm,"data")) ## data used to fit the model
##' summary(e.lmm)

## * lmm (code)
##' @export
lmm <- function(formula, variance, structure, data, df = FALSE, debug = FALSE, ...){
    out <- list()
    
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
    
    ## *** variance 
    if(debug>=2){cat("- variance ")}
    if(!inherits(variance,"formula")){
        stop("Argument \'variance\' must be of class formula, something like: ~ time|id or group ~ time|id. \n")
    }

    name.vcov <- all.vars(variance)
    if(any(name.vcov %in% names(data) == FALSE)){
        invalid <- name.vcov[name.vcov %in% names(data) == FALSE]
        stop("Argument \'formula\' is inconsistent with argument \'data\'. \n",
             "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
             sep = "")
    }

    ## - right hand side
    if(debug>=2){cat(" (rhs ")}

    if(length(all.vars(update(variance,0~.)))==0){ ## t-test

        if("XXclusterXX" %in% names(data)){
            stop("Incorrect specification of argument \'data\'. \n",
                 "The variable \"XXclusterXX\" is used internally but already exists in \'data\' \n")
        }
        data$XXclusterXX <- 1:NROW(data)
        if("XXtimeXX" %in% names(data)){
            stop("Incorrect specification of argument \'data\'. \n",
                 "The variable \"XXtimeXX\" is used internally but already exists in \'data\' \n")
        }
        data$XXtimeXX <- "1"
        var.cluster <- "XXclusterXX"
        formula.var <- ~XXtimeXX
        var.time <- "XXtimeXX"

    }else{

        if(!grepl("|",deparse(variance),fixed = TRUE)){
            stop("Incorrect specification of argument \'variance\'. \n",
                 "No | symbol found so no grouping variable could be defined. \n",
                 "Shoud be something like: ~ time|id or group ~ time|id. \n")
        }

        if(length(grepl("|",deparse(variance),fixed = TRUE))>1){
            stop("Incorrect specification of argument \'variance\'. \n",
                 "The symbol | should only appear once, something like: ~ time|id or group ~ time|id. \n")
        }
        res.split <- strsplit(deparse(variance),"|", fixed = TRUE)[[1]]
        var.cluster <- trimws(res.split[2], which = "both")
        formula.var <- stats::as.formula(res.split[1])
        var.time <- all.vars(update(formula.var,0~.))[1]

        if(length(var.time)==0){
            stop("Incorrect specification of argument \'variance\'. \n",
                 "There should be at least one variable before the grouping symbol (|), something like: ~ time|id or group ~ time|id. \n")
        }
        if(length(var.cluster)!=1){
            stop("Incorrect specification of argument \'variance\'. \n",
                 "Should have exactly one variable after the grouping symbol (|), something like: ~ time|id or group ~ time|id. \n")
        }

        test.duplicated <- tapply(data[[var.time]], data[[var.cluster]], function(iT){any(duplicated(iT))})
        if(any(test.duplicated)){
            stop("Incorrect specification of argument \'variance\'. \n",
                 "The time variable (first variable before |) should contain unique values within clusters \n")
        }
    }

    ## *** left hand side
    if(debug>=2){cat(" lhs) ")}
    if(length(lhs.vars(variance))==0){
        var.strata <- NULL
    }else{
        if(length(lhs.vars(variance))!=1){
            stop("Incorrect specification of argument \'variance\'. \n",
                 "There should be at most one variable on the left hand side, something like: ~ time|id or group ~ time|id. \n")
        }else{
            var.strata <- name.vcov[1]
            U.strata <- as.character(sort(unique(data[[var.strata]])))
            n.strata <- length(U.strata)
            
            tocheck <- setdiff(var.X,var.strata)

            if(var.strata %in% var.X == FALSE){
                stop("When a variable is used to stratify the variance structure, it should also be use to stratify the mean structure. \n",
                     "Consider adding all interactions with \"",var.strata,"\" in the argument \'formula\'. \n")
            }
            
            if(any(rowSums(table(data[[var.cluster]],data[[var.strata]])>0)!=1)){
                stop("When a variable is used to stratify the variance structure, all observations belonging to each cluster must belong to a single strata. \n")
            }

            sapply(tocheck, function(iX){ ## iX <- "age"
                iCoef <- which(attr(formula.terms,"factors")[iX,]==1)
                iInteraction <- attr(formula.terms,"factors")[var.strata,iCoef,drop=FALSE]
                if(length(unique(iInteraction))!=n.strata){
                    stop("When a variable is used to stratify the variance structure, it should also be used to stratify the mean structure. \n",
                         "Consider adding an interaction between \"",iX,"\" and \"",var.strata,"\" in the argument \'formula\'. \n")
                }
            })

            test.unique <- tapply(data[[var.time]], data[[var.strata]], function(iT){list(sort(unique(iT)))})
            
            if(length(unique(sapply(test.unique,length)))>1){
                stop("Incorrect specification of argument \'variance\'. \n",
                    "The time variable should contain the same number of unique values in each strata \n")
            }
            test.unique <- do.call(rbind,test.unique)
            
            if(any(apply(test.unique, MARGIN = 2, FUN = function(x){length(unique(x))})!=1)){
                stop("Incorrect specification of argument \'variance\'. \n",
                     "The time variable should contain the same unique values in each strata \n")
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
    U.time <- as.character(sort(unique(data[[var.time]])))

    out$strata <- list(n = n.strata, levels = U.strata, var = var.strata)
    out$time <- list(n = length(U.time), levels = U.time, var = var.time)
    out$cluster <- list(var = var.cluster)
    out$outcome <- list(var = var.outcome)

    if(debug>=2){cat("\n")}

    ## *** structure
    if(debug>=2){cat("- structure")}
    structure <- match.arg(toupper(structure), c("CS","UN"))
    if(n.strata==1){
        txt.data <- "data"
    }else{
        txt.data <- paste0("data[data$",var.strata,"==\"",U.strata,"\",,drop=FALSE]")
    }

    if(n.strata>1){
        term.exclude <- c(var.strata,paste0(setdiff(var.X, var.strata),":",var.strata))
        formula.design <- update(formula, paste0(".~.-",paste0(term.exclude,collapse="-"))) ## no interaction with the strata variable
        out$formula <- list(mean = formula, ## formula will contain all interactions with strata (cf check)
                            mean.design = formula.design,
                            var = variance,
                            var.design = formula.var)
    }else{
        formula.design <- formula
        out$formula <- list(mean = formula,
                            mean.design = formula.design,
                            var = variance,
                            var.design = formula.var)
    }

    if(debug>=2){cat("\n")}

    ## ** design matrices
    if(debug>=1){cat("- extract design matrices")}
    out$design <- .model.matrix.lmm(formula.mean = out$formula$mean.design,
                                    formula.var = out$formula$var.design,
                                    data = data, var.outcome = var.outcome,
                                    var.strata = var.strata, U.strata = U.strata,
                                    var.time = var.time, U.time = U.time,
                                    var.cluster = var.cluster,
                                    structure = structure
                                    )

    if(debug>=2){cat("\n")}
    
    ## ** Estimate model parameters
    if(debug>=1){cat("2. estimate model parameters")}

    if(max(out$design$cluster$nobs)==1){
        if(structure == "CS"){
            txt.gls <- paste0("nlme::gls(",deparse(formula.design),",
                                         data = ",txt.data,", ..., )")
        }else if(structure == "UN"){
            form.var <- stats::as.formula(paste0("~1|",var.time))
            txt.gls <- paste0("nlme::gls(",deparse(formula.design),",
                                         weights = nlme::varIdent(form = ",deparse(form.var),"),
                                         data = ",txt.data,", ...)")
        }
    }else{
        if(structure == "CS"){
            form.cor <- stats::as.formula(paste0("~1|",var.cluster))
            txt.gls <- paste0("nlme::gls(",deparse(formula.design),",
                                         correlation = nlme::corCompSymm(form = ",deparse(form.cor),"),
                                         data = ",txt.data,", ...)")
        }else if(structure == "UN"){
            form.var <- stats::as.formula(paste0("~1|",var.time))
            form.cor <- stats::as.formula(paste0("~",var.time.index,"|",var.cluster))
            txt.gls <- paste0("nlme::gls(",deparse(formula.design),",
                                          correlation = nlme::corSymm(form = ",deparse(form.cor),"),
                                          weights = nlme::varIdent(form = ",deparse(form.var),"),
                                          data = ",txt.data,", ...)")
        }
    }
    out$gls <- setNames(lapply(txt.gls, function(iTxt){eval(parse(text = iTxt))}),
                        U.strata)
    out$gls.call <- lapply(out$gls, function(iM){
        paste0(gsub(",",",\n    ",gsub(" ","",paste(deparse(iM$call), collapse = ""))),"\n")
    })
    out$data <- data
    out$method.fit <- unique(sapply(out$gls, "[[", "method"))
    out$structure <- structure
    
    if(n.strata==1){
        out$param <- list(mu = coef(out$gls[[1]])[out$design$param$mu],
                          sigma = setNames(sigma(out$gls[[1]]),out$design$param$sigma),
                          k = setNames(coef(out$gls[[1]]$modelStruct$varStruct, unconstrained = FALSE),out$design$param$k),
                          cor = setNames(coef(out$gls[[1]]$modelStruct$corStruct, unconstrained = FALSE),out$design$param$rho)
                          )
        out$param$strata <- rep(U.strata, length(unlist(out$param)))
    }else{
        param.mu <- lapply(U.strata, function(iS){ coef(out$gls[[iS]]) })
        param.sigma <- lapply(U.strata, function(iS){ sigma(out$gls[[iS]]) })
        param.k <- lapply(U.strata, function(iS){  coef(out$gls[[iS]]$modelStruct$varStruct, unconstrained = FALSE) })
        param.cor <- lapply(U.strata, function(iS){ coef(out$gls[[iS]]$modelStruct$corStruct, unconstrained = FALSE) })
        out$param <- list(mu = setNames(unlist(param.mu),out$design$param$mu),
                          sigma = setNames(unlist(param.sigma), out$design$param$sigma),
                          k = setNames(unlist(param.k), out$design$param$k),
                          cor = setNames(unlist(param.cor), out$design$param$rho)
                          )

        out$param$strata <- c(unlist(lapply(1:n.strata, function(iS){rep(U.strata[iS], times = length(param.mu[[iS]]))})),
                              unlist(lapply(1:n.strata, function(iS){rep(U.strata[iS], times = length(param.sigma[[iS]]))})),
                              unlist(lapply(1:n.strata, function(iS){rep(U.strata[iS], times = length(param.k[[iS]]))})),
                              unlist(lapply(1:n.strata, function(iS){rep(U.strata[iS], times = length(param.cor[[iS]]))})))
    }
    out$param$type <- c(rep("mu", length(out$param$mu)), rep("sigma", length(out$param$sigma)), rep("k", length(out$param$k)), rep("cor", length(out$param$cor)))

    param.allnames <- unlist(lapply(out$param[c("mu","sigma","k","cor")], names))
    names(out$param$strata) <- param.allnames
    names(out$param$type) <- param.allnames
    if(debug>=1){cat("\n")}

    ## ** Compute likelihood derivatives and other useful quantities
    if(debug>=1){cat("3. compute likelihood derivatives and other useful quantities \n")}
    out$residuals <- out$design$Y - out$design$X.mean %*% out$param$mu
    out$Omega <- attr(out$design$X.var,"FUN.Omega")(object = out$design$X.var, sigma = out$param$sigma, k = out$param$k, rho = out$param$cor)
    out$OmegaM1 <- lapply(out$Omega,solve)
    out$logLik <- .logLik(X = out$design$X.mean, residuals = out$residuals, precision = out$OmegaM1,
                          index.variance = out$design$index.vargroup, index.cluster = out$design$index.cluster, indiv = FALSE, REML = out$method.fit=="REML")
    out$score <- .score(X = out$design$X.mean, residuals = out$residuals, precision = out$OmegaM1,
                        index.variance = out$design$index.vargroup, index.cluster = out$design$index.cluster, indiv = TRUE, REML = out$method.fit=="REML")
    out$hessian <- .hessian(X = out$design$X.mean, precision = out$OmegaM1,
                            index.variance = out$design$index.vargroup, index.cluster = out$design$index.cluster, indiv = TRUE, REML = object$method.fit=="REML")
    out$information <- -apply(out$hessian, FUN = sum, MARGIN = 2:3)
    out$vcov <- solve(out$information)
    out$df <- NULL

    ## coef(out$gls[[1]])
    ## getVarCov(out$gls[[1]])
    ## head(out$design$X.mean)
    ## intervals(out$gls[[2]])$cor[,"est."]
    ## .logLik(Y = out$Y, X = out$design$X.mean, beta = out$param$mu, sigma = NULL, k = NULL, rho = NULL, precision = out$OmegaM1,
    ##         index.variance = out$design$index.vargroup[out$design$index.vargroup==1], index.cluster = as.numeric(as.factor(out$design$index.cluster[out$design$index.vargroup==1])), indiv = FALSE, REML = out$method.fit=="REML")

    ## out <- list(coef = list(beta = beta, sigma = sigma, k = k, rho = rho,
    ##                         all = c(beta,sigma,k,rho),
    ##                         type = c(rep("mean",length(beta)),"sigma",rep("k",length(k)),rep("rho",length(rho)))),
    ##             design = X,
    ##             resvcov = list(n = length(Omega$variance.vargroup), full = Omega$variance.full, all = Omega$variance.vargroup, all.inverse = precision.vargroup, index = index.vargroup),
    ##             cluster = list(n = n.cluster, name = name.cluster, index = data[index.cluster,"XXindexXX"]),
    ##             time = list(n = length(name.time), name = name.time, index = timegroups),
    ##             logLik = logLik,
    ##             score = score,
    ##             hessian = hessian,
    ##             information = information,
    ##             betavcov = solve(information),
    ##             df = object$dim$N - object$dim$p * (object$method=="REML")
    ##             )

    ## ** Sanity checks
    if(debug>=1){cat("4. sanity check vs. gls \n")}

    test.logLik <- abs(out$logLik - sum(sapply(out$gls, logLik)))
    if(test.logLik>1e-10){
        warning("Mismatch between the gls and lmm log likelihood (largest difference ",test.logLik,"). \n",
                "Consider contacting the package manager. \n")
    }

    lmm.vcov.mu <- out$vcov[names(out$param$mu),names(out$param$mu),drop=FALSE]
    ## lmm.vcov.Omega <- out$vcov[names(out$param$mu),names(out$param$mu),drop=FALSE]
    if(out$method=="ML"){
        if(n.strata==1){
            GS.vcov.mu <- vcov(out$gls[[1]]) * (out$gls[[1]]$dim$N-out$gls[[1]]$dim$p) / out$gls[[1]]$dim$N
        }else{
            browser()
        }
    }else{
        ## sqrt(out$gls[[1]]$apVar["lSigma","lSigma"])*1.96
        ## diff(log(intervals(out$gls[[1]])$sigma))
        GS.vcov.mu <- as.matrix(Matrix::bdiag(lapply(out$gls,vcov)))
    }
    GS.vcov.Omega <- as.matrix(Matrix::bdiag(lapply(out$gls,function(iM){iM$apVar})))
    if(max(abs(lmm.vcov.mu - GS.vcov.mu))>1e-10){
        warning("Mismatch between the gls and lmm variance covariance matrix (mean, largest difference ",max(abs(lmm.vcov.mu - GS.vcov.mu)),"). \n",
                "Consider contacting the package manager. \n")
    }
    ## if(max(abs(lmm.vcov.Omega - GS.vcov.Omega))>1e-10){
    ##     warning("Mismatch between the gls and lmm variance covariance matrix (variance, largest difference ",max(abs(lmm.vcov.Omega - GS.vcov.Omega)),"). \n",
    ##             "Consider contacting the package manager. \n")
    ## }
     
    
    ## ** convert to lmm and export
    class(out) <- "lmm"
    return(out)
}



######################################################################
### lmm.R ends here
