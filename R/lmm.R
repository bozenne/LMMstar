### lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: Jun  7 2021 (12:27) 
##           By: Brice Ozenne
##     Update #: 833
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * lmm (documentation)
##' @title Fit Multivariate Gaussian Model
##' @description Fit a multivariate gaussian model using either a compound symmetry structure or an unstructured covariance matrix.
##' This is essentially an interface to the \code{nlme::gls} function.
##'
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
##' @param trace [interger, >0] Show the progress of the execution of the function.
##' @param ... passed to \code{nlme::gls}.
##'
##' @details \bold{Computation time} the \code{lmm} has not been developped to be a fast function as, by default, it uses REML estimation with the observed information matrix and uses a Satterthwaite approximation to compute degrees of freedom (this require to compute the third derivative of the log-likelihood which is done by numerical differentiation). The computation time can be substantially reduced by using ML estimation with the expected information matrix and no calculation of degrees of freedom: arguments \code{method.fit="ML"}, \code{type.information="expected"}, \code{df=FALSE}. This will, however, lead to less accurate p-values and confidence intervals in small samples.
##' \cr


## * lmm (examples)
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit Multivariate Gaussian Model
##' eCS.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "CS", data = dL)
##'
##' ## output
##' eCS.lmm
##' summary(eCS.lmm)
##' summary(eCS.lmm, ci = TRUE)
## * lmm (code)
##' @export
lmm <- function(formula, repetition, structure, data, method.fit = NULL, df = NULL, type.information = NULL, trace = NULL, ...){
    out <- list(call = match.call())
    options <- LMMstar.options()
    data <- as.data.frame(data)
    
    ## ** check and normalize user imput
    if(is.null(trace)){
        trace <- options$trace
    }
    if(trace>=1){cat("1. Check and normalize user imput \n")}

    ## *** objective function
    if(is.null(method.fit)){
        method.fit <- options$method.fit
    }else{
        method.fit <- match.arg(method.fit, choices = c("ML","REML"))
    }

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
    
    ## *** formula
    if(trace>=2){cat("- formula")}
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
    if(any(grepl("Intercept",var.X))){
        stop("Argument \'formula\' should not contain a variable called \"Intercept\". \n")
    }
    if(any(grepl(":",var.X,fixed=TRUE))){
        stop("Argument \'formula\' should not contain a variable whose name contain \":\". \n")
    }
    formula.terms <- stats::terms(formula)
    if(any(attr(formula.terms,"order")>2)){
        stop("Cannot handle interaction involving more than two variables. \n")
    }

    if(trace>=2){cat("\n")}
    
    ## *** repetition 
    if(trace>=2){cat("- repetition ")}
    if(missing(repetition)){
        if(inherits(structure,"structure")){
            if(grepl("|",deparse(structure$formula), fixed = TRUE)){
                repetition <- structure$formula
            }else{
                if(is.null(structure$formula)){
                    if("XXtimeXX" %in% names(data)){
                        stop("Argument \'data\' should not contain a column named \"XXtimeXX\" as this name is used by the lmm function when the argument \'repetition\' is missing. \n")
                    }                    
                    data$XXtimeXX <- "t1"
                    repetition <-  ~XXtimeXX | XXidXX
                }else{
                    repetition <-  stats::as.formula(paste0("~",paste(all.vars(structure$formula),collapse="*")," | XXidXX"))
                }
                if("XXidXX" %in% names(data)){
                    stop("Argument \'data\' should not contain a column named \"XXidXX\" as this name is used by the lmm function when the argument \'repetition\' is missing. \n")
                }
                data$XXidXX <- 1:NROW(data)
            }
        }else{
            stop("Argument \'repetition\' is misisng. \n")
        }
    }
    if(!inherits(repetition,"formula")){
        stop("Argument \'repetition\' must be of class formula, something like: ~ time|id or group ~ time|id. \n")
    }

    name.vcov <- all.vars(repetition)
    if(any(name.vcov %in% names(data) == FALSE)){
        invalid <- name.vcov[name.vcov %in% names(data) == FALSE]
        stop("Argument \'formula\' is inconsistent with argument \'data\'. \n",
             "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
             sep = "")
    }

    ## - right hand side
    if(trace>=2){cat(" (rhs ")}

    if(!grepl("|",deparse(repetition),fixed = TRUE)){
        stop("Incorrect specification of argument \'repetition\'. \n",
             "No | symbol found so no grouping variable could be defined. \n",
             "Shoud be something like: ~ time|id or group ~ time|id. \n")
    }

    if(length(grepl("|",deparse(repetition),fixed = TRUE))>1){
        stop("Incorrect specification of argument \'repetition\'. \n",
             "The symbol | should only appear once, something like: ~ time|id or group ~ time|id. \n")
    }
    res.split <- strsplit(deparse(repetition),"|", fixed = TRUE)[[1]]
    var.cluster <- trimws(res.split[2], which = "both")
    formula.var <- stats::as.formula(res.split[1])
    var.time <- all.vars(stats::update(formula.var,0~.))

    if(length(var.time)!=1){
        stop("Incorrect specification of argument \'repetition\'. \n",
             "There should be exactly one variable before the grouping symbol (|), something like: ~ time|id or group ~ time|id. \n")
    }
    if(length(var.cluster)!=1){
        stop("Incorrect specification of argument \'repetition\'. \n",
             "Should have exactly one variable after the grouping symbol (|), something like: ~ time|id or group ~ time|id. \n")
    }

    test.duplicated <- tapply(data[[var.time]], data[[var.cluster]], function(iT){any(duplicated(iT))})
    if(any(test.duplicated)){
        stop("Incorrect specification of argument \'repetition\'. \n",
             "The time variable (first variable before |) should contain unique values within clusters \n")
    }

    ## *** left hand side
    if(trace>=2){cat(" lhs) ")}
    if(length(lhs.vars(repetition))==0){
        var.strata <- NULL
    }else{
        if(length(lhs.vars(repetition))!=1){
            stop("Incorrect specification of argument \'repetition\'. \n",
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
                iCoef <- which(attr(formula.terms,"factors")[iX,]>=1)
                iInteraction <- attr(formula.terms,"factors")[var.strata,iCoef,drop=FALSE]
                if(length(unique(iInteraction))!=n.strata){
                    stop("When a variable is used to stratify the variance structure, it should also be used to stratify the mean structure. \n",
                         "Consider adding an interaction between \"",iX,"\" and \"",var.strata,"\" in the argument \'formula\'. \n")
                }
            })

            test.length <- tapply(data[[var.time]], data[[var.strata]], function(iT){list(unique(iT))})
            
            if(length(unique(sapply(test.length,length)))>1){
                stop("Incorrect specification of argument \'variance\'. \n",
                    "The time variable should contain the same number of unique values in each strata \n")
            }
            test.unique <- do.call(rbind,lapply(test.length, sort))
            
            if(any(apply(test.unique, MARGIN = 2, FUN = function(x){length(unique(x))})!=1)){
                stop("Incorrect specification of argument \'variance\'. \n",
                     "The time variable should contain the same unique values in each strata \n")
            }

            test.order <- do.call(rbind,test.length)
            
            if(any(apply(test.order, MARGIN = 2, FUN = function(x){length(unique(x))})!=1)){
                stop("Incorrect specification of argument \'variance\'. \n",
                     "The order of the unique values in the time variable should be the same in each strata \n")
            }
        }
    }


    if(trace>=2){cat("\n")}
    
    ## *** data
    if(trace>=2){cat("- data")}
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
    data[[var.time.index]] <- factor(data[[var.time]], levels = unique(data[[var.time]]))  ## to match with gls which chooses the reference level according to the ordering
    ## not not modify var.time in data as it could be also used in the mean structure and that would mess up the ordering 
    U.time <- levels(data[[var.time.index]])
    data[[var.time.index]] <- as.numeric(data[[var.time.index]])

    out$strata <- list(n = n.strata, levels = U.strata, var = var.strata)
    out$time <- list(n = length(U.time), levels = U.time, var = var.time)
    out$cluster <- list(var = var.cluster)
    out$outcome <- list(var = var.outcome)
    if(trace>=2){cat("\n")}
    
    ## *** structure
    if(trace>=2){cat("- structure")}
    if(inherits(structure,"structure")){
        structure <- structure$type
    }else{
        structure <- match.arg(toupper(structure), c("CS","UN","EXP"))
    }
    if(length(out$time$levels)==1 && structure == "UN"){
        warning("Argument \'structure\' has been set to \"UN\" while there is only a single timepoint. \n",
                "Will be change to \"CS\". \n")
        structure  <- "CS"
    }
    if(structure=="UN"){
        formula.cor <- repetition
        if(n.strata==1){
            formula.var  <- formula.var
        }else{
            terms.var <- stats::delete.response(stats::terms(formula.var))
            formula2.var <- stats::update(terms.var, paste0("~0+",var.strata,"+",var.strata,":.")) ## using ".:var.strata" does not work (it gives the same formula - does not invert . var.strata around the : symbol)
        }
    }else if(structure=="CS"){
        if(n.strata>1){
            formula.var <- stats::as.formula(paste0("~0+",var.strata))
            formula.cor <- stats::as.formula(paste0(var.strata,"~1|",var.cluster))
        }else{
            formula.var <- ~1
            formula.cor <- stats::as.formula(paste0("~1|",var.cluster))
        }
    }
    if(n.strata==1){
        txt.data <- "data.fit"
    }else{
        txt.data <- paste0("data.fit[data.fit$",var.strata,"==\"",U.strata,"\",,drop=FALSE]")
    }

    if(n.strata>1){
        term.exclude <- c(var.strata,paste0(setdiff(var.X, var.strata),":",var.strata))
        formula.design <- stats::update(formula, paste0(".~.-",paste0(term.exclude,collapse="-"))) ## no interaction with the strata variable
        out$formula <- list(mean = formula, ## formula will contain all interactions with strata (cf check)
                            mean.design = formula.design,
                            var = repetition,
                            var.design = formula.var,
                            cor = formula.cor)
    }else{
        formula.design <- formula
        out$formula <- list(mean = formula,
                            mean.design = formula.design,
                            var = repetition,
                            var.design = formula.var,
                            cor = formula.cor)
    }

    ## *** type.information
    if(trace>=2){cat("\n")}

    ## ** design matrices
    if(trace>=1){cat("- extract design matrices")}
    out$design <- .model.matrix.lmm(formula.mean = out$formula$mean.design,
                                    formula.var = out$formula$var.design,
                                    data = data, var.outcome = var.outcome,
                                    var.strata = var.strata, U.strata = U.strata,
                                    var.time = var.time, U.time = U.time,
                                    var.cluster = var.cluster,
                                    structure = structure
                                    )

    ## move from design matrix to dataset + update formula (useful when doing baseline adjustment)
    data.fit <- as.data.frame(out$design$X.mean[,colnames(out$design$X.mean),drop=FALSE])
    colnames(data.fit) <- gsub("(Intercept)","Intercept",gsub(":","_",colnames(data.fit), fixed = TRUE), fixed = TRUE)
    
    txt.formula <- tapply(gsub("(Intercept)","Intercept",gsub(":","_",names(out$design$param$strata.mu), fixed = TRUE), fixed = TRUE),out$design$param$strata.mu, function(iStrata){
        paste0(var.outcome, "~ 0 + ",  paste(iStrata, collapse = " + "))
    })

    ## add outcome,strata,time,id to the dataset
    add.col <- unique(c(var.outcome, var.strata, var.time, var.time.index, var.cluster, all.vars(out$formula$vars)))
    if(any(add.col %in% names(data.fit) == FALSE)){
        data.fit <- cbind(data.fit, data[,add.col[add.col %in% names(data.fit) == FALSE],drop=FALSE])
    }
    out$xfactor <- unique(c(stats::.getXlevels(stats::terms(out$formula$mean.design),data),
                            stats::.getXlevels(stats::terms(out$formula$var.design),data)))
    if(trace>=2){cat("\n")}

    ## ** Estimate model parameters
    if(trace>=1){cat("2. Estimate model parameters")}

    if(max(out$design$cluster$nobs)==1){
        if(structure == "CS"){
            txt.gls <- paste0("nlme::gls(",txt.formula,",
                                         method = ",deparse(method.fit),",
                                         data = ",txt.data,", ...)")
        }else if(structure == "UN"){
            form.var <- stats::as.formula(paste0("~1|",var.time))
            txt.gls <- paste0("nlme::gls(",txt.formula,",
                                         weights = nlme::varIdent(form = ",deparse(form.var),"),
                                         method = ",deparse(method.fit),",
                                         data = ",txt.data,", ...)")
        }
    }else{
        if(structure == "CS"){
            form.cor <- stats::as.formula(paste0("~1|",var.cluster))
            txt.gls <- paste0("nlme::gls(",txt.formula,",
                                         correlation = nlme::corCompSymm(form = ",deparse(form.cor),"),
                                         method = ",deparse(method.fit),",
                                         data = ",txt.data,", ...)")

        }else if(structure == "EXP"){
            form.var <- stats::as.formula(paste0("~1|",var.time))
            form.cor <- stats::as.formula(paste0("~",var.time,"|",var.cluster))
            txt.gls <- paste0("nlme::gls(",txt.formula,",
                                          correlation = nlme::corExp(form = ",deparse(form.cor),"),
                                          method = ",deparse(method.fit),",
                                          data = ",txt.data,", ...)")
        }else if(structure == "UN"){
            form.var <- stats::as.formula(paste0("~1|",var.time))
            form.cor <- stats::as.formula(paste0("~",var.time.index,"|",var.cluster))
            txt.gls <- paste0("nlme::gls(",txt.formula,",
                                          correlation = nlme::corSymm(form = ",deparse(form.cor),"),
                                          weights = nlme::varIdent(form = ",deparse(form.var),"),
                                          method = ",deparse(method.fit),",
                                          data = ",txt.data,", ...)")
        }
    }
    out$gls <- stats::setNames(lapply(txt.gls, function(iTxt){eval(parse(text = iTxt))}),
                        U.strata)
    out$gls.call <- lapply(out$gls, function(iM){
        paste0(gsub(",",",\n    ",gsub(" ","",paste(deparse(iM$call), collapse = ""))),"\n")
    })
    out$data <- data
    out$method.fit <- method.fit
    out$structure <- structure

    ## collect parameters
    param.mu <- lapply(U.strata, function(iS){ coef(out$gls[[iS]]) })
    param.sigma <- lapply(U.strata, function(iS){ stats::sigma(out$gls[[iS]]) })
    param.k <- lapply(U.strata, function(iS){  coef(out$gls[[iS]]$modelStruct$varStruct, unconstrained = FALSE) })
    param.rho <- lapply(U.strata, function(iS){ coef(out$gls[[iS]]$modelStruct$corStruct, unconstrained = FALSE) })

    out$param$value <- c(stats::setNames(unlist(param.mu),out$design$param$mu),
                         stats::setNames(unlist(param.sigma), out$design$param$sigma),
                         stats::setNames(unlist(param.k), out$design$param$k),
                         stats::setNames(unlist(param.rho), out$design$param$rho)
                         )

    out$param$strata <- c(unlist(lapply(1:n.strata, function(iS){rep(U.strata[iS], times = length(param.mu[[iS]]))})),
                          unlist(lapply(1:n.strata, function(iS){rep(U.strata[iS], times = length(param.sigma[[iS]]))})),
                          unlist(lapply(1:n.strata, function(iS){rep(U.strata[iS], times = length(param.k[[iS]]))})),
                          unlist(lapply(1:n.strata, function(iS){rep(U.strata[iS], times = length(param.rho[[iS]]))})))

    out$param$type <- c(rep("mu", length(unlist(param.mu))),
                        rep("sigma", length(unlist(param.sigma))),
                        rep("k", length(unlist(param.k))),
                        rep("rho", length(unlist(param.rho))))
    
    names(out$param$strata) <- names(out$param$value)
    names(out$param$type) <- names(out$param$value)

    if(trace>=1){cat("\n")}

    ## ** Reparametrisation
    if(trace>=1){cat("3. Reparametrization \n")}
    name.allcoef <- names(out$param$value)
    
    index.var <- which(out$param$type %in% c("sigma","k","rho"))
    out$reparametrize <- .reparametrize(p = out$param$value[index.var], type = out$param$type[index.var], strata = out$param$strata[index.var], time.levels = out$time$levels,
                                        time.k = out$design$param$time.k, time.rho = out$design$param$time.rho,
                                        Jacobian = TRUE, dJacobian = 2, inverse = FALSE,
                                        transform.sigma = options$transform.sigma,
                                        transform.k = options$transform.k,
                                        transform.rho = options$transform.rho,
                                        transform.names = TRUE)

    if(out$reparametrize$transform==FALSE){
        out$reparametrize$newname <- NULL
        out$reparametrize$Jacobian <- NULL
        out$reparametrize$dJacobian <- NULL
    }else{
        newname.allcoef <- name.allcoef
        newname.allcoef[match(names(out$reparametrize$p),name.allcoef)] <- out$reparametrize$newname
    }

    ## ** Compute partial derivatives regarding the mean and the variance
    if(trace>=1){cat("4. Compute partial derivatives regarding the mean and the variance \n")}

    if(trace>=2){cat("- residuals \n")}
    out$residuals <- out$design$Y - out$design$X.mean %*% out$param$value[colnames(out$design$X.mean)]
    
    if(trace>=2){cat("- Omega \n")}
    out$Omega <- .calc_Omega(object = out$design$X.var, param = out$param$value, keep.interim = TRUE)
    out$OmegaM1 <- lapply(out$Omega,solve)
    
    if(trace>=2){cat("- dOmega \n")}
    out$dOmega <- .calc_dOmega(object = out$design$X.var, param = out$param$value, type = out$param$type, Omega = out$Omega,
                               Jacobian = out$reparametrize$Jacobian)

    if(trace>=2){cat("- d2Omega \n")}
    out$d2Omega <- .calc_d2Omega(object = out$design$X.var, param = out$param$value, type = out$param$type,
                                 Omega = out$Omega, dOmega = out$dOmega, pair = out$design$param$pair.varcoef,
                                 Jacobian = out$reparametrize$Jacobian, dJacobian = out$reparametrize$dJacobian)

    ## ** Compute likelihood derivatives
    if(trace>=1){cat("5. Compute likelihood derivatives \n")}

    if(trace>=2){cat("- log-likelihood \n")}
    out$logLik <- .logLik(X = out$design$X.mean, residuals = out$residuals, precision = out$OmegaM1,
                          index.variance = out$design$X.var$cluster, time.variance = out$design$index.time, index.cluster = out$design$index.cluster, 
                          indiv = FALSE, REML = method.fit=="REML")

    if(trace>=2){cat("- score \n")}
    out$score <- .score(X = out$design$X.mean, residuals = out$residuals, precision = out$OmegaM1, dOmega = out$dOmega,
                        index.variance = out$design$X.var$cluster, time.variance = out$design$index.time, index.cluster = out$design$index.cluster,
                        name.varcoef = out$design$X.var$param, name.allcoef = name.allcoef,
                        indiv = FALSE, REML = method.fit=="REML", effects = c("mean","variance","correlation"))

    if(trace>=2){cat("- information \n")}
    out$information <- .information(X = out$design$X.mean, residuals = out$residuals, precision = out$OmegaM1, dOmega = out$dOmega, d2Omega = out$d2Omega, robust = FALSE,
                                    index.variance = out$design$X.var$cluster, time.variance = out$design$index.time, index.cluster = out$design$index.cluster,
                                    name.varcoef = out$design$X.var$param, name.allcoef = name.allcoef,
                                    pair.meanvarcoef = out$design$param$pair.meanvarcoef, pair.varcoef = out$design$param$pair.varcoef,
                                    indiv = FALSE, REML = method.fit=="REML", type.information = type.information, effects = c("mean","variance","correlation"))
    attr(out$information, "type.information") <- type.information

    if(trace>=2){cat("- variance-covariance \n")}
    out$vcov <- solve(out$information)

    if(df){
        if(trace>=2){cat("- degrees of freedom \n")}
        out$df <- .df(param = out$param, reparametrize = out$reparametrize, Y = out$design$Y, X.mean = out$design$X.mean, X.var = out$design$X.var,
                      index.variance = out$design$X.var$cluster, time.variance = out$design$index.time, index.cluster = out$design$index.cluster,
                      name.varcoef = out$design$X.var$param, 
                      time.k = out$design$param$time.k, time.rho = out$design$param$time.rho,
                      pair.meanvarcoef = out$design$param$pair.meanvarcoef, pair.varcoef = out$design$param$pair.varcoef, REML = (method.fit=="REML"), type.information = type.information, effects = c("mean","variance","correlation"),
                      transform.sigma = out$reparametrize$transform.sigma, transform.k = out$reparametrize$transform.k, transform.rho = out$reparametrize$transform.rho,
                      vcov = out$vcov, diag = TRUE, method.numDeriv = options$method.numDeriv, robust = FALSE)
        out$dVcov <- attr(out$df,"dVcov")
        attr(out$df,"dVcov") <- NULL
    }    

    ## ** convert to lmm and export
    class(out) <- "lmm"
    return(out)
}



######################################################################
### lmm.R ends here
