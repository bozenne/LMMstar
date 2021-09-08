### lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: sep  8 2021 (13:58) 
##           By: Brice Ozenne
##     Update #: 986
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
##' @param control [glsControl] Control values for gls fit. Passed to gls.
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
##' ## fit Linear Mixed Model
##' eCS.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "CS", data = dL)
##' ## same as
##' eCS.lmm.bis <- lmm(Y ~ X1 + X2 + X5, structure = CS(~visit|id), data = dL)
##'
##' ## output
##' eCS.lmm
##' summary(eCS.lmm)

## * lmm (code)
##' @export
lmm <- function(formula, repetition, structure, data, method.fit = NULL, df = NULL, type.information = NULL, trace = NULL, control = NULL, ...){

    out <- list(call = match.call(), data.original = data)
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
            stop("Argument \'repetition\' is missing. \n")
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
    if(is.factor(data[[var.cluster]])){
        data[[var.cluster]] <- droplevels(data[[var.cluster]])
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
                if(all(iInteraction==0)){
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
    out$data <- data

    ## *** missing values
    var.all <- unname(unique(c(var.strata,var.outcome,var.X,var.time,var.cluster,all.vars(formula.var))))
    index.na <- which(rowSums(is.na(data[,var.all]))>0)
    if(length(index.na)>0){
        out$index.na <- index.na
        attr(out$index.na, "cluster") <- data[[var.cluster]][out$index.na]
        attr(out$index.na, "time") <- data[[var.time]][out$index.na]
        data <- data[-out$index.na,, drop=FALSE]
    }else{
        out$index.na <- NULL
    }
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
    out$structure <- structure

    ## *** formula
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
        var.X.withinStrata <- setdiff(var.X, var.strata)
        if(length(var.X.withinStrata)==0){
            formula.design <- stats::update(formula, paste0(".~1")) ## no interaction with the strata variable
        }else{
            terms.mean <- stats::terms(formula)
            newterm.labels <- gsub(paste0("\\:",var.strata,"$"),"",gsub(paste0("^",var.strata,"\\:"),"",setdiff(attr(terms.mean,"term.labels"),var.strata)))  ## no interaction or main effect with the strata variable
            if(attr(terms.mean,"intercept")>0){
                formula.design <- stats::update(formula, paste0(".~",paste0(unique(newterm.labels),collapse=" + ")))
            }else{
                formula.design <- stats::update(formula, paste0(".~-1+",paste0(unique(newterm.labels),collapse=" + ")))
            }
        }
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

    if(trace>=2){cat("\n")}

    ## ** design matrices
    if(trace>=1){cat("- extract design matrices")}
    out$design <- .model.matrix.lmm(formula.mean = out$formula$mean.design,
                                    formula.var = out$formula$var.design,
                                    data = data, var.outcome = var.outcome,
                                    var.strata = var.strata, U.strata = U.strata,
                                    var.time = var.time, U.time = U.time,
                                    var.cluster = var.cluster,
                                    structure = structure,
                                    precompute.moments = options$precompute.moments
                                    )

    ## move from design matrix to dataset + update formula (useful when doing baseline adjustment)
    data.fit <- as.data.frame(out$design$X.mean[,colnames(out$design$X.mean),drop=FALSE])
    ## colnames(data.fit) <- gsub(" ","_",gsub("(Intercept)","Intercept",gsub(":","_",colnames(data.fit), fixed = TRUE), fixed = TRUE))
    colnames(data.fit) <- paste0("p",1:NCOL(data.fit))
    
    txt.formula <- tapply(paste0("p",1:NCOL(data.fit)),out$design$param$strata.mu, function(iStrata){
        paste0(var.outcome, "~ 0 + ",  paste(iStrata, collapse = " + "))
    })

    ## add outcome,strata,time,id to the dataset
    add.col <- unique(c(var.outcome, var.strata, var.time, var.time.index, var.cluster, all.vars(out$formula$vars)))
    if(any(add.col %in% names(data.fit) == FALSE)){
        data.fit <- cbind(data.fit, data[,add.col[add.col %in% names(data.fit) == FALSE],drop=FALSE])
    }
    out$xfactor <- c(stats::.getXlevels(stats::terms(out$formula$mean.design),data),
                     stats::.getXlevels(stats::terms(out$formula$var.design),data))
    out$xfactor <- out$xfactor[duplicated(names(out$xfactor))==FALSE]
    if(trace>=2){cat("\n")}


    ## ** Estimate model parameters
    if(trace>=1){cat("2. Estimate model parameters")}

    if(options$optimizer=="gls"){
        if(max(out$design$cluster$nobs)==1){
            if(structure == "CS"){
                txt.gls <- paste0("nlme::gls(",txt.formula,",
                                         method = ",deparse(method.fit),",
                                         data = ",txt.data,", control = control)")
            }else if(structure == "UN"){
                form.var <- stats::as.formula(paste0("~1|",var.time))
                txt.gls <- paste0("nlme::gls(",txt.formula,",
                                         weights = nlme::varIdent(form = ",deparse(form.var),"),
                                         method = ",deparse(method.fit),",
                                         data = ",txt.data,", control = control)")
            }
        }else{
            if(structure == "CS"){
                form.cor <- stats::as.formula(paste0("~1|",var.cluster))
                txt.gls <- paste0("nlme::gls(",txt.formula,",
                                         correlation = nlme::corCompSymm(form = ",deparse(form.cor),"),
                                         method = ",deparse(method.fit),",
                                         data = ",txt.data,", control = control)")

            }else if(structure == "EXP"){
                form.var <- stats::as.formula(paste0("~1|",var.time))
                form.cor <- stats::as.formula(paste0("~",var.time,"|",var.cluster))
                txt.gls <- paste0("nlme::gls(",txt.formula,",
                                          correlation = nlme::corExp(form = ",deparse(form.cor),"),
                                          method = ",deparse(method.fit),",
                                          data = ",txt.data,", control = control)")
            }else if(structure == "UN"){
                form.var <- stats::as.formula(paste0("~1|",var.time))
                form.cor <- stats::as.formula(paste0("~",var.time.index,"|",var.cluster))
                txt.gls <- paste0("nlme::gls(",txt.formula,",
                                          correlation = nlme::corSymm(form = ",deparse(form.cor),"),
                                          weights = nlme::varIdent(form = ",deparse(form.var),"),
                                          method = ",deparse(method.fit),",
                                          data = ",txt.data,", control = control)")
            }
        }
        out$gls <- stats::setNames(lapply(txt.gls, function(iTxt){eval(parse(text = iTxt))}),
                                   U.strata)
        out$gls.call <- lapply(out$gls, function(iM){
            paste0(gsub(",",",\n    ",gsub(" ","",paste(deparse(iM$call), collapse = ""))),"\n")
        })

        param.mu <- lapply(U.strata, function(iS){ coef(out$gls[[iS]]) })
        param.sigma <- lapply(U.strata, function(iS){ stats::sigma(out$gls[[iS]]) })
        param.k <- lapply(U.strata, function(iS){  coef(out$gls[[iS]]$modelStruct$varStruct, unconstrained = FALSE) })
        param.rho <- lapply(U.strata, function(iS){ coef(out$gls[[iS]]$modelStruct$corStruct, unconstrained = FALSE) })

        param.value <- c(stats::setNames(unlist(param.mu),out$design$param$mu),
                         stats::setNames(unlist(param.sigma), out$design$param$sigma),
                         stats::setNames(unlist(param.k), out$design$param$k),
                         stats::setNames(unlist(param.rho), out$design$param$rho)
                         )

    }else if(options$optimizer=="FS"){

        outEstimate <- .estimate(design = out$design, time = out$time, method.fit = method.fit, type.information = type.information,
                                 transform.sigma = options$transform.sigma, transform.k = options$transform.k, transform.rho = options$transform.rho,
                                 precompute.moments = options$precompute.moments,
                                 init = control$init, n.iter = control$n.iter, tol.score = control$tol.score, tol.param = control$tol.param, trace = control$trace)
        param.value <- outEstimate$estimate
        out$opt <- outEstimate[c("cv","n.iter","score","previous.estimate")]
        
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
    class(out) <- "lmm"
    return(out)
}



######################################################################
### lmm.R ends here
