### lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: okt  1 2021 (17:00) 
##           By: Brice Ozenne
##     Update #: 1173
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
##' @param control [list] Control values for the optimization method. The element \code{optimizer} indicates which optimizer to use and additional argument will be pass to the optimizer.
##'
##' @details \bold{Computation time} the \code{lmm} has not been developped to be a fast function as, by default, it uses REML estimation with the observed information matrix and uses a Satterthwaite approximation to compute degrees of freedom (this require to compute the third derivative of the log-likelihood which is done by numerical differentiation). The computation time can be substantially reduced by using ML estimation with the expected information matrix and no calculation of degrees of freedom: arguments \code{method.fit="ML"}, \code{type.information="expected"}, \code{df=FALSE}. This will, however, lead to less accurate p-values and confidence intervals in small samples.
##'
##' @return an object of class \code{lmm} containing the estimated parameter values, the residuals, and relevant derivatives of the likelihood.
##' Compatible with standard methods such as \code{summary}, \code{autoplot}, \code{confint}, \code{coef}, \code{anova}, \code{predict}, \code{residuals}.


## * lmm (examples)
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit Linear Mixed Model
##' eCS.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "CS", data = dL)
##' eCS.lmm
##' summary(eCS.lmm)

## * lmm (code)
##' @export
lmm <- function(formula, repetition, structure, data, method.fit = NULL, df = NULL, type.information = NULL, trace = NULL, control = NULL){

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

        if(inherits(structure,"structure")){
            var.cluster <- structure$name$cluster
            var.time <- structure$name$time
            if(is.null(var.time) && structure$type %in% c("CS","UN")){
                stop("Could not identify the time variable based on the \'structure\' argument. \n",
                     "Consider specifying the \'repetition\' argument. \n")
            }
            
        }else if(missing(structure)){

            var.cluster <- NULL
            var.time <- NULL
            structure <- "IND"

        }else{
            stop("Argument \'repetition\' is missing. \n")
        }

        if(is.null(var.time)){
            if("XXtimeXX" %in% names(data)){
                stop("Argument \'data\' should not contain a column named \"XXtimeXX\" as this name is used by the lmm function when the argument \'repetition\' is missing. \n")
            }
            data$XXtimeXX <- "t1"
            var.time <- "XXtimeXX"
        }
        if(is.null(var.cluster)){
            if("XXidXX" %in% names(data)){
                stop("Argument \'data\' should not contain a column named \"XXidXX\" as this name is used by the lmm function when the argument \'repetition\' is missing. \n")
            }
            data$XXidXX <- 1:NROW(data)
            var.cluster <- "XXidXX"
        }
            
        repetition <-  stats::as.formula(paste0("~",var.time," | ",var.cluster))            
    }else{
        if(!inherits(repetition,"formula")){
            stop("Argument \'repetition\' must be of class formula, something like: ~ time|id or group ~ time|id. \n")
        }
        res.split <- strsplit(deparse(repetition),"|", fixed = TRUE)[[1]]
        if(length(res.split)!=2){
            stop("Incorrect specification of argument \'repetition\'. \n",
                 "The symbol | should only exacly once, something like: ~ time|id or group ~ time|id. \n")
        }
        var.cluster <- trimws(res.split[2], which = "both")
        if(length(var.cluster)!=1){
            stop("Incorrect specification of argument \'repetition\'. \n",
                 "Should have exactly one variable after the grouping symbol (|), something like: ~ time|id or group ~ time|id. \n")
        }
        var.time <- all.vars(stats::update(stats::as.formula(res.split[1]),0~.)) 
        if(length(var.time)!=1){
            stop("Incorrect specification of argument \'repetition\'. \n",
                 "There should be exactly one variable before the grouping symbol (|), something like: ~ time|id or group ~ time|id. \n")
        }
    }

    name.vcov <- all.vars(repetition)
    if(any(name.vcov %in% names(data) == FALSE)){
        invalid <- name.vcov[name.vcov %in% names(data) == FALSE]
        if("repetition" %in% names(out$call)){
            stop("Argument \'repetition\' is inconsistent with argument \'data\'. \n",
                 "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
                 sep = "")
        }else{
            stop("Argument \'structure\' is inconsistent with argument \'data\'. \n",
                 "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
                 sep = "")
        }
    }
    if(is.factor(data[[var.cluster]])){
        data[[var.cluster]] <- droplevels(data[[var.cluster]])
    }else{
        data[[var.cluster]] <- factor(data[[var.cluster]], levels = sort(unique(data[[var.cluster]])))
    }
    test.duplicated <- tapply(data[[var.time]], data[[var.cluster]], function(iT){any(duplicated(iT))})
    if(any(test.duplicated)){
        stop("Incorrect specification of argument \'repetition\'. \n",
             "The time variable (first variable before |) should contain unique values within clusters \n")
    }


    if(trace>=2){cat("\n")}

    ## *** structure (residual variance-covariance structure)
    if(trace>=2){cat("- residual variance-covariance structure  ")}

    if(missing(structure)){
        structure <- "UN"
    }
    if(inherits(structure,"structure")){
        if(optimizer=="gls"){
            stop("When using \"gls\" optimizer, the structure should be specified as a character. \n",
                 "Available structures: \"ID\",\"IND\",\"CS\",\"UN\". \n")
        }
        type.structure <- structure$type
    }else if(inherits(structure,"character")){
        type.structure <- match.arg(structure, c("ID","IND","CS","UN"))
        structure <- do.call(type.structure, list(formula = repetition, var.cluster = var.cluster, var.time = var.time))
    }else{
        stop("Argument \'structure\' must either be a character or a structure object. \n")
    }
    var.strata <- structure$name$strata

    ## *** formula (mean structure)
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

    if(!is.na(var.strata) && optimizer == "gls"){
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
    out$formula <- list(mean = formula,
                        mean.design = formula.design,
                        var.design = structure$formula$var,
                        cor.design = structure$formula$cor)
    var.Z <- c(all.vars(out$formula$var.design),all.vars(out$formula$var.design))
    if(trace>=2){cat("\n")}

    ## *** data
    if(trace>=2){cat("- data")}
    ## index
    if("XXindexXX" %in% names(data)){
        stop("Incorrect specification of argument \'data\'. \n",
             "The variable \"XXindexXX\" is used internally but already exists in \'data\' \n")
    }
    data$XXindexXX <- 1:NROW(data)

    ## time
    var.time.index <- paste0("XX",var.time,".indexXX")
    if(var.time.index %in% names(data)){
        stop("Incorrect specification of argument \'data\'. \n",
             "The variable ",var.time.index," is used internally but already exists in \'data\' \n")
    }
    if(is.factor(data[[var.time]])){
        data[[var.time]] <- droplevels(data[[var.time]])
    }else{
        data[[var.time]] <- factor(data[[var.time]], levels = sort(unique(data[[var.time]])))
    }
    U.time <- levels(data[[var.time]])
    data[[var.time.index]] <- as.numeric(data[[var.time]])

    ## strata
    if(is.na(var.strata)){
        var.strata <- "XXstrata.indexXX"
        U.strata <- 1
        n.strata <- 1
        if("XXstrata.indexXX" %in% names(data)){
            stop("Incorrect specification of argument \'data\'. \n",
                 "The variable \"XXstrata.indexXX\" is used internally but already exists in \'data\' \n")
        }
        data$XXstrata.indexXX <- 1
    }else{
        if(is.factor(data[[var.strata]])){
            data[[var.strata]] <- droplevels(data[[var.strata]])
        }else{
            data[[var.strata]] <- factor(data[[var.strata]], levels = sort(unique(data[[var.strata]])))
        }
        U.strata <- levels(data[[var.strata]])
        n.strata <- length(U.strata)
        
        if(any(rowSums(table(data[[var.cluster]],data[[var.strata]])>0)!=1)){
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

        test.length <- tapply(data[[var.time]], data[[var.strata]], function(iT){list(unique(iT))})
            
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

    ## store
    out$strata <- list(n = n.strata, levels = U.strata, var = var.strata)
    out$time <- list(n = length(U.time), levels = U.time, var = var.time)
    out$cluster <- list(var = var.cluster)
    out$outcome <- list(var = var.outcome)
    out$data <- data

    ## *** missing values
    var.all <- unname(unique(c(var.strata,var.outcome,var.X,var.time,var.cluster,var.Z)))
    index.na <- which(rowSums(is.na(data[,var.all]))>0)
    if(length(index.na)>0){
        attr(index.na, "cluster") <- data[[var.cluster]][index.na]
        attr(index.na, "time") <- data[[var.time]][index.na]
        data <- data[-index.na,, drop=FALSE]
    }else{
        index.na <- NULL
    }

    out$index.na <- index.na

    if(trace>=2){cat("\n")}

    ## ** design matrices
    if(trace>=1){cat("- extract design matrices")}
    out$design <- .model.matrix.lmm(formula.mean = out$formula$mean.design,
                                    structure = structure,
                                    data = data, var.outcome = var.outcome,
                                    var.strata = var.strata, U.strata = U.strata,
                                    var.time = var.time, U.time = U.time,
                                    var.cluster = var.cluster,
                                    precompute.moments = options$precompute.moments)

    out$xfactor <- c(stats::.getXlevels(stats::terms(out$formula$mean.design),data),
                     stats::.getXlevels(stats::terms(out$formula$var.design),data))
    if(!is.null(out$formula$cor.design)){
        out$xfactor <- c(out$xfactor,stats::.getXlevels(stats::terms(out$formula$cor.design),data))
    }
    out$xfactor <- out$xfactor[duplicated(names(out$xfactor))==FALSE]

    if(optimizer=="gls"){
        ## move from design matrix to dataset (useful when doing baseline adjustment)
        data.fit <- as.data.frame(out$design$mean)
        ## colnames(data.fit) <- gsub(" ","_",gsub("(Intercept)","Intercept",gsub(":","_",colnames(data.fit), fixed = TRUE), fixed = TRUE))
        colnames(data.fit) <- paste0("p",1:NCOL(data.fit))
    
        ## add outcome,strata,time,id to the dataset
        add.col <- unique(c(var.outcome, var.strata, var.time, var.time.index, var.cluster, all.vars(out$formula$var.design),all.vars(out$formula$cor.design)))
        if(any(add.col %in% names(data.fit) == FALSE)){
            data.fit <- cbind(data.fit, data[,add.col[add.col %in% names(data.fit) == FALSE],drop=FALSE])
        }
        ## order by strata, time, and cluster (strata, cluster, and time does not provide satisfactory results, mixing k-parameters)
        data.fit <- data.fit[order(data[[var.strata]],data[[var.time]],data[[var.cluster]]),,drop=FALSE]

        if(n.strata==1){
            txt.data <- "data.fit"
        }else{
            txt.data <- paste0("data.fit[data.fit$",var.strata,"==\"",U.strata,"\",,drop=FALSE]")
        }
        ##  update formula
        txt.formula <- tapply(paste0("p",1:NCOL(out$design$mean)),out$design$param$strata.mu, function(iStrata){
            paste0(var.outcome, "~ 0 + ",  paste(iStrata, collapse = " + "))
        })
    }

    if(trace>=2){cat("\n")}


    ## ** Estimate model parameters
    if(trace>=1){cat("2. Estimate model parameters")}

    if(optimizer=="gls"){

        if(max(out$design$cluster$nobs)==1 || type.structure == "IND"){
            name.var <- stats::na.omit(setdiff(structure$name$var,structure$name$strata))

            if(type.structure == "CS" || length(name.var)==0){
                txt.gls <- paste0("nlme::gls(",txt.formula,",
                                         method = ",deparse(method.fit),",
                                         data = ",txt.data,", control = control)")
            }else if(type.structure == "UN" || length(name.var)>0){
                form.var <- stats::as.formula(paste0("~1|",var.time))
                txt.gls <- paste0("nlme::gls(",txt.formula,",
                                         weights = nlme::varIdent(form = ",deparse(form.var),"),
                                         method = ",deparse(method.fit),",
                                         data = ",txt.data,", control = control)")
            }
        }else{
            if(type.structure == "CS"){
                form.cor <- stats::as.formula(paste0("~1|",var.cluster))
                txt.gls <- paste0("nlme::gls(",txt.formula,",
                                         correlation = nlme::corCompSymm(form = ",deparse(form.cor),"),
                                         method = ",deparse(method.fit),",
                                         data = ",txt.data,", control = control)")

            }else if(type.structure == "EXP"){
                form.var <- stats::as.formula(paste0("~1|",var.time))
                form.cor <- stats::as.formula(paste0("~",var.time,"|",var.cluster))
                txt.gls <- paste0("nlme::gls(",txt.formula,",
                                          correlation = nlme::corExp(form = ",deparse(form.cor),"),
                                          method = ",deparse(method.fit),",
                                          data = ",txt.data,", control = control)")
            }else if(type.structure == "UN"){
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
        param.value <- c(stats::setNames(unlist(param.mu),out$design$param$mu),
                         stats::setNames(unlist(param.sigma), out$design$param$sigma))
        if(length(out$design$param$k)>0){
            param.k <- lapply(U.strata, function(iS){  coef(out$gls[[iS]]$modelStruct$varStruct, unconstrained = FALSE) })
            param.value <- c(param.value,stats::setNames(unlist(param.k), out$design$param$k))
        }
        if(length(out$design$param$rho)>0){
            param.rho <- lapply(U.strata, function(iS){ coef(out$gls[[iS]]$modelStruct$corStruct, unconstrained = FALSE) })
            param.value <- c(param.value,stats::setNames(unlist(param.rho), out$design$param$rho))
        }        

    }else{
        outEstimate <- .estimate(design = out$design, time = out$time, method.fit = method.fit, type.information = type.information,
                                 transform.sigma = options$transform.sigma, transform.k = options$transform.k, transform.rho = options$transform.rho,
                                 precompute.moments = options$precompute.moments,
                                 optimizer = optimizer, init = control$init, n.iter = control$n.iter, tol.score = control$tol.score, tol.param = control$tol.param, trace = control$trace)
        param.value <- outEstimate$estimate
        out$opt <- outEstimate[c("cv","n.iter","score","previous.estimate")]

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
    class(out) <- "lmm"
    return(out)
}



######################################################################
### lmm.R ends here
