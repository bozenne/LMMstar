### structure.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:28) 
## Version: 
## Last-Updated: jun  1 2023 (15:56) 
##           By: Brice Ozenne
##     Update #: 1022
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * ID (identity)
##' @title identity Structure
##' @description Variance-covariance structure where the residuals are independent and identically distribution.
##' Can be stratified on a categorical variable.
##' 
##' @param formula formula indicating on which variable to stratify the residual variance (left hand side).
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param add.time not used.
##'
##' @details A typical formula would be \code{~1}.
##'
##' @return An object of class \code{IND} that can be passed to the argument \code{structure} of the \code{lmm} function.
##'
##' @examples
##' ID(NULL, var.cluster = "id", var.time = "time")
##' ID(~1, var.cluster = "id", var.time = "time")
##' ID(~gender, var.cluster = "id", var.time = "time")
##' ID(gender~1, var.cluster = "id", var.time = "time")
##' @export
ID <- function(formula, var.cluster, var.time, add.time){

    ## ** normalize input
    outCov <- .formulaStructure(formula, add.X = NULL, strata.X = TRUE, correlation = FALSE)

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = NA,
                                  cor = NA,
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = NULL),
                class = "ID")

    ## ** export
    class(out) <- append("structure",class(out))
    class(out) <- append("ID",class(out))
    return(out)
}

## * IND (independence)
##' @title Independence Structure
##' @description Variance-covariance structure where the residuals are independent.
##' Can be stratified on a categorical variable.
##'
##' @param formula formula indicating variables influencing the residual variance,
##' using either as a multiplicative factor (right hand side) or stratification (left hand side) to model their effect.
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param add.time Should the default formula (i.e. when \code{NULL}) contain a time effect.
##'
##' @details A typical formula would be either \code{~1} indicating constant variance
##' or \code{~time} indicating a time dependent variance.
##' 
##' @return An object of class \code{IND} that can be passed to the argument \code{structure} of the \code{lmm} function.
##'
##' @examples
##' IND(NULL, var.cluster = "id", var.time = "time", add.time = TRUE)
##' IND(~1, var.cluster = "id", var.time = "time")
##' IND(gender~1, var.cluster = "id", var.time = "time")
##' IND(gender~time, var.cluster = "id", var.time = "time")
##' IND(~gender+time, var.cluster = "id", var.time = "time")
##' 
##' @export
IND <- function(formula, var.cluster, var.time, add.time){

    ## ** normalize input
    if(!missing(add.time)){
        if(is.character(add.time)){
            add.X <- list(variance = add.time,
                          correlation = NULL)
        }else if(add.time){
            add.X <- list(variance = var.time,
                          correlation = NULL)
        }else if(!add.time){
            add.X <- NULL
        }else{
            stop("Incorrect argument \'add.time\': should be logical or character. \n")
        }
    }else{
        add.X <- NULL
    }

    outCov <- .formulaStructure(formula, add.X = add.X, strata.X = FALSE, correlation = FALSE)

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  cor = NA,
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = NULL),
                class = "IND")

    ## ** export
    class(out) <- append("structure",class(out))
    class(out) <- append("IND",class(out))
    return(out)
}


## * CS (compound symmetry)
##' @title Compound Symmetry Structure
##' @description Variance-covariance structure where the residuals have constant variance and correlation overall (default) or within covariate levels.
##' Can be stratified on a categorical variable. 
##'
##' @param formula formula indicating on which variable to stratify the residual variance and correlation (left hand side)
##' and variables influencing the residual variance and correlation (right hand side).
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param type [character] \itemize{
##' \item \code{"ho"}, \code{"homo"}, or \code{"homogeneous"}: constant variance and covariate-specific correlation.
##' Analogous to crossed or nested random effects.
##' \item \code{"he"}, \code{"hetero"}, or \code{"heterogeneous"}: variance and correlation specific to the level of covariates.
##' Can be seen as more flexible crossed or nested random effects model.
##' }
##' @param group [integer vector] grouping of the regressor for the correlation structure.
##' A constant value corresponds to nested random effects (default) and a regressor-specific value to crossed random effects
##' @param add.time not used.
##'
##' @details A typical formula would be \code{~1}, indicating a variance constant over time and the same correlation between all pairs of times.
##'
##' @return An object of class \code{CS} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @examples
##' ## no covariates
##' CS(~1, var.cluster = "id", var.time = "time")
##' CS(gender~1, var.cluster = "id", var.time = "time")
##'
##' ## covariates
##' CS(~time, var.cluster = "id", var.time = "time")
##' CS(gender~time, var.cluster = "id", var.time = "time")
##' CS(list(~time,~1), var.cluster = "id", var.time = "time")
##' CS(list(gender~time,gender~1), var.cluster = "id", var.time = "time")
##' 
##' @export
CS <- function(formula, var.cluster, var.time, type = "homogeneous", group.type = NULL, add.time){

    ## ** normalize input
    type.ho <- c("ho","homo","homogeneous")
    type.he <- c("he","hetero","heterogeneous")
    type <- match.arg(type, c(type.ho,type.he))
    if(type %in% type.ho){
        type <- type.ho[3]
    }else if(type %in% type.he){
        type <- type.he[3]
    }
    if(inherits(formula,"formula") && type == "homogeneous"){
        if(attr(terms(formula),"response")==0){
            formula <- list(variance = ~1,
                            correlation = formula)
        }else{
            formula <- list(variance = update(formula,".~0"),
                            correlation = formula)
        }        
    }
    outCov <- .formulaStructure(formula, add.X = NULL, strata.X = FALSE, correlation = TRUE)
    if(length(outCov$X.cor)==0){
        type <- "homogeneous"
    }
    if(is.null(group.type)){
        if(length(outCov$X.cor)==0){
            group.type <- NULL
        }else{
            group.type <- stats::setNames(rep(1,length(outCov$X.cor)), outCov$X.cor)
        }
    }else{
        if(length(group.type) != length(outCov$X.cor)){
            stop("Argument \'group.type\' should have length ",length(outCov$X.cor),", i.e. one value for each variable. \n")
        }
        if(any(duplicated(names(group.type)))){
            stop("Argument \'group.type\' should no have duplicated names. \n")
        }
        if(!is.null(names(group.type))){
            if(any(names(group.type) %in% outCov$X.cor == FALSE)){
                stop("Argument \'group.type\' should contain regressors for the correlation structure. \n",
                     "Valid values: \"",paste(outCov$X.cor, collapse = "\", \""),"\". \n")
            }
        }else{
            names(group.type) <- outCov$X.cor
        }
    }

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  cor = if(length(outCov$X.cor)>0){I(list(outCov$X.cor))}else{NA},
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = outCov$formula.cor),
                type = type,
                group.type = group.type,
                class = "CS")

    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("CS",class(out))
    return(out)
}

## * RE (random effect)
##' @title Random Effect Structure
##' @description Variance-covariance structure parametrized via random effects.
##' Can be stratified on a categorical variable.
##'
##' @param formula formula indicating on which variable to stratify the residual variance and correlation (left hand side)
##' and variables influencing the residual variance and correlation (right hand side).##' 
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param ranef [list] characteristics of the random effects
##' @param add.time not used.
##'
##' @details A typical formula would be \code{~1}, indicating a variance constant over time and the same correlation between all pairs of times.
##'
##' @return An object of class \code{CS} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @examples
##' RE(~1, var.cluster = "id", var.time = "time")
##' RE(~gender, var.cluster = "id", var.time = "time")
##' RE(gender~(1|id), var.cluster = "id", var.time = "time")
##' 
##' @export
RE <- function(formula, var.cluster, var.time, ranef = NULL, add.time){

    ## ** normalize input
    if(!inherits(formula,"formula")){
        stop("Argument \'formula\' must inherits from formula. \n")
    }

    detail.formula <- formula2var(formula)
    test.ranef <- is.null(ranef)

    if(detail.formula$special=="none"){
        if(length(detail.formula$vars$regressor)>0){
            if(length(detail.formula$vars$response)>0){
                stop("The strata variable in argument \'formula\' should be specified on the left or right hand side (not both). \n")
            }
            var.strata <- detail.formula$vars$regressor
        }else{
            var.strata <- detail.formula$vars$response
        }
    }else if(detail.formula$special=="ranef"){
        var.strata <- detail.formula$vars$response
    }else{
        stop("Incorrect argument \'formula\' for structure RE. \n",
             "Should be something like ~strata or strata ~ (1|id).")
    }

    ## ** create structure
    if(length(var.strata)>0){
        ff.var <- stats::as.formula(paste0(var.strata,"~1"))
    }else{
        ff.var <- ~1
    }

    if(test.ranef && detail.formula$special=="ranef"){
        ranef <- detail.formula$special
    }
    if(is.null(ranef)){
        ## single random intercept
        ff.cor <- ff.var
        group <- NULL
    }else{
        ## multiple random intercepts
        if(any(duplicated(unlist(attr(ranef,"hierarchy"))))){
            stop("Crossed random effects should be defined relative to distinct variables. \n")
        }
        if(!missing(var.cluster) && !is.null(attr(var.cluster,"original"))){
            Ovar.cluster <- attr(var.cluster,"original")
            ## remove cluster level from formula 
            ff.cor <- stats::as.formula(paste0(var.strata,"~",paste0(setdiff(attr(ranef,"vars"),Ovar.cluster), collapse = "+")))
            ## add cluster level when nested random effects
            if((attr(ranef,"crossed") == FALSE) && Ovar.cluster %in% attr(ranef,"vars") == FALSE){
                attr(ranef,"nested") <- TRUE
                attr(ranef,"vars") <- c(Ovar.cluster, attr(ranef,"vars"))
                attr(ranef,"hierarchy") <- lapply(attr(ranef,"hierarchy"), function(iVec){c(Ovar.cluster, iVec)})
            }
        }else{
            ff.cor <- stats::as.formula(paste0(var.strata,"~",paste0(attr(ranef,"vars"), collapse = "+")))
        }
        n.group <- length(attr(ranef,"hierarchy"))
        if(n.group==1){
            group <- NULL
        }else{
            group <- unlist(lapply(1:n.group, function(iG){
                stats::setNames(rep(iG, length(attr(ranef,"hierarchy")[[iG]])), attr(ranef,"hierarchy")[[iG]])
            }))
        }
    }
    out <- CS(list(variance = ff.var, correlation = ff.cor),
              var.cluster = var.cluster, var.time = var.time,
              type = "homogeneous", group.type = group)
    
    if(!test.ranef){
        out$ranef <- list(crossed = attr(ranef,"crossed"),
                          nested = attr(ranef,"nested"),
                          formula = attr(ranef,"formula"),
                          vars = attr(ranef,"vars"),
                          terms = attr(ranef,"terms"),
                          hierarchy = attr(ranef,"hierarchy"))
    }

    ## ** export
    out$call <- match.call()
    out$class <- "RE"
    class(out) <- append("RE",class(out))
    return(out)
}


## * TOEPLITZ (Toeplitz)
##' @title Toeplitz Structure
##' @description Variance-covariance structure for stationnary processes.
##' Can be stratified on a categorical variable.
##'
##' @param formula formula indicating on which variable to stratify the residual variance and correlation (left hand side)
##' and variables influencing the residual variance and correlation (right hand side). 
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param heterogeneous [character] degree of flexibility of the correlation structure within covariate (\code{"UN","LAG","CS"}).
##' Will also affect the variance structure when not explicit.
##' @param add.time Should the default formula (i.e. when \code{NULL}) contain a time effect.
##'
##' @details \bold{formula}: there can only be at most one covariate for the correlation structure.
##' A typical formula would be \code{~1}, indicating a variance constant over time and a correlation specific to each gap time.
##'
##' \bold{heterogeneous}: for a binary covariate the correlation matrix can be decomposed into four blocs: A, B, B, C.
##' A correspond the correlation within level 0 of the covariate, C within level 1, and B between level 0 and 1.
##' Different correlation structures can be specified:\itemize{ 
##' \item \code{"UN"}: unstructured matrix except for the diagonal elements of C which are constrained to be equal.
##' \item \code{"LAG"}: Toeplitz structure within A, B, and C, i.e. correlation specific to each time lag and covariate level.
##' \item \code{"CS"}: block-specific value except for C which has a different value for its diagonal elements.
##'}
##' @return An object of class \code{TOEPLITZ} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @examples
##' ## no covariate
##' TOEPLITZ(~1, var.cluster = "id", var.time = "time")
##' TOEPLITZ(gender~1, var.cluster = "id", var.time = "time")
##' TOEPLITZ(list(~time,~1), var.cluster = "id", var.time = "time")
##' TOEPLITZ(list(gender~time,gender~1), var.cluster = "id", var.time = "time")
##'
##' ## with covariates
##' TOEPLITZ(~side, var.cluster = "id", heterogeneous = "UN",
##'          var.time = "time", add.time = TRUE)
##' TOEPLITZ(~side, var.cluster = "id", heterogeneous = "LAG",
##'          var.time = "time", add.time = TRUE)
##' TOEPLITZ(~side, var.cluster = "id", heterogeneous = "CS",
##'          var.time = "time", add.time = TRUE)
##' TOEPLITZ(gender~side, var.cluster = "id", heterogeneous = "CS",
##'          var.time = "time", add.time = TRUE)
##' @export
TOEPLITZ <- function(formula, var.cluster, var.time, heterogeneous = "LAG", add.time){

    ## ** normalize input
    if(is.null(heterogeneous)){
        heterogeneous <- "LAG"
        toeplitz.block <- FALSE
    }else{
        heterogeneous <- match.arg(heterogeneous, c("UN","LAG","CS"))
        toeplitz.block <- NULL
    }
    
    if(!missing(add.time)){
        if(is.character(add.time)){
            if(is.null(heterogeneous) || heterogeneous=="UN"){
                add.X <- list(variance = add.time,
                              correlation = add.time)
            }else if(heterogeneous %in% c("LAG","CS")){
                if(length(add.time)>1){
                    add.X <- list(variance = utils::tail(add.time,1),
                                  correlation = add.time)
                }else{
                    add.X <- list(variance = NULL,
                                  correlation = add.time)
                }
            }
        }else if(add.time){
            if(is.null(heterogeneous) || heterogeneous=="UN"){
                add.X <- list(variance = var.time,
                              correlation = var.time)
            }else if(heterogeneous %in% c("LAG","CS")){
                if(length(var.time)>1){
                    add.X <- list(variance = utils::tail(var.time,1),
                                  correlation = var.time)
                }else{
                    add.X <- list(variance = NULL,
                                  correlation = var.time)
                }
            }
        }else if(!add.time){
            add.X <- NULL
        }else{
            stop("Incorrect argument \'add.time\': should be logical or character. \n")
        }
    }else{
        add.X <- NULL
    }

    outCov <- .formulaStructure(formula, add.X = add.X, strata.X = FALSE, correlation = TRUE)
    if(length(outCov$X.var)==0 && length(outCov$X.cor)==0){
        heterogeneous <- NULL
    }
    if(length(outCov$X.cor)>2){
        stop("TOEPLITZ covariance structure does not support more than 2 covariates for the correlation structure. \n")
    }else if(is.null(toeplitz.block)){
        toeplitz.block <- length(outCov$X.cor) > 1
    }

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(!is.null(outCov$strata)){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  cor = if(length(outCov$X.cor)>0){I(list(outCov$X.cor))}else{NA},
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = outCov$formula.cor),
                heterogeneous = heterogeneous,
                block = toeplitz.block,
                class = "TOEPLITZ")

    ## ** export
    class(out) <- append("structure",class(out))
    class(out) <- append("TOEPLITZ",class(out))
    return(out)
}

## * UN (unstructured)
##' @title Unstructured Structure 
##' @description Variance-covariance structure where the residuals have time-specific variance and correlation.
##' Can be stratified on a categorical variable.
##'
##' @param formula formula indicating on which variable to stratify the covariance structure.
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param add.time Should the default formula (i.e. when \code{NULL}) contain a time effect.
##'
##' @details A typical formula would be \code{~1}, indicating a time-specific variance parameter and a correlation parameter specific to each pair of times.
##'
##' @return An object of class \code{UN} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @examples
##' UN(NULL, var.cluster = "id", var.time = "time", add.time = TRUE)
##' UN(~gender, var.cluster = "id", var.time = "time", add.time = TRUE)
##' UN(gender ~ 1, var.cluster = "id", var.time = "time", add.time = TRUE)
##' UN(list(~gender,~1), var.cluster = "id", var.time = "time", add.time = TRUE)
##' UN(list(gender~age,gender~1), var.cluster = "id", var.time = "time", add.time = TRUE)
##' 
##' @export
UN <- function(formula, var.cluster, var.time, add.time){

    ## ** normalize input
    if(!missing(add.time)){
        if(is.character(add.time)){
            add.X <- list(variance = add.time,
                          correlation = add.time)
        }else if(add.time){
            add.X <- list(variance = var.time,
                          correlation = var.time)
        }else if(!add.time){
            add.X <- NULL
        }else{
            stop("Incorrect argument \'add.time\': should be logical or character. \n")
        }
    }else{
        add.X <- NULL
    }

    outCov <- .formulaStructure(formula, add.X = add.X, strata.X = inherits(formula,"formula"), correlation = 2) ## 2 for fully stratified structure

    ## ** create structure
    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = I(list(outCov$X.var)),
                                  cor = I(list(outCov$X.cor)),
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = outCov$formula.cor),
                class = "UN")

    ## ** export
    class(out) <- append("structure",class(out))
    class(out) <- append("UN",class(out))
    return(out)
}

## * LV (latent variable)

## * EXP (exponential)
##' @title Exponential Structure
##' @description Variance-covariance structure where the residuals have a correlation decreasing exponentially,
##' Can be stratified on a categorical variable.
##'
##' @param formula formula indicating on which variable to stratify the residual variance and correlation (left hand side)
##' and variables influencing the residual variance and correlation (right hand side).
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param nugget [logical] whether a nugget effect is present.
##' @param add.time not used.
##'
##' @details A typical formula would be \code{~1}, indicating a variance constant over time and correlation with exponential decrease over time.
##'
##' Inspired from \code{nlme::corExp} where if \eqn{K} denotes the nugget effect and \eqn{\rho} the time effect,
##' the correlation between two observations with a time gap \eqn{dt} is \eqn{exp(-\rho dt)} when no nugget effect is present and \eqn{(1-K) exp(-\rho dt)} when a nugget effect is assumed. 
##'
##' @return An object of class \code{EXP} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @examples
##' EXP(var.cluster = "id", var.time = "time", add.time = TRUE)
##' EXP(~space, var.cluster = "id", var.time = "time", add.time = TRUE)
##' EXP(list(~space,~space), var.cluster = "id", var.time = "time", add.time = TRUE)
##' 
##' @export
EXP <- function(formula, var.cluster, var.time, nugget = FALSE, add.time){

    if(missing(formula) || is.null(formula)){
        outCov <- .formulaStructure(list(~1,stats::as.formula(paste0("~",var.time))), heterogeneous = nugget)
    }else if(is.list(formula)){
        outCov <- .formulaStructure(formula, heterogeneous = nugget)
    }else if(!missing(add.time) && (is.character(add.time) || identical(add.time,TRUE)) && length(all.vars(stats::update(formula,0~.)))==0){
        if(is.character(add.time)){
            var.time <- add.time
        }
        if(attr(stats::terms(formula),"response")==1){ # with strata
            ff <- stats::as.formula(paste0(all.vars(formula),"~",var.time))
        }else{
            ff <- stats::as.formula(paste0("~",var.time))
        }
        outCov <- .formulaStructure(list(formula,ff), heterogeneous = nugget)
    }else{
        if(attr(stats::terms(formula),"response")==1){ # with strata
            outCov <- .formulaStructure(list(stats::as.formula(paste0(all.vars(formula),"~1")),formula), heterogeneous = nugget)
        }else{
            outCov <- .formulaStructure(list(~1,formula), heterogeneous = nugget)
        }
    }

    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(!is.null(outCov$strata)){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  cor = if(length(outCov$X.cor)>0){I(list(outCov$X.cor))}else{NA},
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = outCov$formula.cor),
                heterogeneous = nugget,
                class = "EXP")

    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("EXP",class(out))
    return(out)
}


## * CUSTOM (user-specified)
##' @title Custom Structure
##' @description Variance-covariance structure specified by the user.
##' 
##' @param formula formula indicating variables influencing the residual variance and correlation (right hand side).
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param FCT.sigma [function] take as argument the model parameters, time, and design matrix.
##' Output the vector of residuals standard deviations.
##' @param dFCT.sigma [list of vectors] list whose elements are the first derivative of argument \code{FCT.sigma}. 
##' @param d2FCT.sigma [list of vectors] list whose elements are the second derivative of argument \code{FCT.sigma} (no cross-terms).
##' @param init.sigma [numeric vector] initial value for the variance parameters.
##' @param FCT.rho [function] take as argument the model parameters, time, and design matrix.
##' Output the matrix of residuals correlation.
##' @param dFCT.rho [list of matrices] list whose elements are the first derivative of argument \code{FCT.rho}. 
##' @param d2FCT.rho [list of matrices] list whose elements are the second derivative of argument \code{FCT.rho} (no cross-terms).
##' @param init.rho [numeric vector] initial value for the correlation parameters.
##' @param add.time not used.
##'
##' @return An object of class \code{CUSTOM} that can be passed to the argument \code{structure} of the \code{lmm} function.
##'
##' @examples
##' 
##' ## Compound symmetry structure
##' CUSTOM(~1,
##'        FCT.sigma = function(p,time,X){rep(p,length(time))},
##'        init.sigma = c("sigma"=1),
##'        dFCT.sigma = function(p,time,X){list(sigma = rep(1,length(time)))},  
##'        d2FCT.sigma = function(p,time,X){list(sigma = rep(0,length(time)))},  
##'        FCT.rho = function(p,time,X){
##'            matrix(p,length(time),length(time))+diag(1-p,length(time),length(time))
##'        },
##'        init.rho = c("rho"=0.5),
##'        dFCT.rho = function(p,time,X){
##'             list(rho = matrix(1,length(time),length(time))-diag(1,length(time),length(time)))
##'        },
##'        d2FCT.rho = function(p,time,X){list(rho = matrix(0,length(time),length(time)))}
##' )
##' 
##' ## 2 block structure
##' rho.2block <- function(p,time,X){
##'    n.time <- length(time)
##'    rho <- matrix(0, nrow = n.time, ncol = n.time)
##'    rho[1,2] <- rho[2,1] <- rho[4,5] <- rho[5,4] <- p["rho1"]
##'    rho[1,3] <- rho[3,1] <- rho[4,6] <- rho[6,4] <- p["rho2"]
##'    rho[2,3] <- rho[3,2] <- rho[5,6] <- rho[6,5] <- p["rho3"]
##'    rho[4:6,1:3] <- rho[1:3,4:6] <- p["rho4"]
##'    return(rho)
##' }
##' drho.2block <- function(p,time,X){
##'    n.time <- length(time)
##'    drho <- list(rho1 = matrix(0, nrow = n.time, ncol = n.time),
##'                 rho2 = matrix(0, nrow = n.time, ncol = n.time),
##'                 rho3 = matrix(0, nrow = n.time, ncol = n.time),
##'                 rho4 = matrix(0, nrow = n.time, ncol = n.time))   
##'    drho$rho1[1,2] <- drho$rho1[2,1] <- drho$rho1[4,5] <- drho$rho1[5,4] <- 1
##'    drho$rho2[1,3] <- drho$rho2[3,1] <- drho$rho2[4,6] <- drho$rho2[6,4] <- 1
##'    drho$rho3[2,3] <- drho$rho3[3,2] <- drho$rho3[5,6] <- drho$rho3[6,5] <- 1
##'    drho$rho4[4:6,1:3] <- drho$rho4[1:3,4:6] <- 1
##'    return(drho)
##' }
##' d2rho.2block <- function(p,time,X){
##'    n.time <- length(time)
##'    d2rho <- list(rho1 = matrix(0, nrow = n.time, ncol = n.time),
##'                  rho2 = matrix(0, nrow = n.time, ncol = n.time),
##'                  rho3 = matrix(0, nrow = n.time, ncol = n.time),
##'                  rho4 = matrix(0, nrow = n.time, ncol = n.time))   
##'    return(d2rho)
##' }
##'
##' CUSTOM(~variable,
##'        FCT.sigma = function(p,time,X){rep(p,length(time))},
##'        dFCT.sigma = function(p,time,X){list(sigma=rep(1,length(time)))},
##'        d2FCT.sigma = function(p,time,X){list(sigma=rep(0,length(time)))},
##'        init.sigma = c("sigma"=1),
##'        FCT.rho = rho.2block,
##'        dFCT.rho = drho.2block,
##'        d2FCT.rho = d2rho.2block,
##'        init.rho = c("rho1"=0.25,"rho2"=0.25,"rho3"=0.25,"rho4"=0.25))
##' 
##' @export
CUSTOM <- function(formula, var.cluster, var.time,
                   FCT.sigma, dFCT.sigma = NULL, d2FCT.sigma = NULL, init.sigma,
                   FCT.rho, dFCT.rho = NULL, d2FCT.rho = NULL, init.rho, add.time){
    if(is.null(formula)){
        outCov <- .formulaStructure(~1)
    }else{
        outCov <- .formulaStructure(formula)
    }

    if(!is.null(outCov$strata)){
        stop("CUSTOM structures not (yet?) compatible with stratification. \n")
    }

    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(!is.null(outCov$strata)){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  cor = if(length(outCov$X.cor)>0){I(list(outCov$X.cor))}else{NA},
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = outCov$formula.cor),
                class = "CUSTOM")

    ## param
    if(!missing(FCT.sigma)){
        if(any(names(formals(FCT.sigma)) %in% c("p","time","X","...") == FALSE)){
            stop("Incorrect argument in \'FCT.sigma\': can be \"p\", \"time\", or \"X\". \n")
        }
        if(!is.null(dFCT.sigma) && any(names(formals(dFCT.sigma)) %in% c("p","time","X","...") == FALSE)){
            stop("Incorrect argument in \'dFCT.sigma\': can be \"p\", \"time\", or \"X\". \n")
        }
        if(!is.null(d2FCT.sigma) && any(names(formals(d2FCT.sigma)) %in% c("p","time","X","...") == FALSE)){
            stop("Incorrect argument in \'d2FCT.sigma\': can be \"p\", \"time\", or \"X\". \n")
        }
        out$FCT.sigma <- FCT.sigma
        out$dFCT.sigma <- dFCT.sigma
        out$d2FCT.sigma <- d2FCT.sigma
        out$init.sigma <- init.sigma
        if(length(out$init.sigma)>0 && is.null(names(out$init.sigma))){
            names(out$init.sigma) <- paste0("sigma",1:length(out$init.sigma))
        }
    }else{
        out$FCT.sigma <- function(p,time,X){rep(p,length(time))}
        out$dFCT.sigma <- function(p,time,X){list(rep(1,length(time)))}
        out$d2FCT.sigma <- function(p,time,X){list(rep(0,length(time)))}
        out$init.sigma <- c(sigma = NA)
    }
    if(!missing(FCT.rho)){
        if(any(names(formals(FCT.rho)) %in% c("p","time","X","...") == FALSE)){
            stop("Incorrect argument in \'FCT.rho\': can be \"p\", \"time\", or \"X\". \n")
        }
        if(!is.null(dFCT.rho) && any(names(formals(dFCT.rho)) %in% c("p","time","X","...") == FALSE)){
            stop("Incorrect argument in \'dFCT.rho\': can be \"p\", \"time\", or \"X\". \n")
        }
        if(!is.null(d2FCT.rho) && any(names(formals(d2FCT.rho)) %in% c("p","time","X","...") == FALSE)){
            stop("Incorrect argument in \'d2FCT.rho\': can be \"p\", \"time\", or \"X\". \n")
        }
        
        out$FCT.rho <- FCT.rho
        out$dFCT.rho <- dFCT.rho
        out$d2FCT.rho <- d2FCT.rho
        out$init.rho <- init.rho
        if(length(init.rho)>0 && is.null(names(out$init.rho))){
            names(out$init.rho) <- paste0("rho",1:length(out$init.rho))
        }
    }else{
        out$FCT.rho <- NULL
        out$dFCT.rho <- NULL
        out$d2FCT.rho <- NULL
        out$init.rho <- NULL
    }

    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("CUSTOM",class(out))
    return(out)
}

## * helper
## ** .formulaStructure
##' @title Extract Variable From Formula For VCOV Structure
##' @description Extract the variables from the variance and correlation formula to be used to initialize the variance-covariance structure.
##' @noRd
##'
##' @param formula A formula or a list of two formulas.
##' @param add.X [character vector] additional covariates to be added to the variance and correlation structure.
##' @param strata.X [logical] all covariates should be used to stratify the variance and/or correlation structure.
##' @param correlation [logical] should a correlation structure be output? Otherwise no correlation parameter is considered.
##' If greater than 1 then all interaction terms are considered when having multiple variables, e.g. time:gender instead of time + gender.
##' 
##' @keywords internal
##' @examples
.formulaStructure <- function(formula, add.X, strata.X, correlation){

    ## ** normalize to formula format
    if(is.null(formula)){
        formula <- ~1
    }

    class.formula <- inherits(formula,"formula")
    if(class.formula){
        detail.formula0 <- formula2var(formula)

        formula <- list(variance = formula,
                        correlation = formula)
        detail.formula <- list(variance = detail.formula0,
                               correlation = detail.formula0)
    }else if(is.list(formula) && length(formula)==2 && all(sapply(formula,inherits,"formula"))){

        if(is.null(names(formula))){
            names(formula) <- c("variance","correlation")
        }else if(all(names(formula) %in% c("var","cor"))){
            names(formula)[names(formula)=="var"] <- "variance"
            names(formula)[names(formula)=="cor"] <- "correlation"
        }else if(all(names(formula) %in% c("variance","correlation"))){
            formula <- formula[c("variance","correlation")]
        }else{
            stop("Incorrect names associated to the formula for the residual variance-covariance structure. \n",
                 "Should be \"variance\" and \"correlation\" (or \"var\" and \"cor\"). \n")
        }
        detail.formula <- lapply(formula, formula2var)
    }else{
        stop("Incorrect argument \'formula\': should be a formula or a list of 2 formula (var, cor).\n")
    }

    ## ** type
    if(any(sapply(detail.formula, function(iDetail){iDetail$special})!="none")){
        stop("Incorrect argument \'formula\': there should be no special character, e.g. no |. \n")
    }

    ## ** left hand side
    ls.var.strata <- lapply(detail.formula, function(iDetail){iDetail$vars$response})
    if(!identical(ls.var.strata[[1]],ls.var.strata[[2]])){
        stop("Incorrect argument \'formula\': strata variable differ between the correlation and variance structure. \n")
    }
    var.strata <- ls.var.strata[[1]]
    if(length(var.strata)==0){
        var.strata <- NULL
    }else if(length(var.strata)>1){
        stop("There should be at most one strata variable. \n")
    }

    ## ** right hand side
    ls.var.X <- list(variance = c(add.X$variance, detail.formula$variance$vars$regressor),
                     correlation = c(add.X$correlation, detail.formula$correlation$vars$regressor))
    test.interaction <- sapply(formula, function(iF){
        any(attr(stats::delete.response(stats::terms(iF)),"order")>1)
    })
    if(any(test.interaction)){
        stop("Does not handle interactions in the formula. \n")
    }

    ## right hand side variables define strata
    if(length(var.strata)==0 && strata.X && length(unlist(ls.var.X))>0){
        if(length(ls.var.X$variance)>0 && length(ls.var.X$correlation)>0){
            if(!identical(detail.formula$variance$vars$regressor,detail.formula$correlation$vars$regressor)){
                stop("The strata variable should not differ between the variance and the correlation structure. \n")
            }else{
                var.strata <- detail.formula$variance$vars$regressor
                ls.var.X$variance <- setdiff(ls.var.X$variance,detail.formula$variance$vars$regressor)
                ls.var.X$correlation <- setdiff(ls.var.X$correlation,detail.formula$variance$vars$regressor)
            }
        }else if(length(ls.var.X$variance)>0){
            var.strata <- detail.formula$variance$vars$regressor
            ls.var.X$variance <- setdiff(ls.var.X$variance,detail.formula$variance$vars$regressor)                
        }else if(length(ls.var.X$correlation)>0){
            var.strata <- detail.formula$correlation$vars$regressor
            ls.var.X$correlation <- setdiff(ls.var.X$correlation,detail.formula$correlation$vars$regressor)                
        }

    }

    ## ** combine left and right hand side
    ## *** variance
    if(length(ls.var.X$variance)==0){
        X.var <- NULL
        if(length(var.strata)==0){            
            if(attr(stats::terms(formula$variance),"intercept")){
                formula.var <- ~1
            }else{
                formula.var <- ~0
            }
        }else{
            formula.var <- stats::as.formula(paste("~0+",var.strata))
        }
    }else{
        X.var <- unname(ls.var.X$variance)
        if(length(var.strata)==0){
            formula.var <- stats::as.formula(paste("~",paste(ls.var.X$variance,collapse=":")))
        }else{
            formula.var <- stats::as.formula(paste("~0+",var.strata,"+",paste(paste(ls.var.X$variance,collapse=":"),var.strata,sep=":")))
            ## stats::update(terms.var, paste0("~0+",out$name$strata,"+",out$name$strata,":."))
            ## using ".:var.strata" does not work (it gives the same formula - does not invert . var.strata around the : symbol)
        }
    }

    ## *** correlation
    if(length(ls.var.X$correlation)==0){
        X.cor <- NULL
        if(correlation==FALSE){
            formula.cor <- NULL
        }else if(length(var.strata)==0){
            formula.cor <- ~1
        }else{
            formula.cor <- stats::as.formula(paste("~0+",var.strata))
        }
    }else{
        X.cor <- unname(ls.var.X$correlation)
        if(correlation==FALSE){
            formula.cor <- NULL
            if(!class.formula){
                warning("Variable(s) \"",paste(ls.var.X$correlation, collapse = "\", \""),"\" in the correlation structure are ignored. \n")
            }
        }else if(length(var.strata)==0){
            if(correlation>1){
                formula.cor <- stats::as.formula(paste("~0+",paste(ls.var.X$correlation,collapse=":")))
            }else{
                formula.cor <- stats::as.formula(paste("~0+",paste(ls.var.X$correlation,collapse="+")))
            }
        }else{
            if(correlation>1){
                formula.cor <- stats::as.formula(paste("~0+",paste(c(ls.var.X$correlation,var.strata),collapse=":")))
                ## formula.cor <- stats::as.formula(paste("~0+",var.strata,"+",paste(paste(ls.var.X$correlation,var.strata,sep=":"),collapse="+")))
            }else{
                formula.cor <- stats::as.formula(paste("~0+",paste(paste(ls.var.X$correlation,var.strata,sep=":"),collapse="+")))
                ## stats::update(terms.var, paste0("~0+",out$name$strata,"+",out$name$strata,":."))
                ## using ".:var.strata" does not work (it gives the same formula - does not invert . var.strata around the : symbol)
            }
        }
    }

    ## ** export
    out <- list(strata = var.strata,
                X.var = X.var,
                X.cor = X.cor,
                formula.var = formula.var,
                formula.cor = formula.cor
                )
    return(out)
}

##----------------------------------------------------------------------
### structure.R ends here
