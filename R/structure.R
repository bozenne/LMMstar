### structure.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:28) 
## Version: 
## Last-Updated: May 30 2022 (22:47) 
##           By: Brice Ozenne
##     Update #: 604
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .formulaStructure
##' @title Extract Variable From Formula For VCOV Structure
##' @description Extract the variables from the variance and correlation formula to be used to initialize the variance-covariance structure.
##' @noRd
##'
##' @param formula A formula or a list of two formulas.
##' @param add.X additional covariates to be added to the variance and correlation structure.
##' @param heterogeneous when \code{FALSE}, main effects are kept in the correlation structure.
##' 
##' @keywords internal
##' @examples
##' .formulaStructure(strata ~ time)
##' .formulaStructure( ~ time)
##' .formulaStructure(list( ~ gender+time,  ~ time))
##' .formulaStructure(strata ~ 1)
.formulaStructure <- function(formula, add.X = NULL, heterogeneous = TRUE){

    ## ** normalize to formula format
    if(!is.list(formula)){
        formula <- list(variance = formula,
                        correlation = formula)
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
    }else if(!inherits(formula,"formula")){
        stop("Incorrect argument \'formula\': should be a formula or a list of 2 formula (var, cor).\n")
    }

    ## ** left hand side
    ls.var.strata <- lapply(formula,lhs.vars)
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
    ls.var.X <- lapply(formula, function(iF){unique(c(rhs.vars(iF),add.X))})

    test.interaction <- sapply(formula, function(iF){
        any(attr(stats::delete.response(stats::terms(iF)),"order")>1)
    })
    if(any(test.interaction)){
        stop("Does not handle interactions in the formula. \n")
    }

    ## ** combine left and right hand side
    ## *** variance
    if(length(ls.var.X$variance)==0){
        if(length(var.strata)==0){
            formula.var <- ~1 
        }else{
            formula.var <- stats::as.formula(paste("~0+",var.strata))
        }
    }else{
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
        if(length(var.strata)==0){
            formula.cor <- ~1 
        }else{
            formula.cor <- stats::as.formula(paste("~0+",var.strata))
        }
    }else{
        if(length(var.strata)==0){
            formula.cor <- stats::as.formula(paste("~",paste(ls.var.X$correlation,collapse="+")))
        }else if(heterogeneous){
            formula.cor <- stats::as.formula(paste("~0+",paste(paste(ls.var.X$correlation,var.strata,sep=":"),collapse="+")))
            ## stats::update(terms.var, paste0("~0+",out$name$strata,"+",out$name$strata,":."))
            ## using ".:var.strata" does not work (it gives the same formula - does not invert . var.strata around the : symbol)
        }else{
            formula.cor <- stats::as.formula(paste("~0+",var.strata,"+",paste(paste(ls.var.X$correlation,var.strata,sep=":"),collapse="+")))
        }
    }
    
    ## ** export
    out <- list(strata = var.strata,
                X.var = unname(ls.var.X$variance),
                X.cor = unname(ls.var.X$correlation),
                formula.var = formula.var,
                formula.cor = formula.cor
                )
    return(out)
}




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

    if(is.null(formula) || length(all.vars(formula))==0){
        outCov <- .formulaStructure(~1)
    }else{
        ## put covariate as strata if in addition to time
        outCov <- .formulaStructure(stats::as.formula(paste0(paste(all.vars(formula), collapse="+"),"~1")))
    }
    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(length(outCov$strata)>0){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = NA,
                                  cor = NA,
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = NULL),
                heterogeneous = TRUE,
                type = "ID")

    ## export
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
##' 
##' IND(~time, var.cluster = "id", var.time = "time")
##' IND(gender~time, var.cluster = "id", var.time = "time")
##' IND(~time+gender, var.cluster = "id", var.time = "time")
##' @export
IND <- function(formula, var.cluster, var.time, add.time){

    if(!missing(add.time)){
        if(is.character(add.time)){
            add.X <- add.time
        }else if(add.time){
            add.X <- var.time
        }else if(!add.time){
            add.X <- NULL
        }else{
            stop("Incorrect argument \'add.time\': should be logical or character. \n")
        }
    }else{
        add.X <- NULL
    }

    if(is.null(formula)){
        outCov <- .formulaStructure(~1, add.X = add.X)
    }else{
        outCov <- .formulaStructure(formula, add.X = add.X)
    }

    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(!is.null(outCov$strata)){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  cor = NA,
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = NULL),
                heterogeneous = TRUE,
                type = "IND")

    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("IND",class(out))
    return(out)
}


## * CS (compound symmetry)
##' @title Compound Symmetry Structure
##' @description Variance-covariance structure where the residuals have constant variance and correlation.
##' Can be stratified on a categorical variable.
##'
##' @param formula formula indicating on which variable to stratify the residual variance and correlation (left hand side)
##' and variables influencing the residual variance (right hand side).
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param heterogeneous [logical] when covariates are used for the correlation structure,
##' should correlation parameters should be specific to each level of the covariate?
##' @param add.time not used.
##'
##' @details A typical formula would be \code{~1}, indicating a variance constant over time and the same correlation between all pairs of times.
##'
##' @return An object of class \code{CS} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @examples
##' CS(~1, var.cluster = "id", var.time = "time")
##' CS(gender~1, var.cluster = "id", var.time = "time")
##' CS(list(~time,~1), var.cluster = "id", var.time = "time")
##' CS(list(gender~time,gender~1), var.cluster = "id", var.time = "time")
##' 
##' @export
CS <- function(formula, var.cluster, var.time, heterogeneous = TRUE, add.time){
    if(is.list(formula)){
        outCov <- .formulaStructure(formula, heterogeneous = heterogeneous)
    }else if(is.null(formula)){
        outCov <- .formulaStructure(~1, heterogeneous = heterogeneous)
    }else if(heterogeneous){
        outCov <- .formulaStructure(formula, heterogeneous = heterogeneous)
    }else{
        if(attr(stats::terms(formula),"response")==1){
            outCov <- .formulaStructure(list(stats::update(formula,.~0),formula), heterogeneous = heterogeneous)
        }else{
            outCov <- .formulaStructure(list(~1,formula), heterogeneous = heterogeneous)
        }
    }
    if(length(outCov$X.var)==0 && length(outCov$X.cor)==0){
        heterogeneous <- FALSE
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
                heterogeneous = heterogeneous,
                type = "CS")

    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("CS",class(out))
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
##' 
##' @export
UN <- function(formula, var.cluster, var.time, add.time){

    if(!missing(add.time)){
        if(is.character(add.time)){
            add.X <- add.time
        }else if(add.time){
            add.X <- var.time
        }else if(!add.time){
            add.X <- NULL
        }else{
            stop("Incorrect argument \'add.time\': should be logical or character. \n")
        }
    }else{
        add.X <- NULL
    }
    
    if(is.null(formula) || length(all.vars(formula))==0){
        outCov <- .formulaStructure(~1, add.X = add.X)
    }else{
        outCov <- .formulaStructure(stats::as.formula(paste0(paste(all.vars(formula), collapse="+"),"~1")),
                                    add.X = add.X)
    }

    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(!is.null(outCov$strata)){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  cor = I(list(outCov$X.cor)),
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = outCov$formula.cor),
                heterogeneous = TRUE,
                type = "UN")

    ## remove intercept
    if(length(all.vars(out$formula$cor)>0)){
        out$formula$cor <- stats::update(out$formula$cor,~0+.)
    }

    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("UN",class(out))
    return(out)
}

## * EXP (exponential)

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
                heterogeneous = TRUE,
                type = "CUSTOM")

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


##----------------------------------------------------------------------
### structure.R ends here
