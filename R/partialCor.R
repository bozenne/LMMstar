### partialCor.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May  1 2022 (17:01) 
## Version: 
## Last-Updated: Mar 10 2024 (15:59) 
##           By: Brice Ozenne
##     Update #: 528
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * partialCor (documentation)
##' @title Partial Correlation
##' @name partialCor
##' @description Estimate the partial correlation based on equation 19 of Lloyd et al 2008 (\code{partialCor.lmm}) or explicitely modeling the correlation via a linear mixed model  (\code{partialCor.list}, \code{partialCor.formula}).
##' The first option is numerically more efficient and exact with a single observation per cluster.
##' With multiple repetitions, what is being estimated with the first option may not be clear and the second option is therefore preferrable.
##' 
##' @param object a formula with in the left hand side the variables for which the correlation should be computed
##' and on the right hand side the adjustment set. Can also be a list of formula for outcome-specific adjustment set.
##' @param data [data.frame] dataset containing the variables.
##' @param repetition [formula] Specify the structure of the data: the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' @param structure [character] Specify the residual variance-covariance structure.
##' Without repetitions, either \code{"UN"} or \code{"CS"}.
##' With repetitions, one of \code{"UN"}, \code{"PEARSON"}, \code{"HLAG"}, \code{"LAG"}, \code{"HCS"}, \code{"CS"}.
##' @param by [character] variable used to stratified the correlation on.
##' @param effects [character or matrix] type of contrast to be used for comparing the correlation parameters. One of \code{"Dunnett"}, \code{"Tukey"}, \code{"Sequen"}, or a contrast matrix.
##' @param rhs [numeric vector] right hand side for the comparison of correlation parameters. 
##' @param method [character] adjustment for multiple comparisons (e.g. \code{"single-step"}).
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param R2 [logical] Should the R2 (coefficient of determination) be computed?
##' @param se [logical] Should the uncertainty about the partial correlation be evaluated? Only relevant for \code{partialCor.lmm}.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the coefficient. Otherwise a normal distribution is used.
##' @param transform.rho [character] scale on which perform statistical inference (e.g. \code{"atanh"})
##' @param name.short [logical vector of length 2] use short names for the output coefficients (omit the name of the by variable, omit name of the correlation parameter)
##' @param ... arguments passed to \code{confint} for \code{partialCor.list} and  \code{partialCor.formula}. Not used for \code{partialCor.lmm}.
##'
##' @details Fit a mixed model to estimate the partial correlation with the following variance-covariance pattern:
##' \itemize{
##' \item \bold{no repetition}: unstructure or compound symmetry structure for M observations, M being the number of variables on the left hand side (i.e. outcomes).
##' \item \bold{repetition}: structure for M*T observations where M being the number of variables (typically 2) and T the number of repetitions. Can be \itemize{
##'       \item \code{"UN"}: unstructured (except the off-diagonal containing the correlation parameter which is constant).
##'       \item \code{"PEARSON"}: same as unstructured except it only uses a single variance parameter per variable, i.e. it assumes constant variance over repetitions.
##'       \item \code{"HLAG"}: toeplitz by block with variable and repetition specific variance.
##'       \item \code{"LAG"}: toeplitz by block, i.e. correlation depending on the gap between repetitions and specific to each variable. It assumes constant variance over repetitions.
##'       \item \code{"HCS"}: heteroschedastic compound symmetry by block, i.e. variable specific correlation constant over repetitions. A specific parameter is used for the off-diagonal crossing the variables at the same repetition (which is the marginal correlation parameter).
##'       \item \code{"CS"}: compound symmetry by block. It assumes constant variance and correlation over repetitions.
##' }}
##'
##' @return A data.frame with the estimate partial correlation (rho), standard error, degree of freedom, confidence interval, and p-value (test of no correlation).
##' When \code{structure="CS"} or \code{structure="HCS"} is used with repeated measurements, a second correlation coefficient (r) is output where the between subject variance has been removed (similar to Bland et al. 1995).
##'
##' @references
##'  Bland J M, Altman D G. Statistics notes: Calculating correlation coefficients with repeated observations: Part 1—correlation within subjects BMJ 1995; 310 :446 doi:10.1136/bmj.310.6977.446
##'  Edwards, L.J., Muller, K.E., Wolfinger, R.D., Qaqish, B.F. and Schabenberger, O. (2008), An R2 statistic for fixed effects in the linear mixed model. Statist. Med., 27: 6137-6157. https://doi.org/10.1002/sim.3429
##' 
##' @keywords models
##' 
##' @examples
##' 
##' #### no repetition ####
##' 
##' ## example from ppcor::pcor 
##' y.data <- data.frame(
##'   hl=c(7,15,19,15,21,22,57,15,20,18),
##'   disp=c(0.000,0.964,0.000,0.000,0.921,0.000,0.000,1.006,0.000,1.011),
##'   deg=c(9,2,3,4,1,3,1,3,6,1),
##'   BC=c(1.78e-02,1.05e-06,1.37e-05,7.18e-03,0.00e+00,0.00e+00,0.00e+00
##',  4.48e-03,2.10e-06,0.00e+00)
##')
##' ## ppcor::pcor(y.data)
##'
##' ## partial correlation based on a formula
##' partialCor(c(hl,disp)~BC+deg, data = y.data)
##' partialCor(hl + disp~BC+deg, data = y.data)
##' ## partial correlation based on a list
##' partialCor(list(hl~BC+deg,disp~BC+deg), data = y.data)
##' ## via an existing model
##' e.lm <- lmm(hl~disp+BC+deg, data = y.data)
##' partialCor(e.lm)
##'
##' ## using a different set of covariates for outcome
##' partialCor(list(hl~BC+deg, disp~BC), data = y.data)
##' 
##' ## statified correlation (using another dataset)
##' data(gastricbypassW, package = "LMMstar")
##' gastricbypassW$weight.bin <- gastricbypassW$weight1>=120
##' partialCor(glucagonAUC1+glucagonAUC2~1, data = gastricbypassW, by = "weight.bin")
##'
##' ## compared correlation between groups
##' partialCor(glucagonAUC1+glucagonAUC2~1, data = gastricbypassW, by = "weight.bin",
##'            effects = "Dunnett") 
##' 
##' #### with repetitions ####
##' \dontrun{
##' data(gastricbypassL, package = "LMMstar")
##' ## via a mixed model
##' eUN.lmm <- lmm(weight ~ glucagonAUC+time, repetition =~time|id,
##'                data = gastricbypassL, structure = "UN")
##' partialCor(eUN.lmm)
##' 
##' ## mean: variable and timepoint specific mean parameter (8)
##' ## variance: variable and timepoint specific variance parameter (8)
##' ## correlation: correlation parameter specific for each variable and time lag (10)
##' e.cor <- partialCor(weight+glucagonAUC~time, repetition =~time|id,
##'                     data = gastricbypassL, structure = "LAG")
##' e.cor
##' coef(attr(e.cor,"lmm"), effects = "correlation")
##' if(require(ggplot2)){
##' autoplot(e.cor)
##' }
##'
##' ## same except for the mean structure: variable specific mean parameter (2)
##' e.cor2 <- partialCor(weight+glucagonAUC~time, repetition =~time|id,
##'                     data = gastricbypassL, structure = "LAG")
##'
##' ## mean: variable and timepoint specific mean parameter (8)
##' ## variance: variable specific variance parameter (2)
##' ## correlation: correlation parameter specific for each variable and some time lag (4)
##' e.cor3 <- partialCor(weight+glucagonAUC~time, repetition =~time|id,
##'                      data = gastricbypassL, structure = "CS")
##' e.cor3
##' coef(attr(e.cor3,"lmm"), effects = "correlation")
##' if(require(ggplot2)){
##' autoplot(e.cor3)
##' }
##' 
##' }
#' @export
`partialCor` <-
  function(object, ...) UseMethod("partialCor")

## * partialCor.list (code)
##' @rdname partialCor
##' @export
partialCor.list <- function(object, data, repetition = NULL, structure = NULL, by = NULL,
                            effects = NULL, rhs = NULL, method = "none", df = NULL, transform.rho = NULL, name.short = c(TRUE,FALSE), ...){

    ## ... confidence level
    
    ## ** normalize arguments
    data <- as.data.frame(data)

    ## *** check object
    if(any(unlist(lapply(object, inherits, "formula"))==FALSE)){
        stop("Argument \'object\' should be a formula or list of formula. \n")
    }
    if(!missing(repetition) && is.null(repetition) && length(object)<2){
        stop("Argument \'object\' should contain at least two formula. \n")
    }else if(!missing(repetition) && !is.null(repetition) && length(object)!=2){
        stop("Argument \'object\' should contain exactly two formula. \n")
    }
    
    ## *** check formula agree with data
    ls.name.XY <- lapply(object,all.vars)
    name.XY <- unique(unlist(ls.name.XY))
    if(any(name.XY %in% names(data) == FALSE)){
        stop("Variable(s) \"", paste(name.XY[name.XY %in% names(data) == FALSE], collapse = "\" \""),"\" are not in argument \'data\'. \n")
    }

    if(any(c("CCvariableCC","CCvalueCC","CCrepetitionCC") %in% names(data))){
        stop("Argument \'data\' should not contain columns \"CCvariableCC\", \"CCvalueCC\", or \"CCrepetitionCC\". \n",
             "Those names are used internally. \n")
    }

    ## *** id and time
    if(!missing(repetition)){
        if(!inherits(repetition,"formula")){
            stop("Argument \'repetition\' must be of class formula, something like: ~ time|cluster or strata ~ time|cluster. \n")
        }

        repetition.split <- strsplit(x = deparse(repetition), split = "|",fixed=TRUE)[[1]]
        if(length(repetition.split)!=2){
            stop("Incorrect specification of argument \'repetition\'. \n",
                 "The symbol | should only exacly once, something like: ~ time|cluster or strata ~ time|cluster. \n")
        }
        name.id <- trimws(repetition.split[2], which = "both")
        if(any(name.id %in% names(data) == FALSE)){
            stop("Cluster variable \"", name.id,"\" is not in argument \'data\'. \n")
        }
        name.time <- all.vars(stats::as.formula(repetition.split[1]))
        if(length(name.time)==0){
            stop("Missing \'time\' variable in the repetition argument. \n",
                 "Should be something of the form ~time|id. \n")
        }
        if(any(name.id %in% names(data) == FALSE)){
            stop("Repetition variable \"", name.time,"\" is not in argument \'data\'. \n")
        }
        
    }else{
        if("CCindexCC" %in% names(data)){
            stop("Argument \'data\' should not contain a variable \"CCindexCC\". \n")
        }
        data$CCindexCC <- 1:NROW(data)
        name.id <- "CCindexCC"
        data$CCrepetitionCC <- 1
        name.time <- "CCrepetitionCC"
    }

    ## *** by
    if(!is.null(by)){
        if(length(by)!=1){
            stop("Argument \'by\' must have length 1. \n")
        }
        if(by %in% names(data) == FALSE){
            stop("Argument \'by\' should correspond to a column name of argument \'data\'. \n")
        }
    }

    ## *** contrast
    valid.contrMat <- c("Dunnett", "Tukey", "Sequen")
    if(length(effects)==1 && (effects %in% valid.contrMat == FALSE)){
        stop("When character, argument \'effects\' should be one of \"",paste(valid.contrMat,collapse="\" \""),"\".\n")
    }

    ## *** df
    options <- LMMstar.options()
    if(is.null(df)){
        df <- options$df
    }

    ## ** reshape    
    ls.name.X <- lapply(object, function(iF){all.vars(stats::delete.response(stats::terms(iF)))})
    name.X <- unique(unlist(ls.name.X))
    ls.name.Y <- lapply(ls.name.XY, function(iF){setdiff(iF,c(name.X,name.id))})
    name.Y <- unlist(ls.name.Y)
    if(any(duplicated(name.Y))){
        stop("Variables in the left hand side of argument should be unique. \n")
    }
    dataL <- stats::reshape(data[, unique(c(name.XY, name.id, name.time, by)),drop=FALSE], direction  = "long",
                            idvar = c(name.id, name.time, by),
                            varying = name.Y,
                            v.names = "CCvalueCC",
                            timevar = "CCvariableCC")
    dataL$CCvariableCC <- factor(dataL$CCvariableCC, labels = name.Y)
    rownames(dataL) <- NULL
    
    ## ** prepare for mixed model
    index.interaction <- which(colSums(1-do.call(rbind,lapply(ls.name.X, function(iX){name.X %in% iX})))==0)

    X.formula <- name.X
    if(length(index.interaction)>0){
        X.formula[index.interaction] <- paste0(name.X[index.interaction],":","CCvariableCC")
    }

    ## set value to reference when not in the ajustment set    
    for(iY in 1:length(name.Y)){ ## iY <- 1
        for(iX in name.X){ ## iX <- name.X[1]
            if(iX %in% ls.name.X[[iY]]==FALSE){
                if(is.numeric(data[[iX]])){
                    dataL[dataL$CCvariableCC == name.Y[iY],iX]  <- 0
                }else if(is.factor(data[[iX]])){
                    dataL[dataL$CCvariableCC == name.Y[iY],iX]  <- levels(data[[iX]])[1]
                }else if(is.character(data[[iX]])){
                    data[[iX]] <- as.factor(data[[iX]])
                    dataL[dataL$CCvariableCC == name.Y[iY],iX]  <- levels(data[[iX]])[1]
                }else{
                    stop("Unknown type of variable for ",iX,". \n")
                }
            }            
        }
    }

    test.duplicated <- tapply(dataL[[name.id]], dataL[["CCvariableCC"]], function(iId){any(duplicated(iId))})
    if(length(X.formula)>0){
        formula.mean <- stats::as.formula(paste("CCvalueCC~CCvariableCC+",paste(X.formula, collapse = "+")))
    }else{
        formula.mean <- CCvalueCC~CCvariableCC
    }

    ## ** fit mixed model
    if(any(test.duplicated)){

        ## *** TOEPLITZ mixed model (repetition)
        dataL <- dataL[order(dataL[[name.id]],dataL$CCvariableCC),]
        dataL$CCrepetitionCC <- unlist(tapply(dataL[[name.id]],dataL[[name.id]],function(iId){1:length(iId)}))
        formula.repetition <- stats::as.formula(paste("~ CCvariableCC + ",name.time,"|",name.id))
        if(is.null(structure)){
            structure <- "CS"
        }else{
            structure <- match.arg(structure, c("UN","PEARSON","HLAG","LAG","HCS","CS"))
        }


        if(structure == "UN"){
            structure2 <- do.call(TOEPLITZ, args = list(formula = list(stats::as.formula(paste("~CCvariableCC+",name.time)),
                                                                       stats::as.formula(paste("~CCvariableCC+",name.time))),
                                                        type = "UN", add.time = FALSE))
        }else if(structure == "PEARSON"){
            structure2 <- do.call(TOEPLITZ, args = list(formula = list(~CCvariableCC,
                                                                       stats::as.formula(paste("~CCvariableCC+",name.time))),
                                                        type = "UN", add.time = FALSE))
        }else if(structure=="HLAG"){
            structure2 <- do.call(TOEPLITZ, args = list(formula = list(stats::as.formula(paste("~CCvariableCC+",name.time)),
                                                                       stats::as.formula(paste("~CCvariableCC+",name.time))),
                                                        type = "LAG", add.time = FALSE))
        }else if(structure=="LAG"){
            structure2 <- do.call(TOEPLITZ, args = list(formula = list(~CCvariableCC,
                                                                       stats::as.formula(paste("~CCvariableCC+",name.time))),
                                                        type = "LAG", add.time = FALSE))
        }else if(structure=="HCS"){
            structure2 <- do.call(TOEPLITZ, args = list(formula = list(stats::as.formula(paste("~CCvariableCC+",name.time)),
                                                                       stats::as.formula(paste("~CCvariableCC+",name.time))),
                                                        type = "CS", add.time = FALSE))
        }else if(structure=="CS"){
            structure2 <- do.call(TOEPLITZ, args = list(formula = list(~CCvariableCC,
                                                                       stats::as.formula(paste("~CCvariableCC+",name.time))),
                                                        type = "CS", add.time = FALSE))
        }
        if(is.null(by)){
            e.lmm <- lmm(formula.mean, df = df, repetition = formula.repetition,
                         data = dataL, structure = structure2,
                         control = list(optimizer = "FS"))
            out.confint <- confint(e.lmm, df = df, columns = c("estimate","se","df","lower","upper","p.value"), effects = "correlation", transform.rho = transform.rho, ...)

            ## identify the marginal correlation coefficient
            name.rho <- e.lmm$design$param[e.lmm$design$param$type=="rho","name"]
            name.rho.marginal <- grep("dt=0", name.rho, fixed=TRUE, value = TRUE)
            out <- out.confint[name.rho.marginal,,drop=FALSE]
            rownames(out) <- "marginal"

            ## compute conditional correlation
            if(structure %in% c("CS","HCS")){
                code.rho <- e.lmm$design$param[e.lmm$design$param$type=="rho","code"]
                name.rhoWithinBlock <- name.rho[grepl("R.",code.rho)]
                name.rhoAcrossBlock <- setdiff(name.rho, c(name.rhoWithinBlock, name.rho.marginal))
                test.atanh <- identical(attr(out,"backtransform")$FUN,"tanh")
                out2 <- estimate(e.lmm, df = df, f = function(p){ ## p <- coef(e.lmm, effects = "all")
                    if(any(p[name.rhoWithinBlock]<0)){
                        iOut <- c(NA,NA)
                    }else{
                        iOut <- c((p[name.rho.marginal]-p[name.rhoAcrossBlock])/sqrt(prod(1-p[name.rhoWithinBlock])),
                                  p[name.rhoAcrossBlock]/sqrt(prod(p[name.rhoWithinBlock])))
                        if(test.atanh){iOut <- atanh(iOut)}
                    }
                    return(iOut)                        
                }, ...)

                if(test.atanh){
                    out2 <- .backtransform(out2, type.param = rep("rho",NROW(out2)), backtransform.names = NULL, backtransform = c(FALSE,FALSE,FALSE,TRUE),
                                           transform.mu = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = "atanh")
                        
                }
                rownames(out2) <- c("conditional","latent")
                out <- rbind(out, out2)
            }            
            attr(out,"name.rho") <- name.rho.marginal
            
        }else{
            e.lmm <- mlmm(formula.mean, df = df, repetition = formula.repetition, data = dataL, structure = structure2, transform.rho = transform.rho, 
                          by = by, effects = "correlation", contrast.rbind = effects, control = list(optimizer = "FS"), trace = FALSE)
            out <- confint(e.lmm, df = df, columns = c("estimate","se","df","lower","upper","p.value"), ...)

        }

        
    }else{
        ## *** UN/CS mixed model (no repetition)
        formula.repetition <- stats::as.formula(paste("~CCvariableCC|",name.id))
        if(is.null(structure)){
            structure <- "UN"
        }else{
            structure <- match.arg(structure, c("UN","CS"))
        }

        if(is.null(by)){
            ## fit model
            e.lmm <- lmm(formula.mean, repetition = formula.repetition,
                         data = dataL, structure = structure, df = df)

            name.param <- e.lmm$design$param$name
            n.param <- length(name.param)
            name.cor <- e.lmm$design$param[e.lmm$design$param$type=="rho","name"]
                
            ## define contrast
            if(length(effects)==0){
                Cmat <- matrix(0, nrow = length(name.cor), ncol = n.param,
                               dimnames = list(name.cor,name.param))
                if(length(name.cor)==1){
                    Cmat[name.cor,name.cor] <- 1
                }else{
                    diag(Cmat[name.cor,name.cor]) <- 1
                }
                
            }else if(length(effects)==1 && is.character(effects)){
                contrast <- multcomp::contrMat(rep(1,length(name.cor)), type = effects)
                Cmat <- matrix(0, nrow = NROW(contrast), ncol = n.param,
                               dimnames = list(NULL,name.param))
                Cmat[,name.cor] <- unname(contrast)
                rownames(Cmat) <- unlist(lapply(strsplit(split = "-",rownames(contrast),fixed=TRUE), function(iVec){
                    paste(name.cor[as.numeric(trimws(iVec))], collapse = " - ")
                }))
            }

            ## run test linear hypothesis
            if(all(rowSums(Cmat!=0)==1) || !is.null(transform.rho)){
                out <- confint(anova(e.lmm, effects = Cmat, transform.rho = transform.rho),
                               method = method, columns = c("estimate","se","df","lower","upper","p.value"), ...)
            }else{
                out0 <- confint(anova(e.lmm, effects = Cmat, transform.rho = "none"),
                                method = "none", columns = c("estimate","se","df","p.value"), ...)
                out <- confint(anova(e.lmm, effects = Cmat, transform.rho = "atanh"),
                               method = method, columns = c("estimate","se","df","p.value"), ...)
                out$estimate <- out0$estimate
                out$se <- out0$se
                out$df <- out0$df
                attr(out,"backtransform")[,"estimate"] <- FALSE
            }
        }else{
            e.lmm <- mlmm(formula.mean, df = df, repetition = formula.repetition, data = dataL, structure = structure, transform.rho = transform.rho,
                          by = by, effects = "correlation", contrast.rbind = effects, name.short = name.short, trace = FALSE)
            out <- confint(e.lmm, columns = c("estimate","se","df","lower","upper","p.value"))
        }
    }

    ## ** export
    attr(out, "call") <- match.call()
    attr(out, "args") <- list(by = by, effects = effects, method = method, df = df, transform.rho = transform.rho, level = attr(out,"level"))
    attr(out, "level") <- NULL
    attr(out, "lmm") <- e.lmm
    class(out) <- append("partialCor", class(out))
    return(out)
}

## * partialCor.formula (code)
##' @rdname partialCor
##' @export
partialCor.formula <- function(object, repetition, ...){

    formula.rhs <- stats::delete.response(stats::terms(object))
    response <- setdiff(all.vars(object),all.vars(formula.rhs))
    if(!missing(repetition) && is.null(repetition) && length(response)<2){
        stop("Argument \'object\' should contain at least two variables on the left hand side of the formula. \n")
    }else if(!missing(repetition) && !is.null(repetition) && length(response)!=2){
        stop("Argument \'object\' should contain exactly two variables on the left hand side of the formula. \n")
    }
    ls.object <- lapply(response, function(iY){stats::as.formula(paste(iY,deparse(formula.rhs)))}) ## iY <- response[1]        
    
    return(partialCor.list(ls.object, repetition = repetition, ...))
}

## * partialCor.lmm (code)
##' @rdname partialCor
##' @export
partialCor.lmm <- function(object, level = 0.95, R2 = FALSE, se = TRUE, df = TRUE, ...){

    ## ** normalize input
    object.df <- object$df
    if(is.null(object.df)){
        stop("Cannot compute the partial correlation when the degrees of freedom are not computed. \n",
             "Consider calling \'lmm\' with the argument df=TRUE. \n")
    }
    mytable <- model.tables(object, columns = c("estimate","statistic","df"))

    name.param <- setdiff(rownames(mytable),"(Intercept)")
    n.param <- length(name.param)

    out <- matrix(NA, nrow = n.param, ncol = 6,
                  dimnames = list(name.param, c("estimate","se","df","lower","upper","p.value")))

    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** partial correlation
    value.meancoef <- coef(object, effects = "mean")
    name.meancoef <- names(value.meancoef)
    keep.meancoef <- setdiff(name.meancoef, "(Intercept)")

    if(se == FALSE){
        SSEstar <- stats::df.residual(object)
        Einfo <- vcov(object, type.information = "expected", effects = "mean")
        SSRstar.indiv <- value.meancoef[keep.meancoef]^2 / diag(Einfo)[keep.meancoef]
        out[,"estimate"] <- sign(value.meancoef[keep.meancoef])*sqrt(abs(SSRstar.indiv) / (SSRstar.indiv + SSEstar))
    }else{
        out <- lava::estimate(object, df = df, level = level, function(p){ ## p <- coef(object, effects = "all")

            iSSEstar <- stats::df.residual(object, p = p)
            iEinfo <- vcov(object, p = p, type.information = "expected", effects = "mean")
            iSSRstar.indiv <- p[keep.meancoef]^2 / diag(iEinfo)[keep.meancoef]
            return(sign(p[keep.meancoef])*sqrt(abs(iSSRstar.indiv) / (iSSRstar.indiv + iSSEstar)))

        })
    }

    ## via another formula
    ## out[,"estimate"] <- sign(mytable[name.param,"statistic"])*sqrt(mytable[name.param,"statistic"]^2/(mytable[name.param,"df"]+mytable[name.param,"statistic"]^2))
    ## out <- estimate(object, df = df, level = level, function(p){ ## p <- coef(object, effects = "all")
    ##     newSigma <- vcov(object, p = p, df = TRUE)
    ##     newStat <- p[name.param]/sqrt(diag(newSigma[name.param,name.param,drop=FALSE]))
    ##     newDf <- attr(newSigma,"df")[name.param]
    ##     return(sign(newStat)*sqrt(newStat^2/(newDf + newStat^2)))
    ## })

    ## ** R2
    if(R2){
        ## linear contrasts
        ls.C <- lapply(anova(object, effects = "mean", df = FALSE, ci = TRUE)$glht[[1]], function(iAOV){
            iAOV$linfct[,name.meancoef,drop=FALSE]
        })
        ls.C[[length(ls.C)+1]] <- matrix(0, nrow = length(keep.meancoef), ncol = NCOL(ls.C[[1]]),
                                         dimnames = list(keep.meancoef, colnames(ls.C[[1]])))
        if(length(keep.meancoef)==1){
            ls.C[[length(ls.C)]][keep.meancoef,keep.meancoef] <- 1
        }else{
            diag(ls.C[[length(ls.C)]][keep.meancoef,keep.meancoef]) <- 1
        }
        names(ls.C)[length(ls.C)] <- "global"

        out2 <- matrix(NA, nrow = length(ls.C), ncol = 6,
                       dimnames = list(names(ls.C), c("estimate","se","df","lower","upper","p.value")))

        if(se == FALSE){

            SSEstar <- stats::df.residual(object)
            Einfo <- vcov(object, type.information = "expected", effects = "mean")
            SSRstar.indiv <- sapply(ls.C, function(iC){t(iC %*% cbind(value.meancoef)) %*% solve(iC %*% Einfo %*% t(iC)) %*% (iC %*% cbind(value.meancoef))})
            out2[,"estimate"] <- SSRstar.indiv / (SSRstar.indiv + SSEstar)

        }else{
            out2 <- estimate(object, df = df, level = level, function(p){ ## p <- coef(object, effects = "all")

                iSSEstar <- stats::df.residual(object, p = p)
                iEinfo <- vcov(object, p = p, type.information = "expected", effects = "mean")
                iSSRstar.indiv <- sapply(ls.C, function(iC){t(iC %*% cbind(p[name.meancoef])) %*% solve(iC %*% iEinfo %*% t(iC)) %*% (iC %*% cbind(p[name.meancoef]))})
                return(iSSRstar.indiv / (iSSRstar.indiv + iSSEstar))

            })
        }

    }else{
        out2 <- NULL
    }
    
    ## ** warning
    if(sum(object$design$param$type=="sigma")>1){
        warning("Formula for the partial correlation",if(R2){" and R2"}," may not lead to meaningful estimates with heteroschedastic residuals. \n",
                "Use the estimated values at your own risk. \n", sep = "")
    }
    
    ## ** export
    attr(out, "call") <- match.call()
    attr(out, "args") <- list(level = level, R2 = R2, se = se, df = df)
    attr(out, "R2") <- out2
    class(out) <- append("partialCor", class(out))
    return(out)
}

## * R2 stuff
## ## ** partial percentage of variance explained (R2)
## index.cluster <- object$design$index.cluster
## index.clusterTime <- attr(object$design$index.clusterTime,"vectorwise")
## X.pattern <- object$design$vcov$Upattern
## name.pattern <- X.pattern$name

## ## design matrix
## X.design <- object$design$mean[,name.param,drop=FALSE]
## assignX.design <- attr(object$design$mean,"assign")[match(name.param,colnames(object$design$mean))]
    
## ## linear predictor per group of variables
## Xbeta.param <- sweep(X.design, FUN = "*", MARGIN = 2, STATS = mytable[name.param,"estimate"])
## Xbeta.assign <- do.call(cbind,lapply(unique(assignX.design), function(iAssign){ ## iAssign <- 1
##     rowSums(Xbeta.param[,name.param[iAssign],drop=FALSE])
## }))
## Xbeta2.assign <- Xbeta.assign^2
    
## ## variance pattern
## Xpattern.var <- object$design$vcov$var$Xpattern
## nXpattern.var <- length(Xpattern.var)
## nameXpattern.var <- names(Xpattern.var)
## ls.nameXpattern.var <- strsplit(nameXpattern.var, split = ".", fixed = TRUE)
    
## M.timevar <- list(pattern = matrix(NA, nrow = nXpattern.var, ncol = object$time$n),
##                   varRes = matrix(NA, nrow = nXpattern.var, ncol = object$time$n),
##                   n = matrix(NA, nrow = nXpattern.var, ncol = object$time$n),
##                   sumExp = array(NA, dim = c(nXpattern.var, object$time$n, NCOL(Xbeta.assign))),
##                   sumExp2 = array(NA, dim = c(nXpattern.var, object$time$n, NCOL(Xbeta.assign)))
##                   )
    
## for(iPattern in 1:nXpattern.var){
##     iIndex.pattern <- which(object$design$vcov$Upattern$var == nameXpattern.var[iPattern])[1]
##     iIndex.time <- object$design$vcov$Upattern$time[[iIndex.pattern]]
##     iIndex.cluster <- object$design$vcov$Upattern$index.cluster[[iIndex.pattern]]
##     iMindex.cluster <- do.call(rbind,index.cluster[iIndex.cluster])

##     M.timevar$pattern[iPattern,iIndex.time] <- ls.nameXpattern.var[[iPattern]]
##     M.timevar$varRes[iPattern,iIndex.time] <- diag(object$Omega[[iIndex.pattern]])
##     M.timevar$n[iPattern,iIndex.time] <- object$design$vcov$X$Upattern$n.cluster[iIndex.pattern]
##     M.timevar$sumExp[iPattern,iIndex.time,] <- apply(iMindex.cluster, 2, function(iCluster.time){
##         colSums(Xbeta.assign[iCluster.time,,drop=FALSE])
##     })
##     M.timevar$sumExp2[iPattern,iIndex.time,] <- apply(iMindex.cluster, 2, function(iCluster.time){
##         colSums(Xbeta2.assign[iCluster.time,,drop=FALSE])
##     })
## }
## ## marginal
## level.pattern.var <- unique(M.timevar$pattern[!is.na(M.timevar$pattern)])
## M.R2marginal <- do.call(rbind,lapply(level.pattern.var, function(iLevel){ ## iLevel <- "1"

    ##     ## set.seed(10)
    ##     ## X <- rnorm(100, mean = 5)
    ##     ## var(X)
    ##     ## sum((X-mean(X))^2)/(length(X)-1)
    ##     ## sum(X^2)/(length(X)-1) - 2*length(X)*mean(X)^2/(length(X)-1) + length(X)*mean(X)^2/(length(X)-1)
    ##     ## sum(X^2)/(length(X)-1) - sum(X)^2/(length(X)*(length(X)-1))

    ##     iMask <- which(M.timevar$pattern==iLevel)
    ##     iN <- sum(M.timevar$n[iMask])
    ##     iSumExp <- apply(M.timevar$sumExp,3,function(iM){sum(iM[iMask])})
    ##     iSumExp2 <- apply(M.timevar$sumExp2,3,function(iM){sum(iM[iMask])})
    ##     return(c(iN,iSumExp2/(iN-1) - iSumExp^2/(iN^2-iN)))
    ## }))

    ## R2weight <- M.R2marginal[,1]/sum(M.R2marginal[,1])
    ## out[[2]][,"estimate"] <-  apply(M.R2marginal[,-1],2,weighted.mean, w = R2weight)

##' \bold{Explained variance and partial correlation}: can be extracted by adding \code{"partial.r"} to the argument \code{columns} for certain types of mixed models,
##' those equivalent to random effect models (see vignette). 
##' Explained variance will be displayed in the column \code{partial.r2} (Multivariate section)
##' and partial correlation will be displayed in the column \code{partial.r} (Univariate section).
##' WARNING: for other type of mixed models, typically heteroschedastic, the values in the columns \code{partial.r2} and \code{partial.r} may not have any intuitive interpretation. 
##             ## *** R2
##             X.design <- object$design$mean
            
##             if(df>0 && all(colnames(outSimp$C)[which(outSimp$C!=0, arr.ind = TRUE)[,"col"]] %in% colnames(X.design))){
                

##                 Mpartial.R2 <- do.call(rbind,lapply(name.pattern, function(iPattern){ ## iPattern <- name.pattern[1]
##                     iIndex.cluster <- X.pattern$index.cluster[[iPattern]]
##                     iMindex.cluster <- do.call(rbind,index.cluster[iIndex.cluster])

##                     iM.Num <- NA*iMindex.cluster
##                     iM.Num[] <- Xbeta.design[iMindex.cluster]
##                     iNum <- apply(iM.Num,2,var)
                    
##                     return(iNum/(iNum+diag(object$Omega[[iPattern]])))
##                 }))

##                 if(object$design$vcov$type == "ID"){

##                     ## R2 is computed as the averaged percentage over variance explained
##                     ## i.e. variance explained in each variance strata, weighted average as a function of the sample size
##                     paramC.design <- (param*colSums(outSimp$C))[colnames(X.design)]
##                     Mpartial.R2 <- do.call(rbind,lapply(name.pattern, function(iPattern){ ## iPattern <- name.pattern[1]
##                         iIndex.cluster <- attr(X.pattern[[iPattern]],"index.cluster")
##                         iVarBeta <- var(rowSums(sweep(X.design[iIndex.cluster,], FUN = "*", MARGIN = 2, STATS = paramC.design)))
##                         iSigma <- as.double(object$Omega[[iPattern]])
##                         return(c(R2 = iVarBeta/(iVarBeta+iSigma),n = length(iIndex.cluster)))
##                     }))
##                     partial.r2 <- sum(Mpartial.R2[,"R2"] * (Mpartial.R2[,"n"] / sum(Mpartial.R2[,"n"])))

##                 }else if(object$design$vcov$type == "CS" && object$design$vcov$heterogeneous == FALSE){
                    
##                   .nestingRanef(object)
  
##                     partial.r2 <- sum(Mpartial.R2[,"R2"] * (Mpartial.R2[,"n"] / sum(Mpartial.R2[,"n"])))

##                 }else{
##                     partial.r2 <- NA
##                 }
##             }else{
##                 partial.r2 <- NA
##             }

##----------------------------------------------------------------------
### partialCor.R ends here
