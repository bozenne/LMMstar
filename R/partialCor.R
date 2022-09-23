### partialCor.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May  1 2022 (17:01) 
## Version: 
## Last-Updated: sep 21 2022 (11:55) 
##           By: Brice Ozenne
##     Update #: 261
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
##' @description Estimate the partial correlation between two variables where the adjustment set may differ between variables.
##' 
##' @param formula a formula with in the left hand side the variables for which the correlation should be computed
##' and on the right hand side the adjustment set. Can also be a list of formula for outcome-specific adjustment set.
##' @param data [data.frame] dataset containing the variables.
##' @param repetition [formula] Specify the structure of the data: the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' @param heterogeneous [logical] Specify whether the variance differ between the two variables (no repetition)
##' or whether the correlation/variance vary over repetitions (repetition).
##' @param by [character] variable used to stratified the correlation on.
##' @param effects [character or matrix] type of contrast to be used for comparing the correlation parameters. One of \code{"Dunnett"}, \code{"Tukey"}, \code{"Sequen"}, or a contrast matrix.
##' @param rhs [numeric vector] right hand side for the comparison of correlation parameters. 
##' @param method [character] adjustment for multiple comparisons (e.g. \code{"single-step"}).
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the coefficient. Otherwise a normal distribution is used.
##' @param transform.rho [character] scale on which perform statistical inference (e.g. \code{"atanh"})
##'
##' @details Fit a mixed model to estimate the partial correlation with the following variance-covariance pattern:
##' \itemize{
##' \item \bold{no repetition}: unstructure or compound symmetry structure for M observations, M being the number of variable on the left hand side (i.e. outcomes).
##' \item \bold{repetition}: toeplitz for M*T observations where T denotes the number of repetitions. In that case heterogeneous can be 1, 0.5, or 0.
##' }
##'
##' @return A data.frame with the estimate partial correlation (rho), standard error, degree of freedom, confidence interval, and p-value (test of no correlation).
##' When \code{heterogeneous=FALSE} is used with repeated measurements, a second correlation coefficient (r) is output where the between subject variance has been removed (similar to Bland et al. 1995).
##'
##' @references
##'  Bland J M, Altman D G. Statistics notes: Calculating correlation coefficients with repeated observations: Part 1—correlation within subjects BMJ 1995; 310 :446 doi:10.1136/bmj.310.6977.446 
##' 
##' @examples
##' #### bivariate (no repetition) ####
##' ## example from ppcor::pcor 
##' y.data <- data.frame(
##'   hl=c(7,15,19,15,21,22,57,15,20,18),
##'   disp=c(0.000,0.964,0.000,0.000,0.921,0.000,0.000,1.006,0.000,1.011),
##'   deg=c(9,2,3,4,1,3,1,3,6,1),
##'   BC=c(1.78e-02,1.05e-06,1.37e-05,7.18e-03,0.00e+00,0.00e+00,0.00e+00
##',  4.48e-03,2.10e-06,0.00e+00)
##')
##' 
##' ## ppcor::pcor(y.data)
##' ## estimate
##' ##              hl       disp        deg        BC
##' ## hl    1.0000000 -0.6720863 -0.6161163 0.1148459
##' ## disp -0.6720863  1.0000000 -0.7215522 0.2855420
##' ## deg  -0.6161163 -0.7215522  1.0000000 0.6940953
##' ## BC    0.1148459  0.2855420  0.6940953 1.0000000
##' 
##' ## $p.value
##' ##              hl       disp        deg         BC
##' ## hl   0.00000000 0.06789202 0.10383620 0.78654997
##' ## disp 0.06789202 0.00000000 0.04332869 0.49299871
##' ## deg  0.10383620 0.04332869 0.00000000 0.05615021
##' ## BC   0.78654997 0.49299871 0.05615021 0.00000000
##'
##' set.seed(10)
##' y.data$gender <- factor(rbinom(10, size = 1, prob = 0.5), labels = c("F","M"))
##' 
##' partialCor(c(hl,disp)~BC+deg, data = y.data)
##' partialCor(hl + disp~BC+deg, data = y.data)
##'
##' partialCor(list(hl~BC+deg, disp~BC), data = y.data)
##' partialCor(list(hl~BC+deg+gender, disp~1), data = y.data)
##' 
##' #### stratified bivariate (no repetition) ####
##' ## partialCor(list(hl~1, disp~1), data = y.data, by = "gender") ## too small dataset
##' 
##' #### bivariate (with repetition) ####
##' \dontrun{
##' data(gastricbypassL)
##' ## mean: variable and timepoint specific mean parameter (8)
##' ## variance: variable and timepoint specific variance parameter (8)
##' ## correlation: correlation parameter specific for each variable and time lag (10)
##' e.cor <- partialCor(weight+glucagonAUC~time, repetition =~time|id,
##'                     data = gastricbypassL)
##' e.cor
##' coef(attr(e.cor,"lmm"), effects = "correlation")
##' if(require(ggplot2)){
##' autoplot(e.cor)
##' }
##'
##' ## same except for the mean structure: variable specific mean parameter (2)
##' e.cor2 <- partialCor(weight+glucagonAUC~time, repetition =~time|id,
##'                     data = gastricbypassL)
##'
##' ## mean: variable and timepoint specific mean parameter (8)
##' ## variance: variable specific variance parameter (2)
##' ## correlation: correlation parameter specific for each variable and some time lag (4)
##' e.cor3 <- partialCor(weight+glucagonAUC~time, repetition =~time|id,
##'                      data = gastricbypassL, heterogeneous = 0.5)
##' e.cor3
##' coef(attr(e.cor3,"lmm"), effects = "correlation")
##' if(require(ggplot2)){
##' autoplot(e.cor3)
##' }
##' 
##' ## mean: variable and timepoint specific mean parameter (8)
##' ## variance: variable specific variance parameter (2)
##' ## correlation: correlation parameter specific for each variable and some time lag (4)
##' e.cor4 <- partialCor(weight+glucagonAUC~time, repetition =~time|id,
##'                      data = gastricbypassL, heterogeneous = 0)
##' e.cor4
##' coef(attr(e.cor3,"lmm"), effects = "correlation")
##' if(require(ggplot2)){
##' autoplot(e.cor3)
##' }##' }

## * partialCor (documentation)
##' @export
partialCor <- function(formula, data, repetition = NULL, heterogeneous = TRUE, by = NULL,
                       effects = NULL, rhs = NULL, method = "none", df = NULL, transform.rho = NULL){

    ## ** normalize arguments
    data <- as.data.frame(data)

    ## *** convert first argument to list of formula
    if(inherits(formula,"formula")){
        formula.rhs <- stats::delete.response(stats::terms(formula))
        response <- setdiff(all.vars(formula),all.vars(formula.rhs))
        if(is.null(repetition) && length(response)<2){
            stop("Argument \'formula\' should contain at least two variables on the left hand side of the formula. \n")
        }else if(!is.null(repetition) && length(response)!=2){
            stop("Argument \'formula\' should contain exactly two variables on the left hand side of the formula. \n")
        }
        formula <- lapply(response, function(iY){stats::as.formula(paste(iY,deparse(formula.rhs)))}) ## iY <- response[1]        
    }else{
        if(!is.list(formula)){
            stop("Argument \'formula\' should be a formula or a list of formula. \n")
        }
        if(any(unlist(lapply(formula, inherits, "formula"))==FALSE)){
            stop("Argument \'formula\' should be a formula or list of formula. \n")
        }
        if(is.null(repetition) && length(formula)<2){
            stop("Argument \'formula\' should contain at least two formula. \n")
        }else if(!is.null(repetition) && length(formula)!=2){
            stop("Argument \'formula\' should contain exactly two formula. \n")
        }
    }


    ## *** check formula agree with data
    ls.name.XY <- lapply(formula,all.vars)
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
    ls.name.X <- lapply(formula, function(iF){all.vars(stats::delete.response(stats::terms(iF)))})
    name.X <- unique(unlist(ls.name.X))

    ls.name.Y <- lapply(ls.name.XY, function(iF){setdiff(iF,c(name.X,name.id))})
    name.Y <- unlist(ls.name.Y)
    if(any(duplicated(name.Y))){
        stop("Variables in the left hand side of argument should be unique. \n")
    }
    dataL <- stats::reshape(data[, unique(c(name.XY, name.id, name.time,by)),drop=FALSE], direction  = "long",
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
        formula.repetition <- stats::as.formula(paste("~",name.time,"+CCvariableCC|",name.id))
        if(heterogeneous>=1){
            structure <- TOEPLITZ(heterogeneous = TRUE)
        }else if(heterogeneous>=0.5){
            structure <- TOEPLITZ(heterogeneous = 0.5)
        }else{
            structure <- TOEPLITZ(heterogeneous = FALSE)
        }

        if(is.null(by)){
            browser()
            e.lmm <- lmm(formula.mean, df = df, repetition = formula.repetition,
                         data = dataL, structure = structure,
                         control = list(optimizer = "FS"))
            out <- confint(e.lmm, df = df, columns = c("estimate","se","df","lower","upper","p.value"), effects = "correlation", transform.rho = transform.rho)

            ## identify the right correlation coefficient
            code.rho <- e.lmm$design$param[e.lmm$design$param$type=="rho","code"]
            name.rho <- e.lmm$design$param[e.lmm$design$param$type=="rho","name"][grepl("D.",code.rho)]

            M.time <- attr(e.lmm$time$levels,"original")
            U.time <- unique(M.time[,1])
            tentative.rho <- sapply(U.time, function(iT){paste0("rho(",paste(interaction(M.time)[M.time[,1]==iT],collapse=","),")")})

            if(any(tentative.rho %in% rownames(out))){
                keep.rho <- intersect(tentative.rho,rownames(out)) 
                out <- out[keep.rho,,drop=FALSE]

                ## compute conditional correlation
                if((length(keep.rho)==1) && (heterogeneous<1)){
                    name.rho2 <- e.lmm$design$param[e.lmm$design$param$type=="rho","name"][grepl("R.",code.rho)]
                    sub.rho <- setdiff(name.rho, keep.rho)

                    test.atanh <- identical(attr(out,"backtransform")$FUN,"tanh")

                    out2 <- estimate(e.lmm, df = df, f = function(p){
                        if(any(p[name.rho2]<0)){
                            iOut <- NA
                        }else{
                            iOut <- (p[keep.rho]-p[sub.rho])/sqrt(prod(1-p[name.rho2]))
                            if(test.atanh){iOut <- atanh(iOut)}
                        }
                        return(iOut)                        
                    })

                    if(test.atanh){
                        out2 <- .backtransform(out2, type.param = "rho", backtransform.names = NULL, backtransform = c(FALSE,FALSE,FALSE,TRUE),
                                               transform.mu = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = "atanh")
                        
                    }

                    rownames(out2) <- gsub("^rho","r",rownames(out))
                    out <- rbind(out, out2)
                }
            }else{
                out <- out[name.rho,,drop=FALSE]
            }

            
        }else{

            e.lmm <- mlmm(formula.mean, df = df, repetition = formula.repetition, data = dataL, structure = structure, control = list(optimizer = "FS"),                          
                          by = by, effects = "correlation", contrast.rbind = effects)
            out <- confint(e.lmm, df = df, columns = c("estimate","se","df","lower","upper","p.value"))

        }

        
    }else{
        ## *** UN/CS mixed model (no repetition)
        formula.repetition <- stats::as.formula(paste("~CCvariableCC|",name.id))
        if(heterogeneous){
            structure <- "UN"
        }else{
            structure <- "CS"
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
                out <- confint(anova(e.lmm, df = df, effects = Cmat, transform.rho = transform.rho),
                               method = method, columns = c("estimate","se","df","lower","upper","p.value"))
            }else{
                out0 <- confint(anova(e.lmm, df = df, effects = Cmat, transform.rho = "none"),
                               method = "none", columns = c("estimate","se","df","p.value"))
                out <- confint(anova(e.lmm, df = df, effects = Cmat, transform.rho = "atanh"),
                                method = method, columns = c("estimate","se","df","p.value"))
                out$estimate <- out0$estimate
                out$se <- out0$se
                out$df <- out0$df
                attr(out,"backtransform")[,"estimate"] <- FALSE
            }
        }else{
            e.lmm <- mlmm(formula.mean, df = df, repetition = formula.repetition, data = dataL, structure = structure,
                          by = by, effects = "correlation", contrast.rbind = effects)
            out <- confint(e.lmm, df = df, columns = c("estimate","se","df","lower","upper","p.value"))
        }
    }

    ## ** export
    attr(out, "lmm") <- e.lmm
    class(out) <- append("partialCor", class(out))
    return(out)
}


##----------------------------------------------------------------------
### partialCor.R ends here
