### partialCor.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May  1 2022 (17:01) 
## Version: 
## Last-Updated: Jun  2 2022 (14:53) 
##           By: Brice Ozenne
##     Update #: 93
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
##' @param repetition [formula] Specify the structure of the data: the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' @param data [data.frame] dataset containing the variables.
##'
##' @details Fit a mixed model to estimate the partial correlation which can be time consuming.
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
##' #### bivariate (with repetition) ####
##' data(gastricbypassL, package = "LMMstar")
##' 
##' partialCor(weight+glucagonAUC~time, data = gastricbypassL)
##' 
##' partialCor(weight+glucagonAUC~time, repetition =~time|id, data = gastricbypassL)

## * partialCor (documentation)
##' @export
partialCor <- function(formula, data, repetition){

    ## ** normalize arguments
    data <- as.data.frame(data)

    ## *** convert first argument to list of formula
    if(inherits(formula,"formula")){
        rhs <- stats::delete.response(stats::terms(formula))
        response <- setdiff(all.vars(formula),all.vars(rhs))

        if(length(response)!=2){
            stop("Argument \'formula\' should contain exactly two variables on the left hand side of the formula. \n")
        }
        formula <- lapply(response, function(iY){stats::as.formula(paste(iY,deparse(rhs)))}) ## iY <- response[1]        
    }

    ## *** check formula agree with data
    if(length(formula)!=2){
        stop("Argument \'formula\' should be a list of two elements. \n")
    }
    if(any(unlist(lapply(formula, inherits, "formula"))==FALSE)){
        stop("Argument \'formula\' should be a list of formula. \n")
    }

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
    
    ## ** reshape    
    ls.name.X <- lapply(formula, function(iF){all.vars(stats::delete.response(stats::terms(iF)))})
    name.X <- unique(unlist(ls.name.X))

    ls.name.Y <- lapply(ls.name.XY, function(iF){setdiff(iF,c(name.X,name.id))})
    name.Y <- unlist(ls.name.Y)
    if(any(duplicated(name.Y))){
        stop("Variables in the left hand side of argument should be unique. \n")
    }
    dataL <- stats::reshape(data[, unique(c(name.XY, name.id, name.time)),drop=FALSE], direction  = "long",
                            idvar = c(name.id, name.time),
                            varying = name.Y,
                            v.names = "CCvalueCC",
                            timevar = "CCvariableCC")
    dataL$CCvariableCC <- factor(dataL$CCvariableCC, labels = name.Y)
    rownames(dataL) <- NULL
    
    ## ** fit mixed model
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
    if(any(test.duplicated)){
        dataL <- dataL[order(dataL$id,dataL$CCvariableCC),]
        dataL$CCrepetitionCC <- unlist(tapply(dataL$id,dataL$id,function(iId){1:length(iId)}))
        formula.repetition <- stats::as.formula(paste("~CCrepetitionCC|",name.id))

        e.lmm <- lmm(formula.mean, repetition = formula.repetition,
                     data = dataL, structure = CS(~CCvariableCC, heterogeneous = TRUE),
                     control = list(optimizer = "FS"))
        out <- model.tables(e.lmm, effects = "correlation")

        name.cor <- paste0("rho(",name.Y,",",rev(name.Y),")")
        index.cor <- which(rownames(out) %in% name.cor)
        out <- out[index.cor,,drop=FALSE]
        attr(out,"old2new") <- attr(out,"old2new")[index.cor]
        attr(attr(out,"old2new"),"names") <- attr(attr(out,"old2new"),"names")
        attr(out,"backtransform.names") <- attr(out,"backtransform.names")[index.cor]
    }else{
        formula.repetition <- stats::as.formula(paste("~CCvariableCC|",name.id))

        e.lmm <- lmm(formula.mean, repetition = formula.repetition,
                     data = dataL, structure = "UN")
        out <- model.tables(e.lmm, effects = "correlation")
    }

    ## ** export
    attr(out, "lmm") <- e.lmm
    return(out)
}


##----------------------------------------------------------------------
### partialCor.R ends here
