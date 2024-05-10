### effects.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 29 2024 (09:47) 
## Version: 
## Last-Updated: maj 10 2024 (16:02) 
##           By: Brice Ozenne
##     Update #: 572
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * effects.lmm (documentation)
##' @title Effects Derived For Linear Mixed Model
##' @description Estimate average counterfactual outcome or contrast of outcome from linear mixed models.
##' @param object a \code{lmm} object. 
##' @param variable [character] the variable relative to which the effect should be computed.
##' @param type [character] should the average counterfactual outcome for each variable level be evaluated (\code{"identity"})?
##' Or the difference in average counterfactual outcome between each pair of  variable level (\code{"difference"})?
##' Can have an second element to consider a transformation of the outcome:
##' the change from baseline (\code{"change"}),
##' area under the outcome curve (\code{"auc"}),
##' or area under the outcome curve minus baseline (\code{"auc-b"}).
##' @param newdata [data.frame] a dataset reflecting the covariate distribution relative to which the average outcome or contrast should be computed.
##' @param conditional [character] variable conditional to which the average conterfactual outcome or treatment effect should be computed.
##' @param rhs [numeric] the right hand side of the hypothesis.
##' @param repetition [character vector] repetition at which the effect should be assessed. By default it will be assessed at all repetitions.
##' @param multivariate [logical] should a multivariate Wald test be used to simultaneously test all null hypotheses.
##' @param prefix.time [character] When naming the estimates, text to be pasted before the value of the repetition variable.
##' Only relevant when \code{type = "aoc"} or \code{type = "ate"}.
##' @param prefix.var [logical] When naming the estimates, should the variable name be added or only the value?
##' @param sep.var [character] When naming the estimates, text to be pasted between the values to condition on.
##' Only relevant when \code{type = "aoc"} or \code{type = "ate"}.
##' @param ... Arguments passed to \code{\link{anova.lmm}}.
##'
##' @keywords htest
##' 
##' @examples
##' #### simulate data in the long format ####
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' #### Linear Mixed Model ####
##' eUN.lmm <- lmm(Y ~ visit + X1 + X2 + X5,
##'                repetition = ~visit|id, structure = "UN", data = dL)
##'
##' ## outcome
##' effects(eUN.lmm, variable = "X1")
##' effects(eUN.lmm, type = "difference", variable = "X1")
##' effects(eUN.lmm, type = "difference", variable = "X1", repetition = "3")
##'
##' ## change
##' effects(eUN.lmm, type = "change", variable = "X1")
##' effects(eUN.lmm, type = "change", variable = "X1", conditional = NULL)
##' effects(eUN.lmm, type = c("change","difference"), variable = "X1")
##'
##' ## auc
##' effects(eUN.lmm, type = "auc", variable = "X1")
##' effects(eUN.lmm, type = c("auc","difference"), variable = "X1")
##' 
##' #### fit Linear Mixed Model with interaction ####
##' dL$X1.factor <- as.factor(dL$X1)
##' dL$X2.factor <- as.factor(dL$X2)
##' eUN.lmmI <- lmm(Y ~ visit * X1.factor + X2.factor + X5,
##'                repetition = ~visit|id, structure = "UN", data = dL)
##'
##' ## average counterfactual conditional to a categorical covariate
##' effects(eUN.lmmI, variable = "X1.factor",
##'         conditional = c("X2.factor"), repetition = "3")
##' effects(eUN.lmmI, type = "change", variable = "X1.factor",
##'         conditional = c("X2.factor"), repetition = "3")
##' effects(eUN.lmmI, type = "auc", variable = "X1.factor",
##'         conditional = c("X2.factor"), repetition = "3")
##' 
##' ## average difference in counterfactual conditional to a categorical covariate
##' effects(eUN.lmmI, type = "difference", variable = "X1.factor",
##'         conditional = c("X2.factor"), repetition = "3")
##' effects(eUN.lmmI, type = c("change","difference"), variable = "X1.factor",
##'         conditional = c("X2.factor"), repetition = "3")
##' effects(eUN.lmmI, type = c("auc","difference"), variable = "X1.factor",
##'         conditional = c("X2.factor"), repetition = "3")
##' 
##' ## average difference in counterfactual conditional to a covariate
##' effects(eUN.lmmI, type = "difference", variable = "X1.factor",
##'         conditional = list(X5=0:2), repetition = "3")
##' effects(eUN.lmmI, type = c("difference","change"), variable = "X1.factor",
##'         conditional = list(X5=0:2))


## * effects.lmm (code)
##' @export
effects.lmm <- function(object, variable, newdata = NULL, type = c("identity","none"), conditional = NULL, rhs = NULL, repetition = NULL, multivariate = FALSE,
                        prefix.time = NULL, prefix.var = TRUE, sep.var = ",", ...){

    call <- match.call()

    ## ** check arguments
    if(!is.character(type)){
        stop("Argument \'type\' must be a character value. \n")
    }
    valid.type1 <- c("identity","difference")
    valid.type2 <- c("none","change","auc","auc-b")
    if(length(type)==1){
        if(type %in% valid.type2){
            type <- c(valid.type1[1],type)
        }else if(type %in% valid.type1){
            type <- c(type,valid.type2[1])
        }else{
            stop("Invalid value for argument \'type\'. It should either indicate: \n",
                 "- a contrast: \"",paste(valid.type1, collapse = "\" or \""),"\", \n",
                 "- a transformation of the outcome: \"",paste(valid.type2, collapse = "\" or \""),"\". \n")
        }
    }else if(length(type)==2){
        if(type[1] %in% valid.type1 && type[2] %in% valid.type2){
        }else if(type[2] %in% valid.type1 && type[1] %in% valid.type2){
            type <- rev(type)
        }else{
            stop("Invalid value(s): \"",paste0(setdiff(type, union(valid.type1,valid.type2)),collapse="\", \""),"\" for argument \'type\'. It should either indicate: \n",
                 "- a contrast: \"",paste(valid.type1, collapse = "\" or \""),"\", \n",
                 "- a transformation of the outcome: \"",paste(valid.type2, collapse = "\" or \""),"\". \n")
        }
        
    }else{
        stop("Argument \'type\' must have length 1 or 2. \n",
             "The first element indicating the contrast: \"",paste(valid.type1, collapse = "\" or \""),"\", \n",
             "and the second element indicating a transformation of the outcome: \"",paste(valid.type2, collapse = "\" or \""),"\". \n")
    }

    if(!is.null(variable) && !is.character(variable)){
        stop("Argument \'variable\' must be a character value or NULL. \n")
    }
    valid.variable <- attr(object$design$mean, "variable")
    rhs <- 0

    if(length(valid.variable) == 0){
        stop("Cannot evaluate the effect as no covariate was specified when fitting the model in the ",ifelse(type %in% c("ate","aco"),"mean",type)," structure. \n")
    }
    if(is.null(variable)){
        if(type[1]!="identity"){
            stop("Argument \'variable\' cannot be NULL when argument \'type\' contains \'difference\'. \n")
        }
    }else if(length(variable)!=1){
        stop("Argument \'variable\' must have length 1 or be NULL. \n")
    }
    if(!is.null(variable) && variable %in% valid.variable == FALSE){
        stop("Argument \'variable\' must be one of \"",paste(valid.variable,collapse ="\", \""),"\". \n")
    }
    ## ** normalized user input    

    ## find unique levels of the exposure variable
    if(!is.null(variable)){
        if(variable %in% names(object$xfactor$mean)){
            Ulevel.variable <- object$xfactor$mean[[variable]]
        }else{
            Ulevel.variable <- sort(unique(object$data[[variable]]))
        }
        if(length(Ulevel.variable)>2 && is.numeric(object$data[[variable]])){
            if(type[1] == "identity"){
                stop("Cannot compute average treatment effects for numeric variables (unless they only have two levels). \n")
            }else if(type[1] == "identity"){
                stop("Cannot compute average counterfactuals outcomes for numeric variables (unless they only have two levels). \n")
            }
        }
    }else{        
        Ulevel.variable <- ""
    }

    ## recover the time variable
    U.time <- object$time$levels
    n.time <- length(U.time)
    if(is.null(newdata)){
        data.augmented <- stats::model.frame(object, type = "add.NA", add.index = TRUE)
        if(any(is.na(data.augmented[manifest(object, effects = "mean")]))){
            warning("Missing value(s) among covariates used to estimate the average effects. \n",
                    "Corresponding lines in the dataset will be removed. \n",
                    "Consider specifying the argument \'newdata\' to specifying the marginal covariate distribution. \n")
        }
    }else{
        data.augmented <- stats::model.frame(object, newdata = newdata, add.index = TRUE, na.rm = FALSE)
    }
    
    ## only consider requested times
    if(!is.null(repetition)){
        if(!all(repetition %in%  U.time) &&  !all(repetition %in%  1:n.time)){
            stop("Argument \'repetition\' incorrect. Should either be: \n",
                 "- any of \"",paste0(U.time, collapse = "\", \""),"\" \n",
                 "- any integer from 1 to ",n.time," \n")
        }
        if(is.numeric(repetition)){
            repetition <- U.time[repetition]
        }
        if(type[2]=="change" && repetition %in% U.time[1] == FALSE){
            repetition <- c(U.time[1],repetition)
        }else if(type[2] %in% c("auc","auc-b") && any(U.time %in% repetition == FALSE)){
            message("Argument \'repetition\' ignored when argument \'type\' contains \"auc\" or \"auc-b\". \n")
            repetition <- U.time
        }
        data.augmented <- data.augmented[data.augmented$XXtimeXX %in% repetition,,drop=FALSE]
    }

    ## define conditioning set (normalize user input)
    if("conditional" %in% names(call) == FALSE){
        conditional <- switch(type[2],
                              "none" = object$time$var,
                              "change" = object$time$var,
                              "auc" = NULL,
                              "auc-b" = NULL)
    }

    if(is.null(conditional) || length(conditional)==0){
        conditional.time <- FALSE
        grid.conditional <- NULL
    }else if(is.character(conditional)){
        valid.conditional <- setdiff(c(object$time$var,names(object$xfactor$mean)),variable)
        if(any(conditional %in% valid.conditional == FALSE)){
            stop("Incorrect value for argument \'conditional\'. \n",
                 "If a character, it should only include the time variable and factor variables influencing the mean structure.\n",
                 "Valid values: \"",paste0(valid.conditional, collapse = "\" \""),"\".\n")
        }
        conditional.time <- object$time$var %in% conditional
        grid.conditional <- expand.grid(object$xfactor$mean[setdiff(conditional,object$time$var)])
    }else if(is.list(conditional)){
        valid.conditional <- setdiff(c(object$time$var,valid.variable),variable)
        if(is.null(names(conditional))){
            stop("Missing names for argument \'conditional\'. \n",
                 "Valid values: \"",paste0(valid.conditional, collapse = "\" \""),"\".\n")
        }
        if(any(names(conditional) %in% valid.conditional == FALSE)){
            stop("Incorrect names for argument \'conditional\'. \n",
                 "It should only include the time variable and variables influencing the mean structure.\n",
                 "Valid values: \"",paste0(valid.conditional, collapse = "\" \""),"\".\n")
        }
        conditional.time <- object$time$var %in% names(conditional)
        if(conditional.time & "repetition" %in% names(call) & !identical(sort(unname(conditional[[object$time$var]])),sort(unname(repetition)))){
            stop("Mismatch between argument \'repetition\' and \'conditional\' regarding timepoint values. \n")
        }
        grid.conditional <- expand.grid(conditional[setdiff(names(conditional),object$time$var)])
    }else if(!is.null(conditional)){
        stop("Argument \'conditional\' must either be NULL, or a character vector indicating factor variables, or a named list of covariate values. \n")
    }
    if(conditional.time && type[2] %in% c("auc","auc-b")){
        stop("Incompatibility between arguments \'conditional\' and \'type\'. \n",
             "Cannot condition on the time variable \"",object$time$var,"\" when type is \"auc\" or \"auc-b\". ")
    }
        

    if(length(grid.conditional)>0){            
        key.conditional <- nlme::collapse(grid.conditional, sep = sep.var)
        key.original <- levels(key.conditional)[match(nlme::collapse(data.augmented[names(grid.conditional)], sep =  sep.var), key.conditional)]
        if(prefix.var){
            if(length(grid.conditional)==1){
                data.augmented$XXstrataXX <- paste0(names(grid.conditional),"=",key.original)
            }else{
                newlevel <- nlme::collapse(as.data.frame(lapply(names(grid.conditional), function(iName){paste0(iName,"=",grid.conditional[,iName])})), sep = sep.var)
                data.augmented$XXstrataXX <- factor(key.original, levels = levels(key.conditional), labels = newlevel)
            }
        }else{
            data.augmented$XXstrataXX <- key.original
        }
        if(any(is.na(key.original))){
            data.augmented <- data.augmented[!is.na(key.original),,drop=FALSE]
        }
        U.strata <- unique(data.augmented$XXstrataXX)
        ## identify combination of timepoints and strata that never occur in the observed data (to be removed from the contrast matrix) 
        attr(U.strata,"table") <- table(data.augmented$XXtimeXX,data.augmented$XXstrataXX)
        if(!is.null(repetition)){
            attr(U.strata,"table") <- attr(U.strata,"table")[repetition,,drop=FALSE]
        }
        n.strata <- length(U.strata)
        
    }

    ## define contrast matrix
    if(is.null(prefix.time)){
        if(conditional.time){
            prefix.time <- switch(type[2],
                                  "none" = "t=",
                                  "change" = "dt=",
                                  "auc" = "auc",
                                  "auc-b" = "auc-b",
                                  NA)
        }else{
            prefix.time <- switch(type[2],
                                  "none" = "t",
                                  "change" = "dt",
                                  "auc" = "auc",
                                  "auc-b" = "auc-b",
                                  NA)
        }                              
    }

    if(n.time==1 && type[2] %in% c("change","auc","auc-b")){
        stop("Cannot evaluate ",type[2]," when there is a single timepoint. \n",
             "Considering setting argument \'type\' to \"static\" or specifying the argument \'repetition\' when calling lmm. \n")
    }

    M.contrast <- .effects_contrast(type = type, grid.conditional = grid.conditional, conditional.time = conditional.time, repetition = repetition, 
                                    U.strata = U.strata, n.strata = n.strata, 
                                    U.time = U.time, n.time = n.time,
                                    prefix.time = prefix.time, sep.var = sep.var)

    ## ** perform average
    if(type[1] == "identity"){
        ls.C <- lapply(Ulevel.variable, function(iLevel){ ## iLevel <- 0
            iData <- data.augmented
            ## iData <- data.augmented[data.augmented$subject==data.augmented$subject[1],]
            if(!is.null(variable)){
                iData[[variable]][] <- iLevel
            }

            iX <- stats::model.matrix(object$formula$mean.design, iData)[,colnames(object$design$mean),drop=FALSE] ## remove uncessary columns in case of (baseline) constraint

            if(length(grid.conditional)==0){
                iStrata <- droplevels(factor(iData$XXtimeXX, 
                                             levels = attr(M.contrast,"time.col")))
            }else{
                iStrata <- droplevels(factor(paste0(iData$XXtimeXX, sep.var, iData$XXstrataXX),
                                             levels = paste0(attr(M.contrast,"time.col"), sep.var, attr(M.contrast,"strata.col"))))
            }

            if(NROW(iX)!=NROW(iData)){
                ## handle missing values in covariates, i.e. rows removed by model.matrix
                iStrata <- iStrata[as.numeric(rownames(iX))]
            }
            iC <- M.contrast[,levels(iStrata),drop=FALSE] %*% do.call(rbind,by(iX,iStrata,colMeans))
            rownames(iC) <- paste0(iLevel,rownames(iC))
            attr(iC, "time") <- attr(M.contrast,"time.row")
            attr(iC, "strata") <- attr(M.contrast,"strata.row")
            attr(iC, "variable") <- rep(iLevel, NROW(M.contrast))
            return(iC)
        })
    }else if(type[1] == "difference"){
        Upair.variable <- unorderedPairs(Ulevel.variable, distinct = TRUE)

        ls.C <- lapply(1:NCOL(Upair.variable), function(iPair){ ## iPair <- 1
            iData1 <- data.augmented
            iData2 <- data.augmented

            iData1[[variable]][] <- Upair.variable[1,iPair]
            iData2[[variable]][] <- Upair.variable[2,iPair]

            iX1 <- stats::model.matrix(object$formula$mean.design, iData1)[,colnames(object$design$mean),drop=FALSE] ## remove uncessary columns in case of (baseline) constraint
            iX2 <- stats::model.matrix(object$formula$mean.design, iData2)[,colnames(object$design$mean),drop=FALSE] ## remove uncessary columns in case of (baseline) constraint
            
            if(length(grid.conditional)==0){
                iStrata <- droplevels(factor(iData1$XXtimeXX,
                                             levels = attr(M.contrast,"time.col")))
            }else{
                iStrata <- droplevels(factor(paste0(iData1$XXtimeXX,sep.var,iData1$XXstrataXX),
                                             levels = paste0(attr(M.contrast,"time.col"),sep.var,attr(M.contrast,"strata.col"))))
            }


            if(NROW(iX1)!=NROW(iData1)){
                ## handle missing values in covariates, i.e. rows removed by model.matrix
                iStrata <- iStrata[as.numeric(rownames(iX1))]
            }
            iC <- M.contrast[,levels(iStrata),drop=FALSE] %*% do.call(rbind,by(iX2-iX1,iStrata,colMeans))
            rownames(iC) <- paste0(Upair.variable[2,iPair],"-",Upair.variable[1,iPair],rownames(M.contrast))
            attr(iC, "time") <- attr(M.contrast,"time.row")
            attr(iC, "strata") <- attr(M.contrast,"strata.row")
            attr(iC, "variable") <- lapply(1:NROW(M.contrast), function(ii){Upair.variable[1:2,iPair]})
            return(iC)
        })
    }

    effect <- do.call(rbind,ls.C)
    effect.time <- unlist(lapply(ls.C,attr,"time"))
    attr(effect.time,"original") <- do.call(rbind,lapply(lapply(ls.C,attr,"time"),attr,"original"))
    effect.strata <- unlist(lapply(ls.C,attr,"strata"))
    effect.variable <- do.call(c,lapply(ls.C,attr,"variable"))

    if(prefix.var && !is.null(variable)){
        rownames(effect) <- paste0(variable,"=",rownames(effect))
    }

    ## ** check arguments
    out <- anova(object, effect = effect, multivariate = multivariate, ...)
    out$args$effect <- list(type)
    out$args$variable <- variable
    out$args$time <- object$time$var
    if(!is.null(effect.strata)){
        out$args$strata <- paste(colnames(grid.conditional), collapse = ", ")
    }
    add <- data.frame(matrix(NA, nrow = NROW(out$univariate)))
    names(add) <- variable
    add[] <- list(effect.variable)
    if(conditional.time){
        add <- cbind(add, stats::setNames(as.data.frame(effect.time), out$args$time))
    }
    if(!is.null(effect.strata)){
        add <- cbind(add, stats::setNames(as.data.frame(effect.strata), out$args$strata))
    }
    out$univariate <- cbind(add, out$univariate)
    rownames(out$univariate) <- rownames(effect)

    ## ** export
    attr(out,"class") <- append("effect_lmm",attr(out,"class"))
    return(out)
}

## * .effects_contrast
##' @description Generate contrast matrix
##' @noRd
.effects_contrast <- function(type, grid.conditional, conditional.time, repetition,
                              U.strata, n.strata, 
                              U.time, n.time,
                              prefix.time, sep.var){


    ## ** generate contrast matrix according to type across repetitions
    attr(U.time,"original") <- NULL
    if(type[2] == "none"){
        M.contrast <- diag(1, n.time, n.time)
        dimnames(M.contrast) <- list(U.time, U.time)
        if(!is.null(repetition)){
            M.contrast <- M.contrast[repetition, ,drop=FALSE]
        }
    }else if(type[2] == "change"){        
        M.contrast0 <- matrix(c(-1,rep(0,n.time-1)), byrow = TRUE, nrow = n.time, ncol = n.time) + diag(1, n.time, n.time)
        dimnames(M.contrast0) <- list(U.time, U.time)
        if(is.null(repetition)){
            M.contrast <- M.contrast0[-1,,drop=FALSE]
        }else{
            M.contrast <- M.contrast0[setdiff(repetition,U.time[1]),,drop=FALSE]
        }
    }else if(type[2] %in% c("auc","auc-b")){
        time.num <- as.numeric(U.time)
        if(any(is.na(time.num))){
            stop("Cannot evaluate the area under the curve of the outcome. \n",
                 "When calling lmm, argument \'repetition\'(=~time|cluster) must contain a numeric time variable. \n",
                 "Or a factor variable whose levels can be converted as numeric")
        }
        if(is.unsorted(time.num)){
            warning("The levels of the time variable do not correspond to numeric values in increasing order. \n",
                    "Can be an issue when evaluating the area under the curve.")
        }
        dtime.num <- diff(c(utils::head(time.num,1),time.num))/2 + diff(c(time.num,utils::tail(time.num,1)))/2
        if(type[2] == "auc-b"){
            dtime.num[1] <- dtime.num[1] - sum(dtime.num)
        }
        M.contrast <- matrix(dtime.num, byrow = TRUE, nrow = 1, ncol = n.time)
        dimnames(M.contrast) <- list(NULL, U.time)
    }

    ## ** possible stratification
    if(length(grid.conditional)>0){
 
        if(any(attr(U.strata,"table")==0)){
            ## in cross-over designs, there will be no observation at period drug=placebo under condition=active
            ## this case should be removed from the contrast matrix
            ls.contractC <- lapply(1:n.strata, function(iS){ ## iS <- 1
                iIndex <- which(attr(U.strata,"table")[,U.strata[iS]]>0)
                if(type[2] == "none"){
                    iM <- M.contrast[iIndex,iIndex,drop=FALSE]
                    attr(iM,"time.row") <- U.time[iIndex]
                }else if(type[2] == "change"){
                    if(1 %in% iIndex){
                        iM <- M.contrast[setdiff(iIndex,1)-1,iIndex,drop=FALSE]
                        attr(iM,"time.row") <- U.time[setdiff(iIndex,1)-1]
                    }else{
                        stop("Cannot evaluate the effect (no baseline measurement under condition ",U.strata[iS],"). \n")
                    }
                }else if(type[2] %in% c("auc","auc-b")){
                    iM <- M.contrast[,iIndex,drop=FALSE]
                }
                attr(iM,"time.col") <- U.time[iIndex]
                attr(iM,"strata.col") <- rep(U.strata[iS],length(iIndex))
                attr(iM,"strata.row") <- rep(U.strata[iS],NROW(iM))
                return(iM)
            })
            M.contrastC <- as.matrix(Matrix::bdiag(ls.contractC))
            attr(M.contrastC,"time.col") <- unlist(lapply(ls.contractC,"attr","time.col")) ## wil be NULL in the case of auc,auc-b
            attr(M.contrastC,"strata.col") <- unlist(lapply(ls.contractC,"attr","strata.col"))
            attr(M.contrastC,"time.row") <- unlist(lapply(ls.contractC,"attr","time.row"))
            attr(M.contrastC,"strata.row") <- unlist(lapply(ls.contractC,"attr","strata.row"))
        }else{
            M.contrastC <- as.matrix(Matrix::bdiag(replicate(n = n.strata, M.contrast, simplify = FALSE)))
            attr(M.contrastC,"time.col") <- rep(U.time, times = n.strata) ## wil be NULL in the case of auc,auc-b
            attr(M.contrastC,"strata.col") <- unlist(lapply(U.strata, rep, times = NCOL(M.contrast)))
            attr(M.contrastC,"time.row") <- rep(rownames(M.contrast), times = n.strata) ## wil be NULL in the case of auc,auc-b
            attr(M.contrastC,"strata.row") <- unlist(lapply(U.strata, rep, times = NROW(M.contrast)))
        }
        colnames(M.contrastC) <- paste0(attr(M.contrastC,"time.col"), sep.var, attr(M.contrastC,"strata.col"))
        
        if(type[2] %in% c("none","change")){
            rownames(M.contrastC) <- paste0("(",prefix.time,attr(M.contrastC,"time.row"),sep.var,attr(M.contrastC,"strata.row"),")")
        }else if(type[2]%in%c("auc","auc-b")){
            if(nchar(prefix.time)==0){
                rownames(M.contrastC) <- paste0("(",attr(M.contrastC,"strata.row"),")")
            }else{
                rownames(M.contrastC) <- paste0("(",prefix.time,sep.var,attr(M.contrastC,"strata.row"),")")
            }
        }

        
    }else{
        M.contrastC <- M.contrast
        attr(M.contrastC,"time.row") <- rownames(M.contrast) ## wil be NULL in the case of auc,auc-b
        attr(M.contrastC,"time.col") <- U.time
        colnames(M.contrastC) <- attr(M.contrastC,"time.col")

        if(type[2] %in% c("none","change")){
            rownames(M.contrastC) <- paste0("(",prefix.time,attr(M.contrastC,"time.row"),")")
        }else if(type[2] %in% c("auc","auc-b")){
            if(nchar(prefix.time)==0){
                rownames(M.contrastC) <- ""
            }else{
                rownames(M.contrastC) <- paste0("(",prefix.time,")")
            }
        }
    }

    ## ** collapse across timepoints
    if(conditional.time || type[2] %in% c("auc","auc-b")){
        out <- M.contrastC
    }else if(length(grid.conditional)==0 && type[2] %in% c("none","change")){
        out <- rbind(colMeans(M.contrastC))
        if(nchar(prefix.time)==0){
            rownames(out) <- ""
        }else{
            rownames(out) <- paste0("(",prefix.time,")")
        }
        attr(out,"time.col") <- attr(M.contrastC,"time.col")
    }else if(length(grid.conditional)>0 && type[2] %in% c("none","change")){
        ls.out <- by(M.contrastC, attr(M.contrastC,"strata.row"), colMeans)
        out <- do.call(rbind,ls.out)
        if(nchar(prefix.time)==0){
            rownames(out) <- paste0("(",names(ls.out),")")
        }else{
            rownames(out) <- paste0("(",prefix.time,sep.var,names(ls.out),")")
        }
        attr(out,"time.col") <- attr(M.contrastC,"time.col")
        attr(out,"strata.row") <- names(ls.out)
        attr(out,"strata.col") <- attr(M.contrastC,"strata.col")
    }
    ## ** export
    return(out)
}
##----------------------------------------------------------------------
### effects.R ends here

