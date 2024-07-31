### effects.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 29 2024 (09:47) 
## Version: 
## Last-Updated: jul 31 2024 (10:40) 
##           By: Brice Ozenne
##     Update #: 1124
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
##' 
##' @param object a \code{lmm} object. 
##' @param variable [character/list] exposure variable relative to which the effect should be computed.
##' Can also be a list with two elements: the first being the variable (i.e. a character) and the second the levels or values for this variable to be considered.
##' @param effects [character] should the average counterfactual outcome for each variable level be evaluated (\code{"identity"})?
##' Or the difference in average counterfactual outcome between each pair of variable level (\code{"difference"})?
##' @param type [character/numeric vector] Possible transformation of the outcome: no transformation (\code{"outcome"}),
##' change from baseline (\code{"change"}),
##' area under the outcome curve (\code{"auc"}),
##' or area under the outcome curve minus baseline (\code{"auc-b"}).
##' Alternatively can be a numeric vector indicating how to weight each timepoint. 
##' @param repetition [character vector] repetition at which the effect should be assessed. By default it will be assessed at all repetitions.
##' @param conditional [character/data.frame] variable(s) conditional to which the average conterfactual outcome or treatment effect should be computed.
##' Alternatively can also be a data.frame where each column correspond to a variable and the rows to the level of the variable(s).
##' @param ref.repetition [numeric or character] index or value of the reference level for the repetition variable.
##' Only relevant when \code{type} equal to \code{"change"}.
##' Can be \code{NA} to evaluate change relative to all possible reference levels.
##' @param ref.variable [numeric or character] index or value of the reference level for the exposure variable.
##' Only relevant when \code{effects} equal to \code{"difference"}.
##' Can be \code{NA} to evaluate the difference relative to all possible reference levels.
##' @param newdata [data.frame] a dataset reflecting the covariate distribution relative to which the average outcome or contrast should be computed.
##' @param rhs [numeric] the right hand side of the hypothesis.
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
##' @details The uncertainty is quantified assuming the contrast matrix to be a-priori known.
##' Said otherwise the standard error does not account for the uncertainty about the covariate distribution.
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
##' e.YbyX1 <- effects(eUN.lmm, variable = "X1")
##' e.YbyX1
##' summary(e.YbyX1)
##' model.tables(e.YbyX1)
##' coef(e.YbyX1, type = "contrast")
##' effects(eUN.lmm, effects = "difference", variable = "X1")
##' effects(eUN.lmm, effects = "difference", variable = "X1", repetition = "3")
##'
##' ## change
##' effects(eUN.lmm, type = "change", variable = "X1")
##' effects(eUN.lmm, type = "change", variable = "X1", ref.repetition = 2)
##' effects(eUN.lmm, type = "change", variable = "X1", conditional = NULL)
##' effects(eUN.lmm, type = "change", effects = "difference", variable = "X1")
##'
##' ## auc
##' effects(eUN.lmm, type = "auc", variable = "X1")
##' effects(eUN.lmm, type = "auc", effects = "difference", variable = "X1")
##' 
##' #### fit Linear Mixed Model with interaction ####
##' dL$X1.factor <- as.factor(dL$X1)
##' dL$X2.factor <- as.factor(dL$X2)
##' eUN.lmmI <- lmm(Y ~ visit * X1.factor + X2.factor + X5,
##'                repetition = ~visit|id, structure = "UN", data = dL)
##'
##' ## average counterfactual conditional to a categorical covariate
##' effects(eUN.lmmI, variable = "X1.factor",
##'         conditional = "X2.factor", repetition = "3")
##' effects(eUN.lmmI, type = "change", variable = "X1.factor",
##'         conditional = "X2.factor", repetition = "3")
##' effects(eUN.lmmI, type = "auc", variable = "X1.factor", conditional = "X2.factor")
##' 
##' ## average difference in counterfactual conditional to a categorical covariate
##' effects(eUN.lmmI, effects = "difference", variable = "X1.factor",
##'         conditional = c("X2.factor"), repetition = "3")
##' effects(eUN.lmmI, effects = "difference", type = "change", variable = "X1.factor",
##'         conditional = c("X2.factor"), repetition = "3")
##' effects(eUN.lmmI, effects = "difference", type = "auc", variable = "X1.factor",
##'         conditional = "X2.factor")
##' 
##' ## average difference in counterfactual conditional to a covariate
##' effects(eUN.lmmI, effect = "difference", variable = "X1.factor",
##'         conditional = data.frame(X5=0:2), repetition = "3")
##' effects(eUN.lmmI, effect = "difference", type = "change", variable = "X1.factor",
##'         conditional = data.frame(X5=0:2))


## * effects.lmm (code)
##' @export
effects.lmm <- function(object, variable, effects = "identity", type = "outcome",
                        repetition = NULL, conditional = NULL, ref.repetition = 1, ref.variable = 1, 
                        newdata = NULL, rhs = NULL, multivariate = FALSE,
                        prefix.time = NULL, prefix.var = TRUE, sep.var = ",", ...){

    call <- match.call()
    time.var <- object$time$var
    alltime.var <- attr(time.var,"original")

    ## ** check arguments
    ## *** type
    if(inherits(type,"data.frame")){
        stop("Argument \'type\' should be a character, numeric vector, or numeric matrix. \n")
    }else if(is.numeric(type)){
        if(is.vector(type)){
            weight.type <- matrix(type, nrow = 1)
        }else if(is.matrix(type)){
            if(is.null(rownames(type)) || any(duplicated(rownames(type)))){
                stop("Argument \'type\' should have (distinct) rownames when a numeric matrix. \n")
            }
            weight.type <- type            
        }else{
            stop("Argument \'type\' should be a character, numeric vector, or numeric matrix. \n")
        }
        type <- "user"
    }else{
        weight.type <- NULL
        if(!is.character(type)){
            stop("Argument \'type\' should be a character or a numeric vector.")
        }
        if(length(type)>1){
            stop("When a character, argument \'type\' should have length 1. \n")
        }
        valid.type <- c("outcome","change","auc","auc-b","user")
        if(any(type %in% valid.type == FALSE)){
            stop("Incorrect value for argument \'type\': \"",paste(setdiff(type,valid.type), collapse ="\", \""),"\". \n",
                 "Valid values: \"",paste(valid.type, collapse ="\", \""),"\". \n")
        }
    }
    
    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character vector. \n")
    }
    if(length(effects)!=1){
        stop("Argument \'effects\' must have length 1. \n")
    }
    valid.effects <- c("identity","difference")
    if(effects %in% valid.effects == FALSE){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }
    
    if(!is.null(variable)){
        if(is.list(variable)){
            if(length(variable)!=2){
                stop("Argument \'variable\' must have length 2 when a list. \n",
                     "The first element being the variable and the second its levels to be considered. \n")
            }
            if(!is.character(variable[[1]])){
                stop("Argument \'variable\' must have length 2 when a list. \n",
                     "The first element being the variable and the second its levels to be considered. \n")
            }
            Ulevel.variable <- variable[[2]]
            variable <- variable[[1]]
        }else if(!is.character(variable)){
            stop("Argument \'variable\' must be a character value, a list, or NULL. \n")
        }else{
            Ulevel.variable <- NULL
        }
    } 
    valid.variable <- attr(object$design$mean, "variable")
    rhs <- 0

    if(length(valid.variable) == 0){
        stop("Cannot evaluate the effect as no covariate was specified when fitting the model in the mean structure. \n")
    }
    if(is.null(variable)){
        if(effects=="difference"){
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
        
        if(!is.null(Ulevel.variable)){
            if(any(is.na(Ulevel.variable))){
                stop("Argument \'variable\' is incorrect: proposed values should not contain NAs. \n")
            }
            if(is.character(object$data[[variable]]) && any(Ulevel.variable %in% unique(object$data[[variable]]) == FALSE)){
                stop("Argument \'variable\' is incorrect: proposed values do not match the values taken by the variable \"",variable,"\" in the dataset. \n")
            }else if(is.factor(object$data[[variable]]) && any(Ulevel.variable %in% levels(object$data[[variable]]) == FALSE)){
                stop("Argument \'variable\' is incorrect: proposed values do not match the values taken by the variable \"",variable,"\" in the dataset. \n")
            }else if(min(Ulevel.variable) < min(object$data[[variable]], na.rm = TRUE) || max(Ulevel.variable) > max(object$data[[variable]], na.rm = TRUE)){
                message("Some values specified in argument \'variable\' are outside the range of values taken by the variable \"",variable,"\" in the dataset. \n")
            }
        }else{
            if(variable %in% names(object$xfactor$mean)){
                Ulevel.variable <- object$xfactor$mean[[variable]]
            }else{
                Ulevel.variable <- sort(unique(object$data[[variable]]))
            }
            if(length(Ulevel.variable)>2 && is.numeric(object$data[[variable]])){
                if(effects == "difference"){
                    message("Cannot compute average treatment effects for numeric variables (unless they only have two levels). \n",
                            "Will contrast 1 to 0. \n")
                }else if(effects == "identity"){
                    message("Cannot compute average counterfactuals outcomes for numeric variables (unless they only have two levels). \n",
                            "Will contrast 1 to 0. \n")
                }
                Ulevel.variable <- 0:1
            }
        }

        if(is.character(ref.variable)){
            if(any(duplicated(ref.variable)) || any(ref.variable %in% Ulevel.variable == FALSE)){
                stop("Incorrect value in argument \'ref.variable\': \"",paste(ref.variable[ref.variable %in% Ulevel.variable == FALSE], collapse = "\" \""),"\". \n",
                     "Valid values: \"",paste(Ulevel.variable, collapse = "\" \""),"\". \n")
            }
            ref.variable <- which(Ulevel.variable %in% ref.variable)
        }else if(is.numeric(ref.variable)){
            if(any(ref.variable %in% 1:length(Ulevel.variable) == FALSE)){
                stop("Incorrect value in argument \'ref.variable\': \"",paste(ref.variable[ref.variable %in% Ulevel.variable == FALSE], collapse = "\" \""),"\". \n",
                     "Valid values: integer from 1 to ",length(Ulevel.variable),". \n")
            }
        }else if(is.null(ref.variable) || all(is.na(ref.variable))){
             ref.variable <- 1:length(Ulevel.variable)
         }else{
             stop("Incorrect specification of argument \'ref.variable\'. Should be numeric, character, or NA. \n")
         }
    }else{        
        Ulevel.variable <- ""
    }

    ## recover the time variable
    U.time <- object$time$levels
    n.time <- length(U.time)
    if(is.null(newdata)){
        data.augmented <- stats::model.frame(object, type = "add.NA", add.index = TRUE)
        if(any(is.na(data.augmented[stats::variable.names(object, effects = "mean")]))){
            warning("Missing value(s) among covariates used to estimate the average effects. \n",
                    "Corresponding lines in the dataset will be removed. \n",
                    "Consider specifying the argument \'newdata\' to specifying the marginal covariate distribution. \n")
        }
    }else{
        data.augmented <- stats::model.frame(object, newdata = newdata, add.index = TRUE, na.rm = FALSE)
    }
    data.augmented$XXstrataXX <- 1

    ## ref.repetition
    if(is.character(ref.repetition)){
        if(any(ref.repetition %in% U.time == FALSE)){
            stop("Incorrect value in argument \'ref.repetition\': \"",paste(ref.repetition[ref.repetition %in% U.time == FALSE], collapse = "\" \""),"\". \n",
                 "Valid values: \"",paste(U.time, collapse = "\" \""),"\". \n")
        }
        ref.repetition <- which(U.time %in% ref.repetition)
    }else if(is.numeric(ref.repetition)){
        if(any(ref.repetition %in% 1:length(U.time) == FALSE)){
            stop("Incorrect value in argument \'ref.repetition\': \"",paste(ref.repetition[ref.repetition %in% U.time == FALSE], collapse = "\" \""),"\". \n",
                 "Valid values: integer from 1 to ",length(U.time),". \n")
        }
    }else if(is.null(ref.repetition) || all(is.na(ref.repetition))){
        ref.repetition <- 1:length(U.time)
    }else{
        stop("Incorrect specification of argument \'ref.repetition\'. Should be numeric, character, or NA. \n")
    }

    ## identify requested times
    if(!is.null(repetition)){
        if(!all(repetition %in%  U.time) &&  !all(repetition %in%  1:n.time)){
            stop("Argument \'repetition\' incorrect. Should either be: \n",
                 "- any of \"",paste0(U.time, collapse = "\", \""),"\" \n",
                 "- any integer from 1 to ",n.time," \n")
        }
        if(is.numeric(repetition)){
            repetition <- U.time[repetition]
        }
        if((type=="change") && any(U.time[ref.repetition] %in% repetition == FALSE)){
            repetition <- union(U.time[ref.repetition],repetition)
        }
    }

    ## define conditioning set (normalize user input)
    if("conditional" %in% names(call) == FALSE){
        if(type %in% c("outcome","change")){
            conditional <- alltime.var
        }else{
            conditional <- NULL
        }
    }

    if(is.null(conditional) || length(conditional)==0){
        conditional.time <- FALSE
        grid.conditional <- NULL
    }else if(is.character(conditional)){
        valid.conditional <- union(c(alltime.var,time.var),setdiff(names(object$xfactor$mean),variable))
        if(any(conditional %in% valid.conditional == FALSE)){
            stop("Incorrect value for argument \'conditional\'. \n",
                 "If a character, it should only include the time variable and factor variables influencing the mean structure.\n",
                 "Valid values: \"",paste0(valid.conditional, collapse = "\" \""),"\".\n")
        }        
        conditional.time <- (type %in% c("outcome","change")) && ((time.var %in% conditional) || all(alltime.var %in% c(variable,conditional)))
        grid.conditional <- expand.grid(object$xfactor$mean[setdiff(conditional,c(alltime.var,time.var))])
    }else if(is.data.frame(conditional)){
        if(is.null(repetition)){
            valid.conditional <- union(c(alltime.var,time.var),setdiff(valid.variable,variable))
        }else{
            ## cannot allow multiple time variable as repetition only takes vectors in
            valid.conditional <- union(time.var,setdiff(valid.variable,variable))
        }
        if(is.null(names(conditional))){
            stop("Missing names for argument \'conditional\'. \n",
                 "Valid values: \"",paste0(valid.conditional, collapse = "\" \""),"\".\n")
        }
        if(any(names(conditional) %in% valid.conditional == FALSE)){
            stop("Incorrect names for argument \'conditional\'. \n",
                 "It should only include the time variable and variables influencing the mean structure.\n",
                 "Valid values: \"",paste0(valid.conditional, collapse = "\" \""),"\".\n")
        }

        conditional.time <- (type %in% c("outcome","change")) && ((time.var %in% names(conditional)) || all(alltime.var %in% c(variable,names(conditional))))
        if((time.var %in% names(conditional)) || all(alltime.var %in% c(variable,names(conditional)))){
            if(!is.null(repetition)){
                if(!identical(sort(unname(conditional[[time.var]])),sort(unname(repetition)))){
                    stop("Mismatch between argument \'repetition\' and \'conditional\' regarding timepoint values. \n")
                }
            }else{
                if(time.var %in% names(conditional)){
                    repetition <- unique(conditional[[time.var]])
                    if(!all(repetition %in%  U.time) &&  !all(repetition %in%  1:n.time)){
                        stop("Argument \'repetition\' incorrect. Should either be: \n",
                             "- any of \"",paste0(U.time, collapse = "\", \""),"\" \n",
                             "- any integer from 1 to ",n.time," \n")
                    }
                    if(is.numeric(repetition)){
                        repetition <- U.time[repetition]
                    }
                }else{
                    repetition <- merge(x = unique(cbind(stats::setNames(list(Ulevel.variable),variable),conditional)[alltime.var]),
                                        y = cbind(attr(object$time$levels,"original"),levels = object$time$levels),
                                        by = alltime.var, all.x = TRUE)$levels
                    if(any(is.na(repetition))){
                        stop("Argument \'repetition\' incorrect: values do not match those stored in the lmm. \n")
                    }
                }
            }
            if((type=="change") && any(U.time[ref.repetition] %in% repetition == FALSE)){
                repetition <- union(U.time[ref.repetition],repetition)
            }
        }
        grid.conditional <- unique(conditional[setdiff(names(conditional),c(alltime.var,time.var))])

        test.positivity <- interaction(grid.conditional) %in% interaction(object$data[names(grid.conditional)])
        if(all(test.positivity==FALSE)){
            stop("No observation has the covariate values specified in argument \'conditional\'. \n")
        }else if(any(test.positivity==FALSE)){
            warning("Some of the covariate values specified in argument \'conditional\' do not exist in the dataset used to fit the linear mixed model (",round(100*mean(test.positivity==FALSE),3),"%) . \n")
        }
    }else if(!is.null(conditional)){
        stop("Argument \'conditional\' must either be NULL, or a character vector indicating factor variables, or a named list of covariate values. \n")
    }

    ## restrict to requested times
    if(!is.null(repetition)){
        data.augmented <- data.augmented[data.augmented$XXtimeXX %in% repetition,,drop=FALSE]
        if(type=="user"){
            if(length(weight.type) != length(repetition)){
                if(length(weight.type) == n.time && all(weight.type[U.time %in% repetition == FALSE]==0)){
                    weight.type <- weight.type[,U.time %in% repetition == FALSE,drop=FALSE]
                }else{
                    stop("Incorrect length for argument \'type\'. \n",
                         "When numeric should have the same length as the argument \'repetition\' (here ",length(repetition),"). \n")
                }
            }
            if(!is.null(colnames(weight.type)) & any(repetition %in% colnames(weight.type) == FALSE)){
                stop("Incorrect column names for argument \'type\'. \n",
                     "When a matrix should either have the repetition levels as column names or no column names. \n")
            }
        }
    }else if(type=="user"){
        if(NCOL(weight.type) != n.time){
            stop("Incorrect length for argument \'type\'. \n",
                 "When numeric should have the same length as the number of repetitions (here ",n.time,"). \n")
        }
        if(!is.null(colnames(weight.type)) & any(U.time %in% colnames(weight.type) == FALSE)){
            stop("Incorrect column names for argument \'type\'. \n",
                 "When a matrix should either have the repetition levels as column names or no column names. \n")
        }
    }

    ## combine strata variables
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
    }
    U.strata <- sort(unique(data.augmented$XXstrataXX))
    n.strata <- length(U.strata)        

    ## identify number of observation at each timepoint for each strata
    ## will be used to exclude from the contrast matrix combination of timepoints and strata that never occur in the observed data
    nobs.time <- table(time = data.augmented$XXtimeXX,strata = data.augmented$XXstrataXX)

    if(type %in% c("change","auc","auc-b","user")){
        if(is.null(repetition)){
            data.augmented.rep <- data.augmented
        }else{
            data.augmented.rep <- data.augmented[data.augmented$XXtimeXX %in% repetition,,drop=FALSE]
        }
        if(any(table(droplevels(data.augmented.rep$XXtimeXX), strata = droplevels(data.augmented.rep$XXclusterXX))!=1)){
            stop("Incorrect argument \'newdata\': no missing data allowed when argument \'type\' is ",ifelse(type=="user","user defined",type),". \n",
                 "(all clusters should have one line for each repetition)")
        }
    }

    ## ** individual specific contrast matrix
    if(is.null(prefix.time)){

        if(length(attr(time.var,"original"))>1){
            alltime.var <- paste(attr(time.var,"original"),collapse=".")
        }else{
            alltime.var <- time.var
        }

        if(conditional.time){
            prefix.time <- switch(type,
                                  "outcome" = paste0(alltime.var,"="),
                                  "change" = paste0("\u0394",alltime.var,"="),
                                  "auc" = "auc",
                                  "auc-b" = "auc-b",
                                  "user" = "",
                                  NA)
        }else{
            prefix.time <- switch(type,
                                  "outcome" = paste0(alltime.var),
                                  "change" = paste0("\u0394",alltime.var),
                                  "auc" = "auc",
                                  "auc-b" = "auc-b",
                                  "user" = "",
                                  NA)
        }                              
    }

    if((n.time==1 || length(repetition)==1) && type %in% c("change","auc","auc-b")){
        stop("Cannot evaluate ",type," when there is a single timepoint. \n",
             "Consider modifying the argument \'repetition\' when calling lmm. \n")
    }

    M.contrast <- .effects_contrast(type = type, weight.type = weight.type, grid.conditional = grid.conditional, conditional.time = conditional.time, repetition = repetition, 
                                    U.strata = U.strata, n.strata = n.strata, 
                                    U.time = U.time, n.time = n.time, nobs.time = nobs.time,
                                    ref.repetition = ref.repetition, prefix.time = prefix.time, sep.var = sep.var)

    ## ** contrast matrix 
    if(effects == "identity"){
        ls.C <- lapply(Ulevel.variable, function(iLevel){ ## iLevel <- 1
            iData <- data.augmented
            ## iData <- data.augmented[data.augmented$subject==data.augmented$subject[1],]
            if(!is.null(variable)){
                if(variable %in% names(object$xfactor$mean)){
                    iData[[variable]] <- factor(iLevel, object$xfactor$mean[[variable]])
                }else{
                    iData[[variable]] <- iLevel
                }
                
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

            ## empirical average of the covariate values and apply contrast
            iC <- M.contrast[,levels(iStrata),drop=FALSE] %*% do.call(rbind,by(iX,iStrata,colMeans))
            ## average over times
            if(conditional.time==FALSE && type %in% c("auc","auc-b","user") == FALSE && NROW(iC)>n.strata){
                iC.strata <- lapply(U.strata, function(iStrata){ ## iStrata <- U.strata[1]
                    iIndex <- which(attr(M.contrast,"strata.row")==iStrata)
                    iiC <- matrix(colSums(iC * attr(M.contrast,"weight")[iStrata,]),
                                  nrow = 1, dimnames = list(attr(M.contrast,"weight.name")[iStrata],colnames(iC)))
                    attr(iiC,"n.row") <- sum(attr(M.contrast,"n.row")[iIndex])
                    attr(iiC,"time.row") <- attr(M.contrast,"time.row")
                    attr(iiC,"strata.row") <- iStrata
                    return(iiC)
                })
                iC <- do.call(rbind,iC.strata)
                attr(iC, "n") <- sapply(iC.strata,attr,"n.row")
                attr(iC, "time") <- lapply(iC.strata,attr,"time.row")
                attr(iC, "strata") <- sapply(iC.strata,attr,"strata.row")
            }else{
                attr(iC, "n") <- attr(M.contrast,"n.row")
                attr(iC, "time") <- attr(M.contrast,"time.row")
                attr(iC, "strata") <- attr(M.contrast,"strata.row")
            }
            if(iLevel!=""){
                if(any(rownames(iC)!="")){
                    rownames(iC) <- paste0(iLevel,"(",rownames(iC),")")
                }else{
                    rownames(iC) <- iLevel
                }
            }
            attr(iC, "variable") <- rep(iLevel, NROW(iC))
            return(iC)
        })
        
    }else if(effects == "difference"){
        
        Upair.variable <- unorderedPairs(Ulevel.variable[union(ref.variable,1:length(Ulevel.variable))], distinct = TRUE)
        if(length(ref.variable)<length(Ulevel.variable)){
            Upair.variable <- Upair.variable[,Upair.variable[1,] %in% Ulevel.variable[ref.variable],drop=FALSE]
        }

        ls.C <- lapply(1:NCOL(Upair.variable), function(iPair){ ## iPair <- 1
            iData1 <- data.augmented
            iData2 <- data.augmented

            if(variable %in% names(object$xfactor$mean)){
                iData1[[variable]] <- factor(Upair.variable[1,iPair], levels = object$xfactor$mean[[variable]])
                iData2[[variable]] <- factor(Upair.variable[2,iPair], levels = object$xfactor$mean[[variable]])
            }else{
                iData1[[variable]] <- Upair.variable[1,iPair]
                iData2[[variable]] <- Upair.variable[2,iPair]
            }

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
            ## empirical average of the covariate values and apply contrast
            iC <- M.contrast[,levels(iStrata),drop=FALSE] %*% do.call(rbind,by(iX2-iX1,iStrata,colMeans))
            ## average over times
            if(conditional.time==FALSE && type %in% c("auc","auc-b","user") == FALSE && NROW(iC)>n.strata){
                iC.strata <- lapply(U.strata, function(iStrata){ ## iStrata <- U.strata[1]
                    iIndex <- which(attr(M.contrast,"strata.row")==iStrata)
                    iiC <- matrix(colSums(iC * attr(M.contrast,"weight")[iStrata,]),
                                  nrow = 1, dimnames = list(attr(M.contrast,"weight.name")[iStrata],colnames(iC)))
                    attr(iiC,"n.row") <- sum(attr(M.contrast,"n.row")[iIndex])
                    attr(iiC,"time.row") <- attr(M.contrast,"time.row")
                    attr(iiC,"strata.row") <- iStrata
                    return(iiC)
                })
                iC <- do.call(rbind,iC.strata)
                attr(iC, "n") <- sapply(iC.strata,attr,"n.row")
                attr(iC, "time") <- lapply(iC.strata,attr,"time.row")
                attr(iC, "strata") <- sapply(iC.strata,attr,"strata.row")
            }else{
                attr(iC, "n") <- attr(M.contrast,"n.row")
                attr(iC, "time") <- attr(M.contrast,"time.row")
                attr(iC, "strata") <- attr(M.contrast,"strata.row")
            }
            if(any(rownames(iC)!="")){
                rownames(iC) <- paste0(Upair.variable[2,iPair],"-",Upair.variable[1,iPair],"(",rownames(iC),")")
            }else{
                rownames(iC) <- paste0(Upair.variable[2,iPair],"-",Upair.variable[1,iPair])
            }
            
            attr(iC, "variable") <- lapply(1:NROW(iC), function(ii){Upair.variable[1:2,iPair]})
            return(iC)
        })
    }

    M.effect <- do.call(rbind,ls.C)
    effect.n <- unlist(lapply(ls.C,attr,"n"))
    effect.time <- unlist(lapply(ls.C,attr,"time"))
    effect.strata <- unlist(lapply(ls.C,attr,"strata"))
    effect.variable <- do.call(c,lapply(ls.C,attr,"variable"))

    if(prefix.var && !is.null(variable)){
        rownames(M.effect) <- paste0(variable,"=",rownames(M.effect))
    }

    ## ** test linear hypotheses
    out <- anova(object, effect = M.effect, multivariate = multivariate, ...)
    out$args$effect <- effects
    out$args$type <- type
    out$args$variable <- variable
    if(effects=="difference"){
        out$args$ref.variable <- list(Ulevel.variable[ref.variable])
    }
    if(type=="change"){
        out$args$ref.repetition <- list(U.time[ref.repetition])
        if(length(attr(time.var,"original"))>1){
            attr(out$args$ref.repetition,"original") <- attr(object$time$levels,"original")[ref.repetition,,drop=FALSE]
        }
    }

    add <- data.frame(matrix(NA, nrow = NROW(out$univariate)))
    names(add) <- variable
    add[] <- list(effect.variable)
    add <- cbind(add, n = effect.n)
    if(conditional.time){
        out$args$time <- time.var
        add <- cbind(add, stats::setNames(as.data.frame(effect.time), out$args$time))
    }
    if(n.strata > 1){
        out$args$strata <- paste(colnames(grid.conditional), collapse = ", ")
        add <- cbind(add, stats::setNames(as.data.frame(effect.strata), out$args$strata))
    }
    out$univariate <- cbind(add, out$univariate)
    rownames(out$univariate) <- rownames(M.effect)

    ## ** export
    attr(out,"class") <- append("effect_lmm",attr(out,"class"))
    return(out)
}

## * .effects_contrast
##' @description Generate contrast matrix
##' @noRd
.effects_contrast <- function(type, weight.type, grid.conditional, conditional.time, repetition,
                              U.strata, n.strata, 
                              U.time, n.time, nobs.time,
                              ref.repetition, prefix.time, sep.var){


    ## ** generate contrast matrix according to type across repetitions
    attr(U.time,"original") <- NULL
    if(type == "outcome"){
        M.contrast <- diag(1, n.time, n.time)
        dimnames(M.contrast) <- list(U.time, U.time)
        if(!is.null(repetition)){
            M.contrast <- M.contrast[repetition, ,drop=FALSE]
        }
    }else if(type == "change"){
        lsM.contrast <- lapply(1:length(ref.repetition), function(iRep){ ## iRep <- 3
            iRow <- rep(0, n.time)
            iRow[ref.repetition[iRep]] <- -1
            iM.contrast <- matrix(iRow, byrow = TRUE, nrow = n.time, ncol = n.time) + diag(1, n.time, n.time)
            dimnames(iM.contrast) <- list(U.time,U.time)
            if(is.null(repetition) && iRep == 1){
                iOut <- iM.contrast[-ref.repetition[iRep],,drop=FALSE]
            }else if(!is.null(repetition)){
                iOut <- iM.contrast[setdiff(repetition,U.time[ref.repetition[1:iRep]]),,drop=FALSE]
            }else{
                iOut <- iM.contrast[setdiff(U.time,U.time[ref.repetition[1:iRep]]),,drop=FALSE]
            }
            if(length(ref.repetition)>1 & NROW(iOut)>0){
                rownames(iOut) <- paste0(rownames(iOut),"-",ref.repetition[iRep])
            }
            return(iOut)
        })
        M.contrast <- do.call(rbind,lsM.contrast)        
    }else if(type %in% c("auc","auc-b")){
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
        if(is.null(repetition)){
            dtime.num <- diff(c(utils::head(time.num,1),time.num))/2 + diff(c(time.num,utils::tail(time.num,1)))/2
            if(type == "auc-b"){
                dtime.num[1] <- dtime.num[1] - sum(dtime.num)
            }
        }else{
            timeRep.num <- time.num[U.time %in% repetition]
            dtime.num <- rep(0, n.time)
            dtime.num[U.time %in% repetition] <- diff(c(utils::head(timeRep.num,1),timeRep.num))/2 + diff(c(timeRep.num,utils::tail(timeRep.num,1)))/2
            if(type == "auc-b"){
                dtime.num[which(U.time %in% repetition)] <- dtime.num[which(U.time %in% repetition)] - sum(dtime.num)
            }
        }
        M.contrast <- matrix(dtime.num, byrow = TRUE, nrow = 1, ncol = n.time)
        dimnames(M.contrast) <- list(NULL, U.time)
    }else if(type == "user"){
        if(is.null(colnames(weight.type))){
            M.contrast <- weight.type
            if(is.null(repetition)){
                colnames(M.contrast) <- U.time
            }else{
                colnames(M.contrast) <- repetition
            }
        }else{
            if(is.null(repetition)){
                M.contrast <- weight.type[,U.time,drop=FALSE]
            }else{
                M.contrast <- weight.type[,repetition,drop=FALSE]
            }
        }
    }

    ## ** possible stratification
    if(length(grid.conditional)>0){
 
        if(any(attr(U.strata,"table")==0)){
            ## in cross-over designs, there will be no observation at period drug=placebo under condition=active
            ## this case should be removed from the contrast matrix
            ls.contractC <- lapply(1:n.strata, function(iS){ ## iS <- 1
                iIndex <- which(attr(U.strata,"table")[,U.strata[iS]]>0)
                if(type == "outcome"){
                    iM <- M.contrast[iIndex,iIndex,drop=FALSE]
                    attr(iM,"time.row") <- U.time[iIndex]
                }else if(type == "change"){
                    if(1 %in% iIndex){
                        iM <- M.contrast[setdiff(iIndex,1)-1,iIndex,drop=FALSE]
                        attr(iM,"time.row") <- U.time[setdiff(iIndex,1)-1]
                    }else{
                        stop("Cannot evaluate the effect (no baseline measurement under condition ",U.strata[iS],"). \n")
                    }
                }else if(type %in% c("auc","auc-b","user")){
                    iM <- M.contrast[,iIndex,drop=FALSE]
                }
                attr(iM,"time.col") <- U.time[iIndex]
                attr(iM,"strata.col") <- rep(U.strata[iS],length(iIndex))
                attr(iM,"strata.row") <- rep(U.strata[iS],NROW(iM))
                return(iM)
            })
            out <- as.matrix(Matrix::bdiag(ls.contractC))
            attr(out,"time.col") <- unlist(lapply(ls.contractC,"attr","time.col")) ## wil be NULL in the case of auc,auc-b
            attr(out,"strata.col") <- unlist(lapply(ls.contractC,"attr","strata.col"))
            attr(out,"time.row") <- unlist(lapply(ls.contractC,"attr","time.row"))
            attr(out,"strata.row") <- unlist(lapply(ls.contractC,"attr","strata.row"))
        }else{
            out <- as.matrix(Matrix::bdiag(replicate(n = n.strata, M.contrast, simplify = FALSE)))
            attr(out,"time.col") <- rep(U.time, times = n.strata) ## wil be NULL in the case of auc,auc-b
            attr(out,"strata.col") <- unlist(lapply(U.strata, rep, times = NCOL(M.contrast)))
            attr(out,"time.row") <- rep(rownames(M.contrast), times = n.strata) ## wil be NULL in the case of auc,auc-b
            attr(out,"strata.row") <- unlist(lapply(U.strata, rep, times = NROW(M.contrast)))
        }
        colnames(out) <- paste0(attr(out,"time.col"), sep.var, attr(out,"strata.col"))

        if(type %in% c("outcome","change")){
            rownames(out) <- paste0(prefix.time,attr(out,"time.row"),sep.var,attr(out,"strata.row"))
        }else if(type %in%c("auc","auc-b","user")){
            if(nchar(prefix.time)==0){
                rownames(out) <- paste0(attr(out,"strata.row"))
            }else{
                rownames(out) <- paste0(prefix.time,sep.var,attr(out,"strata.row"))
            }
        }

        
    }else{
        out <- M.contrast
        attr(out,"time.col") <- U.time
        attr(out,"strata.col") <- rep(U.strata, times = NROW(out))
        attr(out,"time.row") <- rownames(M.contrast) ## will be NULL in the case of auc,auc-b
        attr(out,"strata.row") <- rep(U.strata, times = NROW(out))
        colnames(out) <- attr(out,"time.col")

        if(type %in% c("outcome","change")){
            rownames(out) <- paste0(prefix.time,attr(out,"time.row"))
        }else if(type %in% c("auc","auc-b")){
            if(nchar(prefix.time)==0){
                rownames(out) <- ""
            }else{
                rownames(out) <- paste0(prefix.time)
            }
        }else if(type == "user"){
            if(is.null(attr(out,"time.row"))){
                if(nchar(prefix.time)==0){
                    rownames(out) <- ""
                }else{
                    rownames(out) <- paste0(prefix.time)
                }
            }else{
                rownames(out) <- paste0(prefix.time,attr(out,"time.row"))
            }
        }

    }

    ## ** number of observations
    ## nobs.time is a matrix (n.rep,n.strata) where strata is defined by the conditional argument (number of unique covariate sets)

    if(length(grid.conditional)==0){ ## single strata
        attr(out,"n.col") <- nobs.time[,1]            
    }else{ ## multiple strata
        attr(out,"n.col") <- mapply(x = attr(out,"time.col"), y = attr(out,"strata.col"), FUN = function(x,y){nobs.time[x,y]})
    }
    
    attr(out,"n.row") <- apply(out, MARGIN = 1, FUN = function(iRow){ ## iRow <- out[1,]
        ## NOTE: attr(out,"time.col") %in% repetition to handle n.col calculated only on few timepoints, e.g. partial AUC.
        if(is.null(repetition)){
            iN <- unique(attr(out,"n.col")[iRow!=0])
        }else{
            iN <- unique(attr(out,"n.col")[iRow!=0 & attr(out,"time.col") %in% repetition])
        }
        if(length(iN)==1){return(iN)}else{return(NA)}
    })

    ## ** weight for averaging and new names
    if(conditional.time == FALSE && NROW(out)>n.strata){        
        if(type %in% c("outcome","change")){
            attr(out,"weight") <- do.call(rbind,lapply(U.strata, function(iStrata){ ## iStrata <- U.strata[1]
                iN <- attr(out,"n.row")*(attr(out,"strata.row")==iStrata)
                return(iN/sum(iN))
            }))
            rownames(attr(out,"weight")) <- U.strata

            if(nchar(prefix.time)==0){
                if(n.strata>1){
                    attr(out,"weight.name") <- stats::setNames(paste0(U.strata), U.strata)
                }else{
                    attr(out,"weight.name") <- stats::setNames("", U.strata)
                }                
            }else{
                if(n.strata>1){
                    attr(out,"weight.name") <- stats::setNames(paste0(prefix.time,sep.var,U.strata), U.strata)
                }else{
                    attr(out,"weight.name") <- stats::setNames(paste0(prefix.time), U.strata)
                }
            }
        } ## nothing to average over time for auc and auc-b
    }

    ## ** export
    return(out)
}
##----------------------------------------------------------------------
### effects.R ends here

