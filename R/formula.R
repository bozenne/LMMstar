### formula.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:53) 
## Version: 
## Last-Updated: jul  9 2024 (16:55) 
##           By: Brice Ozenne
##     Update #: 281
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * formula.lmm (code)
##' @export
formula.lmm <- function(x, effects = "mean", ...){

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    if(identical(effects,"all")){
        effects <- c("mean","variance")
    }
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)

    ## ** extract formula
    if("mean" %in% effects && "variance" %in% effects){
        return(list(mean = x$formula$mean,
                    variance = x$formula$var))
    }else if("mean" %in% effects){
        return(x$formula$mean)
    }else if("variance" %in% effects){
        return(x$formula$var)
    }
}

## * formula2var (examples)
##' @description Extract information from a formula (outcome, regressors, special terms)
##' @noRd
##' @examples
##'
##' formula2var(~X)
##' formula2var(Y~1)
##' formula2var(Y~0)
##' formula2var(Y~X)
##' formula2var(Y~X1*X2+s(age))
##' formula2var(Y+Z~1)
##' formula2var(~time|id)
##' formula2var(group~time|id)
##' formula2var(group~time1+time2|id)
##' formula2var(Y~(time|id))
##' formula2var(Y+Z~X+(1|id)+(1|region))
##' formula2var(Y+Z~s(X1)+X2*X3 + (X1|id:baseline))
##'
##' df <- cbind(Y=1, as.data.frame(matrix(1,ncol = 5, nrow = 1)))
##' formula2var(Y ~ ., data = df)

## * formula2var (code)
formula2var <- function(formula, data = NULL, specials = NULL, name.argument  = "formula",
                        suggestion = ""){

    ## ** normalize user input
    if(!inherits(formula,"formula")){
        stop("Incorrect type for argument \'",name.argument,"\': it must inherit from \"formula\". \n",
             suggestion)
    }
    if(!is.null(specials)){
        stop("Incorrect value for argument \'",name.argument,"\' \n",
             "Current version does not allow specials operator. \n")
    }

    ## ** process
    out <- list(formula = list(all = formula), ## note for regressor Y~X1+X2*X3
                vars = list(), ## X1, X2, X3, time, id
                terms = list(), ## X1, X2*X3, (time|id)
                index.terms = list() ## position of the term in the right hand side of the formula
                )                                   
    out$vars$all <- all.vars(formula)
    ff.terms <- stats::terms(formula, specials = specials, data = data)
    out$terms$intercept <- as.logical(attr(ff.terms,"intercept"))

    ## *** identify response variables
    ff.formulaLHS <- stats::update(formula, ~1)
    ff.varLHS <- all.vars(ff.formulaLHS)
    n.outcome <- length(ff.varLHS)

    if(n.outcome>0){
        ## any is necessary when the formula is long when desparse will be a vector with pieces of the formula
        if(any(grepl("(",deparse(ff.formulaLHS), fixed = TRUE))){
            stop("The left hand side of the formula should not contain any parenthesis. \n",
                 "Consider adding a transformed variable in the dataset instead of via the formula interface. \n")
        }
        out$vars$response <- ff.varLHS
    }

    ## *** restrict to the right hand side
    ff.termsRHS <- stats::delete.response(ff.terms)
    txt.termsRHS <- deparse(ff.termsRHS[[2]])
    ff.formulaRHS <- formula(ff.termsRHS)
    ff.varRHS <- all.vars(ff.termsRHS)

    ## *** identify regressor variables from variables involved in the structure/random effects based on |
    ## i.e. ~ X0 + X1*X2 + (time|id) ---> X0, X1, X2, time|id
    ff.gvarRHS <- as.character(attr(ff.termsRHS,"variables")[-1])
    test.special <- grepl(pattern = "|", x = ff.gvarRHS, fixed = TRUE)
    n.special <- sum(test.special)
    ff.varREG <- ff.gvarRHS[test.special==FALSE] 
    n.regressor <- length(ff.varREG)
    ff.varSPECIAL <- ff.gvarRHS[test.special]

    if(any(grepl(pattern = "||", x = ff.varSPECIAL, fixed = TRUE))){
        stop("Incorrect value for argument \'",name.argument,"\' \n",
             "Current version does not handle || for random effects. Use something like Y~X+(1|id). \n")
    }
    if(any(grepl(pattern = ":", x = ff.varSPECIAL, fixed = TRUE))){
        stop("Incorrect value for argument \'",name.argument,"\' \n",
             "Current version does not handle : for random effects. Use something like Y~X+(1|region/city). \n")
    }

    ## *** distinguish special (between repetition and ranef)
    ## repetition: single |, no regressor, e.g. ~time|id
    ## ranef: possibly several | wrapped with parentheses, e.g. ~ X + (1|id) + (1|region)
    split.termsRHS <- trimws(strsplit(txt.termsRHS, split = "+", fixed = TRUE)[[1]]) ## may not be correct with repetition, e.g. ~ time1+time2|id (is fixed below)

    if(n.special == 0){
        type.special <- "none"
    }else if(length(ff.varSPECIAL)==1 && ff.varSPECIAL==txt.termsRHS){
        ## parenthesis are removed by terms so equality indicates no parentheses
        type.special <- "repetition"
    }else if(all(paste0("(",ff.varSPECIAL,")") %in% split.termsRHS)){
        type.special <- "ranef"
    }else{
        test.start <- grepl(split.termsRHS, pattern = "(", fixed = TRUE)
        test.stop <- grepl(split.termsRHS, pattern = ")", fixed = TRUE)

        if(any((test.start==TRUE)*(test.stop==FALSE))){
            stop("Incorrect value for argument \'",name.argument,"\' \n",
                 "Cannot handle random effects with several covariates (symbol + within parentheses) \n")
        }else{
            stop("Incorrect value for argument \'",name.argument,"\' \n",
                 "Could not recognize special operators. \n")
        }
    }

    ## *** reconstruct terms
    if(type.special == "repetition"){
        out$terms$all <- txt.termsRHS
    }else{
        out$terms$all <- split.termsRHS
    }

    ## *** reconstruct regressor and differentiate repetition from random effects
    if((n.regressor == 0)  && (n.outcome > 0) && (type.special != "repetition")){

        if(attr(ff.terms,"intercept")==1){
            out$formula$regressor <- ff.formulaLHS
            out$formula$design <- ~1
        }else{ ## NOTE: when updating formula not using stats::drop.terms or stats::udpate as it re-writes the interaction X1*X2 ---> X1 + X2 + X1:X2
            out$formula$regressor <- stats::as.formula(paste0(deparse(ff.formulaLHS[[2]]),"~0"))
            out$formula$design <- ~0
        }

    }else if(type.special == "none"){

        out$index.terms$regressor <- 1:n.regressor
        out$terms$regressor <- ff.varREG
        out$vars$regressor <- ff.varRHS
        out$formula$regressor <- formula        
        out$formula$design <- ff.formulaRHS        

    }else if(type.special == "ranef"){
        out$index.terms$regressor <- which(split.termsRHS %in% paste0("(",ff.varSPECIAL,")") == FALSE)
        out$terms$regressor <- out$terms$all[out$index.terms$regressor]
        out$formula$regressor <- updateFormula(formula, drop.x = out$terms$all[-out$index.terms$regressor])
        if(n.outcome>0){
            out$formula$design <- stats::as.formula(paste("~",deparse(out$formula$regressor[[3]]))) ## remove outcome
        }
        out$vars$regressor <- all.vars(out$formula$regressor)        
    }

    ## *** extract time and cluster variables from special
    if(type.special=="repetition"){

        out$vars$repetition <- ff.varRHS
        out$terms$repetition <- out$terms$all
        out$index.terms$repetition <- 1
        out$formula$repetition <- ff.formulaRHS

        
        if(countChar(ff.varSPECIAL, pattern = "|", fixed = TRUE)!=1){
            stop("Incorrect value for argument \'",name.argument,"\' \n",
                 "Right hand side of the formula should only contain one term, e.g. ~time|cluster.")
        }

        ## time
        var.time <- all.vars(ff.formulaRHS[[2]][[2]])
        if(length(var.time)>0){
            out$vars$time <- var.time
        }
        ## cluster
        out$vars$cluster <- all.vars(ff.formulaRHS[[2]][[3]])
        
    }else if(type.special=="ranef"){

        out$ranef <- list(formula = updateFormula(ff.formulaRHS, drop.x = out$terms$regressor))
        out$ranef$vars <- all.vars(out$ranef$formula)
        out$ranef$terms <- setdiff(out$terms$all[out$terms$all %in% out$terms$regressor == FALSE], c("0","1")) ## drop intercept
        out$ranef$index.terms <- which(out$terms$all %in% out$terms$regressor == FALSE)

        ls.timeCluster <- lapply(out$ranef$term, function(iRanef){ ## iRanef <- out$terms$ranef[[1]]
            iForm <- stats::as.formula(paste0("~",iRanef))
            iForm2 <- iForm[[2]][[2]] ## select after ~ and inside ()
            iOut <- list(time = all.vars(iForm2[[2]]),
                         hierarchy = all.vars(iForm2[[3]]),
                         cluster = all.vars(iForm2[[3]])[1])
            if(length(iOut$time)==0){
                iOut$name <- "(Intercept)"
            }else if(length(iOut$time)==1){
                iOut$name <- iOut$time
            }else{
                stop("Cannot handle multiple covariates in a single random effect. \n")
            }
            return(iOut)
        })
        var.time <- unique(unlist(lapply(ls.timeCluster,"[[","time")))
        if(length(var.time)>0){
            out$vars$time <- var.time
        }
        out$ranef$cluster <- unique(unlist(lapply(ls.timeCluster,"[[","cluster")))
        out$ranef$hierarchy <- stats::setNames(lapply(ls.timeCluster,"[[","hierarchy"),
                                               sapply(ls.timeCluster,"[[","name"))
        out$ranef$crossed <- length(out$ranef$cluster)>1
        out$ranef$nested <- any(lengths(out$ranef$hierarchy)>1)
    }

    ## ** export
    out$special <- type.special
    return(out)
}

## * formula2repetition
##' @description Normalize the formula and repetition arguments.
##' In particular if repetition is missing but formula contains a cluster variable, it will try to re-create the repetition variable
##'
##' @param formula [formula]  Y ~ X1 + X2 
##' @param data [data.frame]
##' @param repetition [formula] ~ time | cluster
##' @param keep.time [logical] if formula = Y ~ time|cluster should the new formula be Y ~ time or Y ~ 1
##' 
##' @noRd
##'
##' @examples
##' 
formula2repetition <- function(formula, data, repetition, keep.time, filter){

    ## ** formula
    if(!inherits(formula,"formula")){
        stop("Argument \'formula\' should inherit from formula. \n")
    }
    detail.formula <- formula2var(formula, data = data)
    if(any(setdiff(detail.formula$var$all,".") %in% names(data) == FALSE)){
        stop("Argument \'formula\' incompatible with argument \'data\': could not find column(s) \"",paste(setdiff(detail.formula$var$all,names(data)), collapse ="\", \""),"\".\n")
    }
    if(detail.formula$special=="ranef"){
        stop("Argument \'formula\' should not contain any random effect. \n",
             "Should be something like Y ~ time or Y ~ time + group, i.e., without ",detail.formula$ranef$terms[1],". \n")
    }
    if(detail.formula$special=="repetition"){
        if(missing(repetition) || is.null(repetition)){

            if(length(detail.formula$var$time)==0){
                stop("Argument \'formula\' contain a cluster but no repetition variable. \n",
                     "Should be something like Y ~ time or Y ~ time + group, i.e., without |",detail.formula$var$cluster[1],". \n",
                     "Use the argument \'repetition\' to define the repetition and cluster variables. \n")
            }
            terms.formula <- stats::terms(formula)
            
            ## check variables that are constant within individuals: should be excluded from the time variables
            Mduplicated <- do.call(rbind,by(data[detail.formula$var$time], data[detail.formula$var$cluster], function(iDF){
                apply(iDF, MARGIN = 2, FUN = function(iX){sum(!duplicated(iX))})
            }, simplify = FALSE))
            
            test.duplicated <- colSums(Mduplicated!=1)
            if(all(test.duplicated==0)){
                stop("Argument \'formula\' contain a cluster but all possible repetition variables are constant within cluster. \n",
                     "Should be something like Y ~ time or Y ~ time + group, i.e., without |",detail.formula$var$cluster[1],". \n",
                     "Use the argument \'repetition\' to define the repetition and cluster variables. \n")
            }else if(any(test.duplicated==0)){
                term.rm <- names(test.duplicated)[test.duplicated==0]
                repetition <- stats::reformulate(termlabels = paste(paste(setdiff(detail.formula$vars$time,term.rm),collapse="+"),"|",detail.formula$vars$cluster))
            }else{
                term.rm <- NULL
                repetition <- stats::reformulate(attr(stats::delete.response(terms.formula),"term.labels"))
            }
            
            if(keep.time){
                formula <- stats::reformulate(termlabels = detail.formula$vars$time, response = paste(detail.formula$vars$response, collapse = " + "))
            }else{
                formula <- stats::reformulate(termlabels = "1", response = paste(detail.formula$vars$response, collapse = " + "))
            }
            detail.formula <- formula2var(formula, data = data)
        }else{
            stop("Argument \'formula\' should not contain any grouping variable. \n",
                 "The argument \'repetition\' should be used instead. \n")
        }
    }

    ## ** repetition
    if(missing(repetition) || !is.null(repetition)){
        if(!inherits(repetition,"formula")){
            stop("Argument \'repetition\' should inherit from formula. \n")
        }
        detail.repetition <- formula2var(repetition, data = data)
        if(any(detail.repetition$var$all %in% names(data) == FALSE)){
            stop("Argument \'repetition\' incompatible with argument \'data\': could not find column(s) \"",paste(setdiff(detail.repetition$var$all,names(data)), collapse ="\", \""),"\".\n")
        }
        if(detail.repetition$special!="repetition"){
            stop("Could not identifying the cluster variable in argument \'repetition\'. \n",
                 "Should be a formula such as Y~time|id. \n")
        }
        if(length(detail.repetition$vars$cluster)==0){
            stop("Could not identifying the cluster variable in argument \'repetition\'. \n",
                 "Should be a formula such as Y~time|id. \n")
        }
        if(length(detail.repetition$vars$cluster)>1){
            stop("Cannot handle multiple cluster variable in argument \'repetition\'. \n",
                 "Should be a formula such as Y~time|id. \n")
        }
        if(length(detail.repetition$vars$time)==0){
            stop("Could not identifying the time variable in argument \'repetition\'. \n",
                 "Should be a formula such as Y~time|id. \n")
        }

        if(length(detail.repetition$vars$time)==1){
            if(any(unlist(tapply(data[[detail.repetition$vars$time]],data[[detail.repetition$vars$cluster]],duplicated, simplify = FALSE)))){
                stop("Mismatch between argument \'repetition\' and argument \'data\'. \n",
                     "There should no duplicated \"",detail.repetition$vars$time,"\" value (timepoint) within observations with the same \"",detail.repetition$vars$cluster,"\" value (cluster). \n")
            }
        }else{
            if(any(unlist(tapply(interaction(data[,detail.repetition$vars$time,drop=FALSE],drop=TRUE),data[[detail.repetition$vars$cluster]],duplicated, simplify = FALSE)))){
                stop("Mismatch between argument \'repetition\' and argument \'data\'. \n",
                     "There should no duplicated \"",detail.repetition$vars$time,"\" value (timepoint) within observations with the same \"",detail.repetition$vars$cluster,"\" value (cluster). \n")
            }
        }

    }else{
        detail.repetition <- NULL
    }

    ## ** handle dots (.~Group or Y~.)
    browser()
    if(identical(detail.formula$vars$response,".") & identical(detail.formula$vars$regressor,".")){
        stop("Argument \'formula\' cannot be .~. as the left or right hand side need to be explicit \n",
             "Consider for instance using .~1. \n",sep = "")
    }else if(identical(detail.formula$vars$response,".")){
        name.X <- detail.formula$var$regressor
        name.Y <- setdiff(names(data),c(name.X,detail.repetition$var$all))
        if(length(name.Y)==0){
            stop("Incorrect argument \'formula\': no variable on the left hand side of the formula. \n")
        }
        termlabels <- ifelse(is.null(detail.formula$vars$regressor),"1",paste(detail.formula$vars$regressor,collapse="+"))
        formula <- stats::reformulate(response = paste0(name.Y,collapse="+"), termlabels = termlabels)
        detail.formula <- formula2var(formula)$detail.formula
    }else if(identical(detail.formula$vars$regressor,".")){
        name.Y <- detail.formula$var$response
        formula <- stats::formula(stats::terms(formula, data = data))
        detail.formula <- formula2var(formula)
    }
    
    ## ** return
    return(list(formula = formula,
                detail.formula = detail.formula,
                repetition = repetition,
                detail.repetition = detail.repetition))

}

## * updateFormula
##' @description Remove or add a term in a formula while keeping interaction term untouched
##' @noRd
##' @details when updating formula not using stats::drop.terms or stats::udpate as it re-write the interaction X1*X2 ---> X1 + X2 + X1:X2
##' @examples
##' ## issue 
##' stats::drop.terms(terms(Y~X1*X2+(1|id)), dropx = 3)
##' stats::update(terms(Y~X1*X2+(1|id)), .~.-(1|id))
##'
##' ## solution
##' updateFormula(Y~X1*X2+(1|id), drop.x = "(1|id)")
##' updateFormula(Y~0+X1*X2+(1|id), drop.x = "(1|id)")
##' updateFormula(Y~X1*X2, add.x = "(1 | id)")
##'
##' ## other cases
##' updateFormula(Y~X1+X2+X3, drop.x = "X1")
##' updateFormula(Y~X1*X2+X3, drop.x = "X1") ## WARNING DO NOT DROP THE INTERACTION
##' updateFormula(Y~X1*X2+X3, drop.x = "X1*X2")
##' 
##' updateFormula(Y~X1, drop.x="X1")
##' updateFormula(Y~X1, drop.y = TRUE, drop.x="X1")
##' 
updateFormula <- function(formula, add.x = NULL, drop.x = NULL, add.y = NULL, drop.y = FALSE){

    drop.x <- gsub(" ","",drop.x)
    if(!inherits(formula,"formula") || (length(formula) %in% 2:3 == FALSE)){
        stop("Argument \'formula\' should be a formula with length 2 or 3. \n")
    }
    test.response <- length(formula) == 3

    txt.formula <- as.character(utils::tail(formula,1))
    term.formula <- gsub(" ","",strsplit(txt.formula, split = "+", fixed = TRUE)[[1]])
    ## if(any(drop.x %in% term.formula == FALSE)){
    ##     stop("Mismatch between argument \'formula\' and \'drop.x\', \n",
    ##          "Could not find \"",paste(drop.x[drop.x %in% term.formula == FALSE], collapse = "\" \""),"\". \n")
    ## }
    if(!is.null(drop.x)){
        term.formula <- term.formula[term.formula %in% drop.x == FALSE]
    }
    if(!is.null(add.x)){
        term.formula <- c(term.formula, add.x)
    }
    if(length(term.formula)==0){
        term.formula <- 1
    }
    if(!test.response || drop.y){
        txt.new <- paste0(deparse(formula[[1]]),paste0(term.formula, collapse="+"))
    }else{
        txt.new <- paste0(deparse(formula[[2]]),deparse(formula[[1]]),paste0(term.formula, collapse="+"))
    }
    if(!is.null(add.y) && length(add.y)>0){
        txt.new <- paste0(paste0(add.y,collapse="+"),txt.new)
    }
    return(stats::as.formula(txt.new))
}

##----------------------------------------------------------------------
### formula.R ends here
