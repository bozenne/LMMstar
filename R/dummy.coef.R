### dummy.coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2024 (11:48) 
## Version: 
## Last-Updated: sep 30 2024 (14:15) 
##           By: Brice Ozenne
##     Update #: 70
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * dummy.coef.lmm (documentation)
##' @title Extract Mean Coefficients in Original Coding From a Linear Mixed Model
##' @description This expands the mean coefficients of the linear mixed model into one coefficient per level of the original variable,
##' i.e., including the reference level(s) where the fitted coefficients are 0.
##'
##' @param object  a \code{lmm} object.
##' @param use.na [logical] Should \code{NA} or \code{0} be used to represent excluded coefficients in singular models.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' dL$group <- as.factor(paste0("d",dL$X1))
##'
##' ## fit mixed model with interaction
##' eUN.lmm <- lmm(Y ~ visit*group, repetition =~visit|id, data = dL, df = FALSE)
##' dummy.coef(eUN.lmm)
##' 
##' ## fit mixed model with baseline constraint
##' dL$drug <- dL$group
##' dL[dL$visit==1,"drug"] <- "d0"
##' eUN.clmm <- lmm(Y ~ visit + visit:drug, repetition =~visit|id, data = dL, df = FALSE)
##' dummy.coef(eUN.clmm)
##' 
##' dL$drug <- factor(dL$group, levels = c("none","d0","d1"))
##' dL[dL$visit==1,"drug"] <- "none"
##' eUN.clmm.none <- lmm(Y ~ visit:drug, repetition =~visit|id, data = dL, df = FALSE)
##' dummy.coef(eUN.clmm.none)
##'

## * dummy.coef.lmm (code)
##' @export
dummy.coef.lmm <- function(object, use.na = FALSE, ...){

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    
    ## ** extract information
    ## mean coefficient
    e.mucoef <- stats::coef(object, effects = "mean")

    ## design matrix
    design.mean <- object$design$mean
    dataClasses <- attr(attr(object$design$mean, "terms"),"dataClasses")
    if(is.null(attr(design.mean, "M.level"))){
        design.mean <- .model.matrix_regularize(formula = stats::formula(object), data = stats::model.frame(object), augmodel = TRUE, type = "mean", drop.X = object$design$drop.X)
        term.labels <- attr(design.mean,"term.labels")
    }else{
        term.labels <- c("(Intercept)",attr(design.mean,"term.labels"))[attr(design.mean,"assign")+1]
    }
    M.level <- attr(design.mean, "M.level")
    M.term2Var <- attr(attr(design.mean,"formula"), "factors")
    intercept <- attr(attr(design.mean,"formula"), "intercept")

    ## level(s) associated to each variable
    xlevels  <- object$xfactor$mean
    ## numeric variables
    numeric.variable <- setdiff(colnames(M.level),names(xlevels))
    if(length(numeric.variable)>0){
        xlevels[numeric.variable] <- list(TRUE)
    }
    
    ## ** gather information
    Ulabel <- colnames(M.term2Var)

    ls.term <- lapply(Ulabel, function(iLabel){ ## iLabel <- colnames(M.term2Var)[3]
        iCol <- stats::setNames(M.term2Var[,iLabel], rownames(M.term2Var)) ## name may disappear when only a single column
        iVar <- names(which(iCol>0)) ## find active covariates for this regression term
        iGrid <- expand.grid(xlevels[intersect(iVar, names(xlevels))]) ## generate all combinations between covariate levels
        iM.level <- M.level[term.labels==iLabel,,drop=FALSE] ## find covariate levels corresponding to the regression term
        iM.level$value <- e.mucoef[rownames(iM.level)] ## find corresponding regression coefficient
        iGridA <- merge(x = iGrid, y = iM.level, all.x = TRUE, by = names(iGrid), sort = TRUE) ## merge with grid (may not be identical in case of constraint, e.g., no baseline difference)
        if(use.na==FALSE && any(is.na(iGridA$value))){           
            iGridA$value[is.na(iGridA$value)] <- 0 ## set missing coefficients to 0 instead of NA
        }
        iName <- iGridA[,colnames(iGridA) %in% setdiff(colnames(iGrid),numeric.variable),drop=FALSE]
        if(length(iName)>0){
            return(stats::setNames(iGridA$value, interaction(iName)))
        }else{
            return(unname(iGridA$value))
        }
    })
    names(ls.term) <- Ulabel

    if(intercept){
         ls.term <- c(list("(Intercept)" = e.mucoef["(Intercept)"]), ls.term)
    }

    ## ** export
    return(structure(ls.term, class = "dummy_coef", matrix = FALSE))
}

##----------------------------------------------------------------------
### dummy.coef.R ends here
