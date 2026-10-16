### levels.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 18 2021 (10:32) 
## Version: 
## Last-Updated: okt 16 2026 (14:41) 
##           By: Brice Ozenne
##     Update #: 21
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Contrasts and Reference Level in a Linear Mixed Model
##' @description Contrasts and reference level used when modeling the mean in a linear mixed modek.
##' 
##' @param x an \code{lmm} object.
##' 
##' @return a list with two elements \itemize{
##' \item all: contrast matrix for each categorical or factor variable
##' \item reference: reference level: one value for each categorical variable
##' }
##'
##' @keywords methods
##' 
##' @export
levels.lmm <- function(x){

    x.formula <- formula(x, effects = "mean")
    
    data.X <- x$data[all.vars(stats::delete.response(x.formula))]
    C <- lapply(data.X, function(iCol){
        if(inherits(iCol,"factor")){stats::contrasts(iCol)}else if(inherits(iCol,"character")){stats::contrasts(as.factor(iCol))}
    })
    if(length(C)>0){
        C <- C[!unlist(lapply(C, is.null))]
        if(length(C)>0 && "(Intercept)" %in% x$design$param$name){
            ref.level <- unlist(lapply(names(C), function(iC){
                stats::setNames(rownames(C[[iC]])[1],iC)
            }))  
        }else{
            ref.level <- NULL
        }
    }else{
        ref.level <- NULL
    }
    return(list(all = C, reference = ref.level))
}

##----------------------------------------------------------------------
### levels.R ends here
