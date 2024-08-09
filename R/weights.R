### weights.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct 11 2022 (10:56) 
## Version: 
## Last-Updated: aug  8 2024 (15:26) 
##           By: Brice Ozenne
##     Update #: 46
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * weights.rbindWald_lmm (documentation)
##' @title Extract Weights Used to Pool Estimates
##' @description Extract weights used to pool estimates.
##' 
##' @param object a \code{Wald_lmm} object, output of \code{anova.lmm}, or \code{rbind.lmm}, or \code{mlmm}.
##' @param method [character] method for combining the estimates, one of \code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.gls1"}, \code{"pool.rubin"}.
##' @param effects [character] should the weights relative to the Wald test be output (\code{"Wald"}),
##' or the one relative to the linear mixed model parameters (\code{"all"})?
##' @param transform.names [logical] should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return a numeric vector whose elements sum to 1.
##'
##' @keywords methods
##'
##' @examples
##' set.seed(10)
##' dL <- sampleRem(250, n.times = 3, format = "long")
##'
##' e.mlmm <- mlmm(Y~X1+X2+X6, repetition = ~visit|id, data = dL,
##'                by = "X4", effects = "X1=0", structure = "CS")
##' weights(e.mlmm, method = c("average","pool.se","pool.gls"))

## * weights.mlmm (code)
##' @export
weights.rbindWald_lmm <- function(object, method, effects = "Wald", transform.names = TRUE, ...){

    ## ** check and normalize user input

    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    pool.method <- setdiff(options$pool.method,"p.rejection")


    ## *** method
    if(missing(method)){
        stop("Argument \'method\' is missing.\n",
             "Should be one of \"",paste(pool.method, collapse = "\" \""),"\".\n")
    }
    if(any(method %in% pool.method == FALSE)){
        stop("Argument \'method\' is missing.\n",
             "Should be one of \"",paste(pool.method, collapse = "\" \""),"\".\n")
    }

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character. \n")
    }
    valid.effects <- c("Wald","all")
    if(any(effects %in% valid.effects == FALSE)){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }
    if(length(effects)!=1){
        stop("Argument \'effects\' should have length 1. \n")
    }

    
    ## ** Extract contrast
    object.confint <- stats::confint(object, columns = "estimate", method = method, df = FALSE, options = options)
    if(effects == "all"){
        out <- attr(object.confint,"contrast") %*% stats::model.tables(object, effects = "contrast", transform.names = transform.names, simplify = FALSE)[[1]]
    }else{
        out <- attr(object.confint,"contrast")
    }
    
    ## ** Export
    return(out)

}

##----------------------------------------------------------------------
### weights.R ends here
