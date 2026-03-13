### time.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 21 2025 (16:18) 
## Version: 
## Last-Updated: mar 13 2026 (14:48) 
##           By: Brice Ozenne
##     Update #: 52
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * time.ID
## Format for the covariates and the time variable
## factor: covariate and time as factor
## factor0: covariates as factor
## numeric: covariate as factor, time as numeric
## original: nothing

##' @export
time.ID <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(variance = "original") ## format for the time variables 
    return(out)
}

## * time.IND
##' @export
time.IND <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(variance = "factor") ## format for the time variables
    return(out)
}

## * time.CS
##' @export
time.CS <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(variance = "factor", correlation = "factor0") ## format for the time variables
    if(!is.na(x$class["correlation.cross"])){ ## add format of the sub-blocks
        attr(out,"format")$correlation.cross <- attr(time(do.call(x$class["correlation.cross"], args = list(formula=reformulate(termlabels="1")))), "format")$correlation
    }
    return(out)
}

## * time.EXP
##' @export
time.EXP <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(variance = "factor", correlation = "numeric") ## format for the time variables
    if(!is.na(x$class["correlation.cross"])){ ## add format of the sub-blocks
        attr(out,"format")$correlation.cross <- attr(time(do.call(x$class["correlation.cross"], args = list(formula=reformulate(termlabels="1")))), "format")$correlation
    }
    return(out)
}

## * time.TOEPLITZ
##' @export
time.TOEPLITZ <- time.EXP

## * time.UN
##' @export
time.UN <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(variance = "factor", correlation = "factor") ## format for the time variables 
    if(!is.na(x$class["correlation.cross"])){ ## add format of the sub-blocks
        attr(out,"format")$correlation.cross <- attr(time(do.call(x$class["correlation.cross"], args = list(formula=reformulate(termlabels="1")))), "format")$correlation
    }
    return(out)
}

## * time.DUN
##' @export
time.DUN <- time.UN

## * time.AR1
##' @export
time.AR1 <- time.UN

## * time.CUSTOM
##' @export
time.CUSTOM <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(variance = "original", correlation = "original") ## format for the time variables
    return(out)
}

##----------------------------------------------------------------------
### time.R ends here
