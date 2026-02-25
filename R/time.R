### time.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 21 2025 (16:18) 
## Version: 
## Last-Updated: feb 25 2026 (18:16) 
##           By: Brice Ozenne
##     Update #: 44
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * time.ID
##' @export
time.ID <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(variance = "original", correlation = "original") ## format for the time variables (numeric/factor/original format) 
    return(out)
}

## * time.IND
##' @export
time.IND <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(variance = "factor", correlation = "original") ## format for the time variables (numeric/factor/original format) 
    return(out)
}

## * time.CS
##' @export
time.CS <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(variance = "factor", correlation = "factor") ## format for the time variables (numeric/factor/original format)    
    if(!is.na(x$class["correlation.cross"])){ ## add format of the sub-blocks
        attr(out,"format")$correlation.cross <- attr(time(do.call(x$class["correlation.cross"], args = list(formula=reformulate(termlabels="1")))), "format")$correlation
    }
    return(out)
}

## * time.EXP
##' @export
time.EXP <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(variance = "factor", correlation = "numeric") ## format for the time variables (numeric/factor/original format)
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
time.UN <- time.CS

## * time.CUSTOM
##' @export
time.CUSTOM <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(variance = "original", correlation = "original") ## format for the time variables (numeric/factor/original format) 
    return(out)
}

##----------------------------------------------------------------------
### time.R ends here
