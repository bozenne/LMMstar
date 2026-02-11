### time.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 21 2025 (16:18) 
## Version: 
## Last-Updated: feb 10 2026 (15:52) 
##           By: Brice Ozenne
##     Update #: 26
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
    attr(out,"format") <- list(var = "orignal", cor = "original") ## format for the time variables (numeric/factor/original format) 
    return(out)
}

## * time.IND
##' @export
time.IND <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(var = "factor", cor = "original") ## format for the time variables (numeric/factor/original format) 
    return(out)
}

## * time.CS
##' @export
time.CS <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(var = "factor", cor = "factor") ## format for the time variables (numeric/factor/original format) 
    if(any(!is.na(x$name$cor[[1]]))){ ## add format of the sub-blocks
        format.cor <- attr(time(do.call(x$cross, args = list(formula=reformulate(termlabels="1")))), "format")
        attr(out,"format") <- mapply(attr(out,"format"), format.cor, FUN = union, SIMPLIFY = FALSE)
    }
    return(out)
}

## * time.EXP
##' @export
time.EXP <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(var = "factor", cor = "numeric") ## format for the time variables (numeric/factor/original format)
    if(any(!is.na(x$name$cor[[1]]))){ ## add format of the sub-blocks
        format.cor <- attr(time(do.call(x$cross, args = list(formula=reformulate(termlabels="1")))), "format")
        attr(out,"format") <- mapply(attr(out,"format"), format.cor, FUN = union, SIMPLIFY = FALSE)
    }
    return(out)
}

## * time.TOEPLITZ
##' @export
time.TOEPLITZ <- time.CS

## * time.UN
##' @export
time.UN <- time.CS

## * time.CUSTOM
##' @export
time.CUSTOM <- function(x, ...){
    out <- x$name$time
    attr(out,"format") <- list(var = "orignal", cor = "original") ## format for the time variables (numeric/factor/original format) 
    return(out)
}

##----------------------------------------------------------------------
### time.R ends here
