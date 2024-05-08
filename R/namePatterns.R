### namePatterns.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 11 2024 (09:49) 
## Version: 
## Last-Updated: maj  7 2024 (16:24) 
##           By: Brice Ozenne
##     Update #: 77
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Group Covariance Patterns
##' @description Group covariance patterns 
##'
##' @param structure covariance pattern.
##' @param xfactor [list] list with the levels of each categorical variable for the variance and correlation structure.
##' @param ignore.time [logical] should the time variable be ignored when forming the groups?
##' @param sep [categorical vector of length 2] characters used to aggregate variables and levels in the pattern name.
##' 
##' @noRd

`.nameUpatterns` <-function(structure, xfactor, ignore.time, sep){
    UseMethod(".nameUpatterns")
}

## * .nameUpatterns.ID
.nameUpatterns.ID <- function(structure, xfactor, ignore.time, sep){
    if(all(is.na(structure$name$strata))){
        return(structure$Upattern$index.strata)
    }else{
        return(xfactor$var[[structure$name$strata]][structure$Upattern$index.strata])
    }

}

## * .nameUpatterns.IND
.nameUpatterns.IND <- function(structure, xfactor, ignore.time = FALSE, sep){

        
    ## ** find groups
    Upattern <- structure$Upattern
    out <- groupSet(Upattern$param, strata = Upattern$strata)
    if(length(unique(out))==1){
        return(out)
    }

    ## ** name groups
    pattern2lp <- structure$var$pattern2lp
    lp2data <- structure$var$lp2data
    if(ignore.time){
        lp2data <- lp2data[names(lp2data) %in% attr(structure$name$time,"original") == FALSE]
        if(NCOL(lp2data)==0){
            return(rep(1,length(out)))
        }
    }
    
    vec.name <- try(by(Upattern, out, function(iUpattern){ ## iUpattern <- Upattern[out == 1,]
        iUlp.var <- unique(unlist(pattern2lp[unique(iUpattern$var)]))
        iData <- lp2data[iUlp.var,,drop=FALSE]
        iLs.name <- lapply(colnames(iData), function(iCol){
            if(all(xfactor$var[[iCol]] %in% iData[[iCol]])){
                return(NULL)
            }else{
                paste(xfactor$var[[iCol]][xfactor$var[[iCol]] %in% iData[[iCol]]], collapse = sep[2])
            }
        })
        return(paste(unlist(iLs.name), collapse = sep[1]))
    }, simplify = FALSE), silent = FALSE)

    ## ** udpate
    if(!inherits(vec.name, "try-error") && all(lengths(vec.name)==1)){
        if(ignore.time){
            out <- factor(out, labels = unlist(vec.name))
        }else{
            out <- factor(out, labels = make.unique(unlist(vec.name)))
        }
        
    }

    ## ** export
    return(out)
}

## * .nameUpatterns.CS
.nameUpatterns.CS <- function(structure, xfactor, ignore.time = FALSE, sep){

    ## ** find groups
    Upattern <- structure$Upattern
    out <- groupSet(Upattern$param, strata = Upattern$strata)
    if(length(unique(out))==1){
        return(out)
    }

    ## ** name groups
    pattern2lp.var <- structure$var$pattern2lp
    lp2data.var <- structure$var$lp2data
    
    pattern2lp.cor <- structure$cor$pattern2lp
    lp2data.cor <- structure$cor$lp2data
    if(ignore.time){
        lp2data.var <- lp2data.var[names(lp2data.var) %in% attr(structure$name$time,"original") == FALSE]
        lp2data.cor <- lp2data.var[names(lp2data.cor) %in% attr(structure$name$time,"original") == FALSE]
        if(NCOL(lp2data.var)==0 && NCOL(lp2data.cor)==0){
            return(rep(1,length(out)))
        }
    }
    
    vec.name <- try(by(Upattern, out, function(iUpattern){ ## iUpattern <- Upattern[out == 1,]
        iUlp.var <- unique(unlist(pattern2lp.var[unique(iUpattern$var)]))
        iData.var <- lp2data.var[iUlp.var,,drop=FALSE]
        iName.var <- lapply(colnames(iData.var), function(iCol){
            if(all(xfactor$var[[iCol]] %in% iData.var[[iCol]])){
                return(NULL)
            }else{
                paste(xfactor$var[[iCol]][xfactor$var[[iCol]] %in% iData.var[[iCol]]], collapse = sep[2])
            }
        })        
        iUlp.cor <- unique(unlist(pattern2lp.cor[unique(iUpattern$cor)]))
        iData.cor <- lp2data.cor[iUlp.cor,,drop=FALSE]
        iName.cor <- lapply(colnames(iData.cor), function(iCol){ ## iCol <- colnames(iData.cor)
            if(iCol %in% names(xfactor$cor)){
                if(all(xfactor$cor[[iCol]] %in% iData.cor[[iCol]])){
                    return(NULL)
                }else{
                    paste(xfactor$cor[[iCol]][xfactor$cor[[iCol]] %in% iData.cor[[iCol]]], collapse = sep[2])
                }
            }else if(length(iData.cor)==1){
                return(unlist(iData.cor))
            }else{
                return(NULL)
            }
        })
        
        return(paste(unlist(union(iName.var,iName.cor)), collapse = sep[1]))
    }, simplify = FALSE), silent = FALSE)

    ## ** udpate
    if(!inherits(vec.name, "try-error") && all(lengths(vec.name)==1)){
        if(ignore.time){
            out <- factor(out, labels = unlist(vec.name))
        }else{
            out <- factor(out, labels = make.unique(unlist(vec.name)))
        }
    }
    
    ## ** export
    return(out)
}

## * .nameUpatterns.UN
.nameUpatterns.UN <- .nameUpatterns.IND

## * .nameUpatterns.RE
.nameUpatterns.RE <- .nameUpatterns.IND

## * .findUpatterns.TOEPLITZ
.nameUpatterns.TOEPLITZ <- .nameUpatterns.IND

## * .nameUpatterns.EXP
.nameUpatterns.EXP <- .nameUpatterns.IND

## * .findUpatterns_CUSTOM
.nameUpatterns.CUSTOM <- .nameUpatterns.CS

##----------------------------------------------------------------------
### namePatterns.R ends here
