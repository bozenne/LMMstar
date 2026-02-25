### structure-crossBlock.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 10 2026 (15:39) 
## Version: 
## Last-Updated: feb 25 2026 (18:09) 
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

## * .crossBlockID
.crossBlockID <- function(structure, XpairPattern, U.strata, name.cov, sep){

    ## ** run 'normal' CS
    out <- .crossBlockCS(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, sep = sep)

    ## ** add constraint
    lp2data <- structure$cor$lp2data
    lpB <- XpairPattern$lpB

    diffS.nesting <- rowSums(lp2data[lpB[,1],name.cov,drop=FALSE]==lp2data[lpB[,2],name.cov,drop=FALSE], as.factor = FALSE)
    out$constraint <- unname(ifelse(diffS.nesting==0,0,NA))

    ## ** export
    return(out)

}


## * .crossBlockCS
.crossBlockCS <- function(structure, XpairPattern, U.strata, name.cov, sep){
    
    ## ** extract data
    lpB <- XpairPattern$lpB
    strataB <- XpairPattern$strataB
    lp2data <- structure$cor$lp2data
    name.cor <- unlist(structure$name$cor)
    name.strata <-  structure$name$strata
    name.cluster <- attr(structure$name$cluster,"original")
    n.strata <- length(U.strata)
    
    ## ** identify parameters
    if(structure$twin){ ## same correlation for all blocks at the same nesting level
        diff.nesting <- (lp2data[lpB[,1],name.cov,drop=FALSE]==lp2data[lpB[,2],name.cov,drop=FALSE])
        diffS.nesting <- collapse(diff.nesting, sep = sep["lp"], as.factor = FALSE)
        index.rhoB <- paste(diffS.nesting, strataB, sep = sep["rho.strata"])
    }else{ ## correlation specific to each pair of linear predictors
        index.rhoB <- paste(collapse(lp2data[lpB[,1],,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                            collapse(lp2data[lpB[,2],,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                            strataB, sep = sep["rho.strata"])
    }
    ## internal identifier for the parameter        
    code.rhoB <- tapply(paste("B",lpB[,1],lpB[,2], sep=sep["rho.name"]), index.rhoB, FUN = identity, simplify = FALSE)
    ## strata
    strata.rhoB <- tapply(strataB, index.rhoB, unique)
    ## provide a name for the parameter
    if(structure$twin){
        ## order the covariates to find the ordering of the nesting
        diffE.nesting <- cbind(diff.nesting, matrix(TRUE, ncol = 1, nrow = NROW(diff.nesting), dimnames = list(NULL, name.cluster)))
        name.nesting <- names(sort(colSums(diffE.nesting), decreasing = FALSE)) ## outer nesting is when fewest pairs have equal levels
        ## sanity check
        lsDiff.nesting <- apply(cbind(FALSE, diffE.nesting[,name.nesting,drop=FALSE]),1,diff, simplify = FALSE)
        if(any(unlist(lsDiff.nesting)<0)){
            ## if one specifies id/baseline/hour and that some observations from the same subject have same baseline but different hours
            ## while others have same hours and different baseline
            ## hour should be changed to be relative to baseline, e.g. (no,0h) (no,1h) (yes,0h) (yes,1h) ----> (no,no.0h) (no,no.1h) (yes,no.0h) (yes,no.1h)
            warning("Variables \"",paste(name.nesting,collapse="\", \""),"\" do not define a proper nesting. \n")
        }
        ## naming
        level.rhoB <- paste0("(",sapply(lsDiff.nesting, function(iVec){paste(rev(names(iVec[iVec>0])),collapse="/")}),")")
        if(n.strata>1){
            level.rhoB <- paste(level.rhoB, U.strata[strata.rhoB], sep = sep["rho.strata"])
        }
    }else{ 
        level.rhoB <- paste0("(",
                             nlme::collapse(lp2data[lpB[,1],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                             ",",
                             nlme::collapse(lp2data[lpB[,2],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                             ")")
        name.nesting <- NULL
    }


    ## ** constraint
    constraint.rhoB <- rep(NA, length(code.rhoB))

    ## ** export
    return(list(code = code.rhoB,
                level = level.rhoB,
                strata = unname(strata.rhoB),
                nesting = name.nesting,
                constraint = constraint.rhoB))
}

## * .crossBlockTOEPLITZ
.crossBlockTOEPLITZ <- function(structure, XpairPattern, U.strata, name.cov, sep){

    ## ** extract data
    lpB <- XpairPattern$lpB
    strataB <- XpairPattern$strataB
    lp2data <- structure$cor$lp2data.cross
    name.strata <-  structure$name$strata
    name.time <-  structure$name$time
    n.strata <- length(U.strata)

    ## ** identify parameters
    diff.time <- abs(lp2data[lpB[,1],name.time]-lp2data[lpB[,2],name.time])
    
    if(structure$twin){ ## shared correlation parameters for all blocks
        index.rhoB <- paste(diff.time, strataB, sep = sep["rho.strata"])
    }else{ ## specific correlation parameters for each block
        index.rhoB <- paste(diff.time,
                            collapse(lp2data[lpB[,1],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                            collapse(lp2data[lpB[,2],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                            strataB, sep = sep["rho.strata"])
    }
    ## internal identifier for the parameter
    code.rhoB <- tapply(paste("B",lpB[,1],lpB[,2], sep=sep["rho.name"]), index.rhoB, FUN = identity, simplify = FALSE)
    ## strata
    strata.rhoB <- tapply(strataB, index.rhoB, unique)
    ## provide a name for the parameter
    if(structure$twin){
        level.rhoB <- paste0("(",tapply(diff.time, index.rhoB, unique),")")
    }else{
        level.rhoB <- tapply(paste0("(",diff.time,"|",collapse(lp2data[lpB[,1],name.cov,drop=FALSE], sep = sep["lp"]),",",collapse(lp2data[lpB[,2],name.cov,drop=FALSE], sep = sep["lp"]),")"),
                             index.rhoB, unique)
    }
    if(n.strata>1){
        level.rhoB <- paste(strata.rhoB, U.strata[strata.rhoB], sep = sep["rho.strata"])
    }
    
    ## ** constraint
    constraint.rhoB <- rep(NA, length(code.rhoB))

    ## ** export
    return(list(code = unname(code.rhoB),
                level = unname(level.rhoB),
                strata = unname(strata.rhoB),
                nesting = NULL,
                constraint = constraint.rhoB))

}

## * .crossBlockDUN
.crossBlockDUN <- function(structure, XpairPattern, U.strata, name.cov, sep){

    ## ** extract data
    lpB <- XpairPattern$lpB
    strataB <- XpairPattern$strataB
    lp2data <- structure$cor$lp2data.cross
    name.strata <-  structure$name$strata
    name.time <-  structure$name$time
    n.strata <- length(U.strata)

    ## ** identify parameters
    diff.time <- abs(lp2data[lpB[,1],name.time]-lp2data[lpB[,2],name.time])
    
    if(structure$twin){ ## shared correlation parameters for all blocks
        index.rhoB <- paste(diff.time, strataB, sep = sep["rho.strata"])
    }else{ ## specific correlation parameters for each block
        index.rhoB <- paste(diff.time,
                            collapse(lp2data[lpB[,1],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                            collapse(lp2data[lpB[,2],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                            strataB, sep = sep["rho.strata"])
    }
    ## internal identifier for the parameter
    code.rhoB <- tapply(paste("B",lpB[,1],lpB[,2], sep=sep["rho.name"]), index.rhoB, FUN = identity, simplify = FALSE)
    ## strata
    strata.rhoB <- tapply(strataB, index.rhoB, unique)
    ## provide a name for the parameter
    if(structure$twin){
        level.rhoB <- paste0("(",tapply(diff.time, index.rhoB, unique),")")
    }else{
        level.rhoB <- tapply(paste0("(",diff.time,"|",collapse(lp2data[lpB[,1],name.cov,drop=FALSE], sep = sep["lp"]),",",collapse(lp2data[lpB[,2],name.cov,drop=FALSE], sep = sep["lp"]),")"),
                             index.rhoB, unique)
    }
    if(n.strata>1){
        level.rhoB <- paste(strata.rhoB, U.strata[strata.rhoB], sep = sep["rho.strata"])
    }
    
    ## ** constraint
    constraint.rhoB <- rep(NA, length(code.rhoB))

    ## ** export
    return(list(code = unname(code.rhoB),
                level = unname(level.rhoB),
                strata = unname(strata.rhoB),
                nesting = NULL,
                constraint = constraint.rhoB))
    
}

## * .crossBlockUN
.crossBlockUN <- function(structure, XpairPattern, U.strata, name.cov, sep){
    
    ## ** extract data
    lpB <- XpairPattern$lpB
    strataB <- XpairPattern$strataB
    lp2data <- structure$cor$lp2data.cross
    name.strata <-  structure$name$strata
    name.time <-  structure$name$time
    n.strata <- length(U.strata)

    ## ** identify parameters
    if(structure$twin){ ## shared correlation parameters for all blocks
        index.rhoB <- paste(lp2data[lpB[,1],name.time], lp2data[lpB[,2],name.time], strataB, sep = sep["rho.strata"])
    }else{ ## specific correlation parameters for each block
        index.rhoB <- paste(lp2data[lpB[,1],name.time], lp2data[lpB[,2],name.time],
                            collapse(lp2data[lpB[,1],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                            collapse(lp2data[lpB[,2],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                            strataB, sep = sep["rho.strata"])
    }
    ## internal identifier for the parameter
    code.rhoB <- tapply(paste("B",lpB[,1],lpB[,2], sep=sep["rho.name"]), index.rhoB, FUN = identity, simplify = FALSE)
    ## strata
    strata.rhoB <- tapply(strataB, index.rhoB, unique)
    ## provide a name for the parameter
    if(structure$twin){
        level.rhoB <- tapply(paste0("(",paste(lp2data[lpB[,1],name.time], lp2data[lpB[,2],name.time], sep = ","), ")"), index.rhoB, unique)
    }else{
        level.rhoB <- tapply(paste0("(",lp2data[lpB[,1],name.time],",", lp2data[lpB[,2],name.time],
                                    "|",collapse(lp2data[lpB[,1],name.cov,drop=FALSE], sep = sep["lp"]),",",collapse(lp2data[lpB[,2],name.cov,drop=FALSE], sep = sep["lp"]), ")"),
                             index.rhoB, unique)
    }
    if(n.strata>1){
        level.rhoB <- paste(strata.rhoB, U.strata[strata.rhoB], sep = sep["rho.strata"])
    }
browser()    
    ## ** constraint
    constraint.rhoB <- rep(NA, length(code.rhoB))

    ## ** export
    return(list(code = unname(code.rhoB),
                level = unname(level.rhoB),
                strata = unname(strata.rhoB),
                nesting = NULL,
                constraint = constraint.rhoB))
}

##----------------------------------------------------------------------
### structure-crossBlock.R ends here
