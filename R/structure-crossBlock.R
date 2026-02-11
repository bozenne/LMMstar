### structure-crossBlock.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 10 2026 (15:39) 
## Version: 
## Last-Updated: feb 10 2026 (17:32) 
##           By: Brice Ozenne
##     Update #: 16
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .crossBlockID
.crossBlockID <- function(structure, XpairPattern, n.strata, name.cov, sep){

    out <- .crossBlockCS(structure = structure, XpairPattern = XpairPattern, n.strata = n.strata, name.cov = name.cov, sep = sep,
                         constraint0to0 = TRUE)
    return(out)

}


## * .crossBlockCS
.crossBlockCS <- function(structure, XpairPattern, n.strata, name.cov, sep,
                          constraint0to0 = FALSE){
    
    ## ** extract data
    lpB <- XpairPattern$lpB
    strataB <- XpairPattern$strataB
    lp2data <- structure$cor$lp2data
    name.cor <- unlist(structure$name$cor)
    name.strata <-  structure$name$strata
    name.cluster <- attr(structure$name$cluster,"original")
    
    ## ** identify parameters
    if(structure$twin){ ## same correlation for all blocks at the same nesting level

        diff.nesting <- lp2data[lpB[,1],name.cov,drop=FALSE]==lp2data[lpB[,2],name.cov,drop=FALSE]
        diffS.nesting <- rowSums(diff.nesting)
        ## internal identifier for the parameter        
        code.rhoB <- tapply(paste("B",lpB[,1],lpB[,2], sep=sep[2]), paste(diffS.nesting, strataB, sep = sep[3]), FUN = identity, simplify = FALSE)
        ## order the covariates to find the ordering of the nesting
        diffE.nesting <- cbind(diff.nesting, matrix(TRUE, ncol = 1, nrow = NROW(diff.nesting), dimnames = list(NULL, name.cluster)))
        name.nesting <- names(sort(colSums(diffE.nesting), decreasing = FALSE)) ## outer nesting is when fewest pairs have equal levels
        ## strata
        strata.rhoB <- tapply(strataB, paste(diffS.nesting, strataB, sep = sep[3]), unique)
        ## provide a name for the parameter
        lsDiff.nesting <- apply(cbind(FALSE, diffE.nesting[,name.nesting,drop=FALSE]),1,diff, simplify = FALSE)
        if(any(unlist(lsDiff.nesting)<0)){
            ## if one specifies id/baseline/hour and that some observations from the same subject have same baseline but different hours
            ## while others have same hours and different baseline
            ## hour should be changed to be relative to baseline, e.g. (no,0h) (no,1h) (yes,0h) (yes,1h) ----> (no,no.0h) (no,no.1h) (yes,no.0h) (yes,no.1h)
            warning("Variables \"",paste(name.nesting,collapse="\", \""),"\" do not define a proper nesting. \n")
        }
        
        if(n.strata==1){
            level.rhoB <- paste0("(",sapply(lsDiff.nesting, function(iVec){paste(rev(names(iVec[iVec>0])),collapse="/")}),")")
        }else{
            level.rhoB <- paste0("(",sapply(lsDiff.nesting, function(iVec){paste(rev(names(iVec[iVec>0])),collapse="/")}),")")
        }
        
    }else{ ## correlation specific to each pair of linear predictors
        
        ## internal identifier for the parameter
        code.rhoB <- as.list(paste("B",lpB[,1],lpB[,2], sep=sep[2]))
        ## provide a name for the parameter
        level.rhoB <- paste0("(",
                             nlme::collapse(lp2data[lpB[,1],name.cov,drop=FALSE], sep = sep[1], as.factor = FALSE),
                             ",",
                             nlme::collapse(lp2data[lpB[,2],name.cov,drop=FALSE], sep = sep[1], as.factor = FALSE),
                             ")")
        name.nesting <- NULL
        ## strata
        strata.rhoB <- strataB
        
    }


    ## *** constraint
    if(constraint0to0){
        if(structure$twin){
            diffS.nesting <- rowSums(lp2data[lpB[,1],name.cov,drop=FALSE]==lp2data[lpB[,2],name.cov,drop=FALSE])
        }
        constraint.rhoB <- unname(ifelse(diffS.nesting==0,0,NA))
    }else{
        constraint.rhoB <- rep(NA, length(code.rhoB))
    }

    ## *** export
    return(list(code = code.rhoB,
                level = level.rhoB,
                strata = unname(strata.rhoB),
                nesting = name.nesting,
                constraint = constraint.rhoB))
}

## * .crossBlockTOEPLITZ
.crossBlockTOEPLITZ <- function(structure, XpairPattern, n.strata, name.cov, sep){

    ## *** export
    return(list(code = code.rhoB,
                level = level.rhoB,
                strata = unname(strata.rhoB),
                nesting = name.nesting,
                constraint = constraint.rhoB))

}

## * .crossBlockDUN
.crossBlockDUN <- function(structure, XpairPattern, n.strata, name.cov, sep){

    ## *** export
    return(list(code = code.rhoB,
                level = level.rhoB,
                strata = unname(strata.rhoB),
                nesting = name.nesting,
                constraint = constraint.rhoB))
    
}

## * .crossBlockUN
.crossBlockUN <- function(structure, XpairPattern, n.strata, name.cov, sep){

    ## ** extract data
    lpB <- XpairPattern$lpB
    strataB <- XpairPattern$strataB
    lp2data <- structure$cor$lp2data
    name.cor <- unlist(structure$name$cor)
    name.strata <-  structure$name$strata
    name.cluster <- attr(structure$name$cluster,"original")

    ## *** export
    return(list(code = code.rhoB,
                level = level.rhoB,
                strata = unname(strata.rhoB),
                nesting = name.nesting,
                constraint = constraint.rhoB))
}

##----------------------------------------------------------------------
### structure-crossBlock.R ends here
