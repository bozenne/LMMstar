### structure-crossBlock.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 10 2026 (15:39) 
## Version: 
## Last-Updated: mar  6 2026 (14:28) 
##           By: Brice Ozenne
##     Update #: 85
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

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
