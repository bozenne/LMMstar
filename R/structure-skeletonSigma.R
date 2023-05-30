### structure-skeletonSigma.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 11 2023 (11:02) 
## Version: 
## Last-Updated: maj 30 2023 (18:19) 
##           By: Brice Ozenne
##     Update #: 45
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .skeletonSigma
##' @title Parametrization of Baseline Variance Structure
##' @description Parametrization of the variance structure for the reference level.
##' @noRd
##'
##' @param structure [structure]
##' @param U.strata [character vector] strata levels
##' @param sep [character] character used to name the variance parameters, used in between \code{"sigma"} and the strata level.
##' Only relevant when more than a single strata.
##'
##' @keywords internal

`.skeletonSigma` <-
    function(structure, U.strata, sep) UseMethod(".skeletonSigma")

## * .skeletonSigma.ID
.skeletonSigma.ID <- function(structure, U.strata, sep = ":"){

    ## ** extract information
    strata.var <- structure$name$strata
    X.var <- structure$X$var
    if(NCOL(X.var)==0){
        n.strata <- 0
    }else{
        n.strata <- length(U.strata)
    }

    ## ** identify and name parameters
    if(n.strata==0){
        if(!is.null(structure$param)){
            warning("Erase existing parameter structure")
        }
        structure$param <- NULL
        return(structure)
    }else if(n.strata==1){
        param.sigma <- "sigma"
        index.sigma <- which(attr(X.var,"assign")==0)
        if(!identical(colnames(X.var)[index.sigma],"(Intercept)")){
            stop("Could not find the intercept in the design matrix for the variance.\n")
        }
        colnames(X.var)[index.sigma] <- param.sigma
        strata.sigma <- 1
    }else if(n.strata>1){
        param.sigma <- paste0("sigma:",U.strata)
        index.sigma <- which(attr(X.var,"assign")==1)
        if(!identical(colnames(X.var)[index.sigma],paste0(strata.var,U.strata))){
            stop("Could not find the strata-specific intercepts in the design matrix for the variance.\n")
        }
        colnames(X.var)[index.sigma] <- param.sigma
        strata.sigma <- 1:n.strata
        attr(X.var,"order")[index.sigma] <- 0
        
    }
    n.sigma <- length(param.sigma)

    ## ** find code associated to each parameter
    ## subset rows corresponding only to the 'baseline' variance
    if(NCOL(X.var)>length(param.sigma)){
        index.keep <- rowSums(abs(X.var[,-index.sigma,drop=FALSE]))==0
        X.Usigma <- unique(X.var[index.keep,,drop=FALSE])
    }else{
        X.Usigma <- unique(X.var)
    }
        
    ## generate code
    code.sigma <- stats::setNames(as.character(interaction(as.data.frame(X.Usigma), drop = TRUE, sep = sep)),
                                  apply(X.Usigma, 1, function(iRow){names(which(iRow==1))}))
    code.sigma <- as.character(code.sigma[param.sigma])

    ## ** find level
    M.level <- attr(X.var,"M.level")
    if(length(M.level)==0){ ## no strata nor covariate
        level.sigma <- ""
    }else if(n.strata==1 || NCOL(M.level)>1){ ## no strata or strata with covariate
        level.sigma <- paste0(".",as.character(interaction(M.level[index.sigma,,drop=FALSE],sep=sep, drop = TRUE)))
    }else if(n.strata>1){ ## only strata
        level.sigma <- paste0(sep,as.character(interaction(M.level[index.sigma,,drop=FALSE],sep=sep, drop = TRUE)))
    }

    ## ** update
    structure$X$var <- X.var
    param.sigma <- data.frame(name = param.sigma,
                              index.strata = strata.sigma,
                              type = rep("sigma",length=n.sigma),
                              index.level = index.sigma,
                              level = level.sigma,
                              code = code.sigma,
                              code.x = as.numeric(NA),
                              code.y = as.numeric(NA),
                              sigma = as.character(NA),
                              k.x = as.character(NA),
                              k.y = as.character(NA),                                  
                              stringsAsFactors = FALSE)
    if(!is.null(structure$param)){
        warning("Erase existing parameter structure")
    }
    structure$param <- param.sigma

    ## ** export
    return(structure)

}

## * skeletonSigma.IND
.skeletonSigma.IND <- .skeletonSigma.ID

## * skeletonSigma.CS
.skeletonSigma.CS <- .skeletonSigma.ID

## * skeletonSigma.RE
.skeletonSigma.RE <- .skeletonSigma.ID

## * skeletonSigma.TOEPLITZ
.skeletonSigma.TOEPLITZ <- .skeletonSigma.ID

## * skeletonSigma.UN
.skeletonSigma.UN <- .skeletonSigma.ID

## * skeletonSigma.EXP
.skeletonSigma.UN <- .skeletonSigma.ID

##----------------------------------------------------------------------
### structure-skeletonSigma.R ends here
