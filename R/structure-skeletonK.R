### structure-skeletonK.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 11 2023 (11:55) 
## Version: 
## Last-Updated: jul 20 2023 (17:16) 
##           By: Brice Ozenne
##     Update #: 73
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .skeletonK
##' @title Parametrization of Variance Multiplier Structure.
##' @description Parametrization of the ratio between the variance for the reference level and the other levels.
##' @noRd
##'
##' @param structure [structure]
##' @param U.strata [character vector] strata levels
##' @param sep [character vector of length 2] characters used to name the variance parameters,
##' the first is used to specify the covariate level whereas
##' the second is used to specify the strata level (when more than 1 strata).
##'
##' @keywords internal

`.skeletonK` <-
    function(structure, U.strata) UseMethod(".skeletonK")

## * .skeletonK.ID
.skeletonK.ID <- function(structure, U.strata){

    ## no variance multipler
    return(structure)

}

## * .skeletonK.IND
.skeletonK.IND <- function(structure, U.strata){

    ## ** extract information
    strata.var <- structure$name$strata
    X.var <- structure$var$X
    lp2X.var <- structure$var$lp2X
    if(NCOL(X.var)==0){
        n.strata <- 0
    }else{
        n.strata <- length(U.strata)
    }
    param.sigma <- structure$param[structure$param$type=="sigma","name"]
    strata.sigma <- structure$param[structure$param$type=="sigma","index.strata"]
    
    sep <- LMMstar.options()$sep[c("k.cov","k.strata")]

    ## ** identify and name parameters
    index.k <- which(attr(X.var,"assign")>(n.strata>1)) ## complement with structure$param$index.level
    if(length(index.k)==0){return(structure)} ## no variance multiplier

    M.level.k <- attr(X.var,"M.level")[index.k,,drop=FALSE]
    vec.level.k <- nlme::collapse(lapply(1:NCOL(M.level.k), function(iVar){
        paste0(names(M.level.k)[iVar],M.level.k[,iVar])
    }), sep = sep[2], as.factor = FALSE) 

    if(n.strata==1){ ## no strata
        level.k <- paste0(sep[1],nlme::collapse(M.level.k, sep = sep[2], as.factor = FALSE))
    }else if(any(colnames(M.level.k) %in% strata.var == FALSE)){ ## strata with covariate(s) - make sure the strata variable is at the end
        level.k <- paste0(sep[1],nlme::collapse(M.level.k[,c(setdiff(colnames(M.level.k),strata.var),strata.var),drop=FALSE], sep = sep[2], as.factor = FALSE))
    }else{ ## only strata
        level.k <- paste0(sep[2],nlme::collapse(M.level.k, sep = sep[2], as.factor = FALSE))
    }
    if(!identical(colnames(X.var)[index.k],vec.level.k)){
        stop("Could not find the k parameters in the design matrix for the variance.\n")
    }    
    colnames(X.var)[index.k] <- paste0("k",level.k)
    colnames(lp2X.var)[index.k] <- paste0("k",level.k)
    param.k <- colnames(X.var)[index.k]
    n.k <- length(param.k)

    ## ** find code associated to each parameter
    ## subset rows corresponding only to the multipliers i.e. non-0 intercept and covariate value in the design matrix
    X.Uk <- lp2X.var[rowSums(abs(lp2X.var))>1,,drop=FALSE]

    ## generate code
    code.k <- stats::setNames(rownames(X.Uk),
                              apply(X.Uk[,param.k,drop=FALSE], 1, function(iRow){names(which(iRow==1))}))[param.k]
    code.k <- as.character(code.k[param.k])

    if(n.strata == 1){
        strata.k <- stats::setNames(rep(1,n.k),param.k)
    }else{
        strata.k <- stats::setNames(match(M.level.k[,strata.var],U.strata),param.k)
    }

    ## ** update
    structure$var$X <- X.var
    structure$var$lp2X <- lp2X.var
    structure.k <- data.frame(name = param.k,
                              index.strata = strata.k,
                              type = rep("k",length=n.k),
                              constraint = as.numeric(NA),
                              level = level.k,
                              code = code.k,
                              lp.x = as.numeric(NA),
                              lp.y = as.numeric(NA),
                              sigma = param.sigma[match(strata.k,strata.sigma)],
                              k.x = as.character(NA),
                              k.y = as.character(NA),                                  
                              stringsAsFactors = FALSE)
    rownames(structure.k) <- NULL
    structure$param <- rbind(structure$param, structure.k)
    ## attr(structure$param,"level.var") <- c(attr(structure$param,"level.var"), outK$code)

    ## ** export
    return(structure)

}


## * skeletonK.CS
.skeletonK.CS <- .skeletonK.IND

## * skeletonK.RE
.skeletonK.RE <- .skeletonK.IND

## * skeletonK.TOEPLITZ
.skeletonK.TOEPLITZ <- .skeletonK.IND

## * skeletonK.UN
.skeletonK.UN <- .skeletonK.IND

## * skeletonK.EXP
.skeletonK.UN <- .skeletonK.IND

##----------------------------------------------------------------------
### structure-skeletonK.R ends here
