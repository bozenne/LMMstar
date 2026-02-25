### structure-skeletonRho.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 11 2023 (13:27) 
## Version: 
## Last-Updated: feb 25 2026 (18:24) 
##           By: Brice Ozenne
##     Update #: 1529
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .skeletonRho
##' @title Parametrization of Correlation Structure.
##' @description Parametrization of the correlation structure.
##' @noRd
##'
##' @param structure [structure]
##' @param data [data.frame] dataset.
##' @param U.cluster [character vector] cluster levels.
##' @param index.cluster [list of integer vector] position of each cluster in the dataset.
##' @param U.time [character vector] time levels.
##' @param index.clusterTime [list of integer vector] time index for each cluster.
##' @param U.strata [character vector] time levels.
##' @param index.clusterStrata [list of integer vector] strata index for each cluster.
##' @param sep [character vector of length 2] characters used to name the variance parameters,
##' the first is used to specify the covariate level whereas
##' the second is used to specify the strata level (when more than 1 strata).
##'
##' @keywords internal

`.skeletonRho` <-
    function(structure, data, 
             U.cluster, index.cluster,
             U.time, index.clusterTime, 
             U.strata, index.clusterStrata,
             sep) UseMethod(".skeletonRho")

## * .skeletonRho.ID
.skeletonRho.ID <- function(structure, data, 
                            U.cluster, index.cluster,
                            U.time, index.clusterTime, 
                            U.strata, index.clusterStrata,
                            sep){

    ## no correlation parameter
    return(structure)

}

## * .skeletonRho.IND
.skeletonRho.IND <- .skeletonRho.ID

## * .skeletonRho.CS
.skeletonRho.CS <- function(structure, data, 
                            U.cluster, index.cluster,
                            U.time, index.clusterTime, 
                            U.strata, index.clusterStrata,
                            sep){
    ## ** prepare
    n.strata <- length(U.strata)

    ## use to name the correlation coefficient, e.g. rho(id)
    cluster.var <- stats::na.omit(c(attr(structure$name$cluster,"original"),structure$name$cluster))[1]
    if(cluster.var=="XXcluster.indexXX"){
        cluster.varNULL <- NULL
        cluster.varPNULL <- NULL
    }else{
        cluster.varNULL <- cluster.var
        cluster.varPNULL <- paste0("(",cluster.var,")")
    }

    ## sigma parameter(s)
    param.sigma <- structure$param[structure$param$type=="sigma","name"]
    strataIndex.sigma <- structure$param[structure$param$type=="sigma","index.strata"]
    
    ## covariates that are not time nor strata
    name.cov <- na.omit(setdiff(structure$name$cor[[1]], c(structure$name$strata,structure$name$time)))
    
    ## ** handle special cases

    ## ***  no repetition
    if(all(lengths(index.cluster)==1)){return(structure)}

    ## ***  default pattern (no covariate, only strata)
    if(length(name.cov)==0){

        ## identify strata with at least on pair of observations
        vec.strataRho <- sort(unique(index.clusterStrata[lengths(index.cluster)>1]))
        if(length(U.strata)==1){
            level.strataRho <- ""
        }else{
            level.strataRho <- paste0(sep["rho.name"],U.strata[vec.strataRho])
        }
        
        ## table of correlation coefficients
        structure.rho <- data.frame(name = paste0("rho",cluster.varPNULL,level.strataRho),
                                    index.strata = vec.strataRho,
                                    type = "rho",
                                    constraint = as.numeric(NA),
                                    level = paste0(cluster.varPNULL,level.strataRho),
                                    code = paste("W",structure$cor$lp[match(vec.strataRho,index.clusterStrata)],structure$cor$lp[match(vec.strataRho,index.clusterStrata)], sep =  sep["rho.name"]),
                                    sigma = param.sigma[match(vec.strataRho, strataIndex.sigma)],
                                    k.x = NA,
                                    k.y = NA,                                  
                                    stringsAsFactors = FALSE)        
        structure$param <- rbind(structure$param, structure.rho)
        rownames(structure$param) <- NULL
        return(structure)
    }
browser()
    ## ** identify unique pairs of linear predictors
    XpairPattern <- .pairPatternX(structure$cor, 
                                  U.cluster = U.cluster, index.cluster = index.cluster,
                                  U.strata = U.strata, index.clusterStrata = index.clusterStrata)

    ## ** within-block (pairs with the same linear predictor)
    ## NOTE: XpairPattern$lpW should be constant row-wise: all(lengths(apply(XpairPattern$lpW,1,unique,simplify=FALSE))==1)    
    if(NROW(XpairPattern$lpW) == 0){ ## no within block: typically the case with crossed random effects
        code.rhoW <- NULL
        strata.rhoW <- NULL
        level.rhoW <- NULL
        constraint.rhoW <- NULL        
    }else if(structure$twin){ ## identical blocks on the diagonal or single block (no covariates)
        if(n.strata==1){
            ## internal indentifier for the parameter (single parameter)
            code.rhoW <- list(paste("W",XpairPattern$lpW[,1],XpairPattern$lpW[,2],sep=sep["rho.name"]))
            ## strata relative to the parameter
            strata.rhoW <- XpairPattern$strataW[1]
            ## provide a name for the parameter
            level.rhoW <- paste0("(",cluster.varNULL,")")
            ## constraint for the correlation parameter
            constraint.rhoW <- NA
        }else{
            ## internal indentifier for the parameter (one parameter per strata)
            code.rhoW <- tapply(1:NROW(XpairPattern$lpW), XpairPattern$strataW, FUN = function(iVec){
                paste("W",XpairPattern$lpW[iVec,1,drop=FALSE],XpairPattern$lpW[iVec,2,drop=FALSE],sep=sep[2])
            }, simplify = FALSE)
            ## strata relative to the parameter
            strata.rhoW <- tapply(XpairPattern$strataW, XpairPattern$strataW, FUN = unique)
            ## provide a name for the parameter
            level.rhoW <- paste(paste0("(",cluster.varNULL,")"),U.strata[strata.rhoW], sep = sep[3])
            ## constraint for the correlation parameter
            constraint.rhoW <- rep(NA, length(code.rhoW))
        }
    }else{ ## distinct blocks on the diagonal
        ## internal indentifier for the parameter (one parameter per strata and covariate, i.e. per linear predictor)
        code.rhoW <- as.list(paste("W",XpairPattern$lpW[,1],XpairPattern$lpW[,2],sep=sep["rho.name"]))
        ## strata relative to the parameter
        strata.rhoW <- XpairPattern$strataW
        ## provide a name for the parameter
        level.rhoW <- paste(paste0("(",nlme::collapse(structure$cor$lp2data[XpairPattern$lpW[,1],name.cov,drop=FALSE], sep = sep[1], as.factor = FALSE),")"),U.strata[strata.rhoW], sep = sep[3])
        ## constraint for the correlation parameter
        constraint.rhoW <- rep(NA, length(code.rhoW))
    }

browser()
    ## ** between-block (pairs with distinct linear predictors)
    if(NROW(XpairPattern$lpB)>0){
        outCrossBlock <- switch(structure$class["correlation.cross"],
                                ID = .crossBlockID(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, sep = sep),
                                CS = .crossBlockCS(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, sep = sep),
                                TOEPLITZ = .crossBlockTOEPLITZ(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, sep = sep),
                                DUN = .crossBlockDUN(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, sep = sep),
                                UN = .crossBlockUN(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, sep = sep)
                                )
        if(!is.null(outCrossBlock$nesting)){ ## same as structure$twin == TRUE
            ## instead of rho(id) which should be for the cross block element, uses rho(id/...) for the diagonal block
            level.rhoW <- paste0("(",paste(rev(outCrossBlock$nesting),collapse = "/"),")")
            if(n.strata>1){ 
                level.rhoW <- paste(level.rhoW, U.strata[strata.rhoW], sep = sep["rho.strata"])            
            }
        }
    }

    ## ***  retrieve k
    if("k" %in% structure$param$type){
        
        param.k <- structure$param[structure$param$type=="k","name"]
        ## k parameters
        param.k.augmented <- structure$param$name[match(rownames(structure$var$lp2X)[structure$var$lp],structure$param$code)]
        param.k.augmented[param.k.augmented %in% param.sigma] <- NA

        ## variables
        strata.var <- structure$name$strata
        reg.var <- stats::na.omit(structure$name$cor[[1]])
        n.reg <- length(reg.var)
        n.strata <- length(U.strata)
        reg.var.unequal <- intersect(reg.var,names(test.active.covariate)[test.active.covariate>0])    

        indexObs <- do.call(rbind,lapply(XpairPattern$LpU.strata,attr,"index"))
        
        k.x <- param.k.augmented[indexObs[,1]]
        k.y <- param.k.augmented[indexObs[,2]]

    }else{

        k.x <- NA
        k.y <- NA

    }

    ## ***  collect    
    code.rho <- unname(c(code.rhoW, outCrossBlock$code))
    level.rho <- unname(c(level.rhoW, outCrossBlock$level))
    strata.rho <- unname(c(strata.rhoW, outCrossBlock$strata))

    structure.rho <- data.frame(name = paste0("rho",level.rho),
                                index.strata = strata.rho,
                                type = rep("rho",length=length(level.rho)),
                                constraint = NA_real_,
                                level = level.rho,
                                code = NA,
                                sigma = param.sigma[match(strata.rho,strataIndex.sigma)],
                                k.x = k.x,
                                k.y = k.y,                                  
                                stringsAsFactors = FALSE)
    structure.rho$code <- code.rho

    ## ** export
    rownames(structure.rho) <- NULL
    structure$param <- rbind(structure$param, structure.rho)
    return(structure)
}

## * .skeletonRho.RE
.skeletonRho.RE <- function(structure, data, 
                            U.cluster, index.cluster,
                            U.time, index.clusterTime, 
                            U.strata, index.clusterStrata,
                            sep){

    ## ** handle special case (no repetition)
    if(all(lengths(index.cluster)==1)){return(structure)}

    ## ** identify variance and correlation parameters
    structure <- .skeletonRho.CS(structure = structure, data = data, 
                                 U.cluster = U.cluster, index.cluster = index.cluster,
                                 U.time = U.time, index.clusterTime = index.clusterTime, 
                                 U.strata = U.strata, index.clusterStrata = index.clusterStrata,
                                 sep = sep)

    ## ** relate correlation parameters to random effects
    n.strata <- length(U.strata)

    ## *** extract hierarchy
    hierarchy <- structure$ranef$hierarchy
    n.hierarchy <- length(hierarchy) 
    name.ranef <- unname(unlist(hierarchy))
    n.ranef <- length(name.ranef) 
    if(any(duplicated(name.ranef))){
        stop("Duplicated name(s) for the random effects: \"",paste(name.ranef[duplicated(name.ranef)], collapse = "\", \""),"\". \n",
             "Cannot associate correlation parameters and random effects. \n")
    }

    if(n.ranef==1){
        param.rho <- structure$param[structure$param$type=="rho" & is.na(structure$param$constraint),,drop=FALSE]
        structure$ranef$param <- list(matrix(param.rho$name[param.rho$index.strata], nrow = n.ranef, ncol = n.strata,
                                             dimnames = list(name.ranef,U.strata)))
    }else{

        ## *** prepare output
        structure$ranef$param <- lapply(hierarchy, function(iH){
            matrix(as.character(NA), nrow = length(iH), ncol = n.strata,
                   dimnames = list(iH,U.strata))
        })

        
        ## *** identify the variable constant outside the hierarchy and varying inside the hierarchy
        structure.activeRho <- structure$param[which(structure$param$type == "rho" & is.na(structure$param$constrain)),c("name","index.strata","code")]
        activeRho2var <- structure$rho$param2var[structure.activeRho$code,,drop=FALSE]
        code2name <- stats::setNames(structure.activeRho$name, structure.activeRho$code)
        var2col <- attr(structure$rho$param2var, "var2col")
        
        for(iStrata in 1:n.strata){ ## iStrata <- 1
            iActiveRho2var <- activeRho2var[structure.activeRho$index.strata==iStrata,,drop=FALSE]

            for(iH in 1:n.hierarchy){ ## iH <- 1                
                iHierarchy <- hierarchy[[iH]]
                iIndex.equal <- which(rowSums(iActiveRho2var[,unlist(var2col[iHierarchy]),drop=FALSE])>0)
                if(length(iIndex.equal)!=length(iHierarchy)){
                    stop("Something when wrong when associating the correlation parameters and random effects. \n",
                         "hierarchy: ",paste(iHierarchy, collapse="\", \""),"\n",
                         "correlation parameter(s): \"",paste(rownames(iActiveRho2var)[iIndex.equal], collapse="\", \""),"\"\n")
                }

                ## with multiple variable inside a hierarchy, random effects are ordered according to how many variables vary within the hierarchy
                iOrder.ranef <- sort(rowSums(iActiveRho2var[iIndex.equal,,drop=FALSE]))
                structure$ranef$param[[iH]][,iStrata] <-  code2name[names(iOrder.ranef)]
            }
        }
    }

    ## ** export
    return(structure)
}

## * .skeletonRho.TOEPLITZ
.skeletonRho.TOEPLITZ <- function(structure, data, 
                                  U.cluster, index.cluster,
                                  U.time, index.clusterTime, 
                                  U.strata, index.clusterStrata,
                                  sep){

    ## ** handle special case (no repetition)
    if(all(lengths(index.cluster)==1)){return(structure)}

    ## ** extract information
    ## parameters
    param.sigma <- structure$param[structure$param$type=="sigma","name"]
    strataIndex.sigma <- structure$param[structure$param$type=="sigma","index.strata"]

    param.k <- structure$param[structure$param$type=="k","name"]

    ## design matrix (reduce sample size to unique replicates)
    XpairPattern <- .pairPatternX(structure$cor,
                                  U.cluster = U.cluster, index.cluster = index.cluster,
                                  U.strata = U.strata, index.clusterStrata = index.clusterStrata)

    ## k parameters
    param.k.augmented <- structure$param$name[match(rownames(structure$var$lp2X)[structure$var$lp],structure$param$code)]
    param.k.augmented[param.k.augmented %in% param.sigma] <- NA

    ## variables
    strata.var <- structure$name$strata
    reg.var <- setdiff(structure$name$cor[[1]], NA)
    n.reg <- length(reg.var)
    n.strata <- length(U.strata)
    
    ## levels
    M.level.cor <- structure$cor$lp2data
    lp2X.cor <- structure$cor$lp2X
    lp2data.cor <- structure$cor$lp2data
    
    ## type
    type <- structure$type
    
    ## ** identify and name parameters

    ## *** combine linear predictor accross strata
    diffLp <- do.call(rbind,XpairPattern$diffU.strata) ## difference in linear predictor index for pair of observation
    indexLp <- do.call(rbind,XpairPattern$LpU.strata) ## linear predictor index for each observation of each pair
    lp.x <- lp2X.cor[indexLp[,1],,drop=FALSE] ## design matrix for one observation of each pair of observations
    data.x <- lp2data.cor[indexLp[,1],,drop=FALSE] ## data for one observation of each pair of observations
    data.y <- lp2data.cor[indexLp[,2],,drop=FALSE] ## data for the other observation of each pair of observations
    if(n.strata>1){
        strataLp <- data.x[[strata.var]]
    }else{
        strataLp <- rep(1, NROW(data.x))
    }
    n.Lp <- NROW(lp.x)

    ## *** identify pairs with equal or non-equal linear predictors
    index.equal <- which(rowSums(diffLp!=0)==0)
    index.unequal <- setdiff(1:n.Lp, index.equal)
    index.0 <- NULL ## cells in the design matrix constrained to be 0

    ## *** generate code
    code.rho <- rep(NA, length = n.Lp)
    level.rho <- rep(NA, length = n.Lp)
    ## pairs containing the same linear predictor
    if(length(index.equal)>0){
        stop("Something went wrong when identifying the linear predictors for the TOEPLITZ structure. \n",
             "Any two pairs of observations should have distinct linear predictors")        
    }

    ## pairs with distinct linear predictors
    if(length(index.unequal)>0){
        time.var <- utils::tail(reg.var,1)
        block.var <- setdiff(reg.var, time.var)

        if(n.strata==1){
            coltime.var <- time.var
            colblock.var <- block.var
        }else{
            ## in case of strata, possible duplicated time column (one for each strata)
            factor.corterm <- attr(attr(structure$cor$X,"formula"),"factor")
            ## identify column time
            index.coltime <- attr(structure$cor$X,"term.labels") %in% colnames(factor.corterm)[factor.corterm[time.var,]>0]
            ## redefine time.var
            coltime.var <- colnames(structure$cor$X)[index.coltime]
            colblock.var <- setdiff(colnames(structure$cor$X), coltime.var)
        }

        if(structure$bloc){ ## block
            test.diffBlock <- rowSums(diffLp[,colblock.var,drop=FALSE]!=0)

            ## same block
            index.sameBlock <- index.unequal[which(test.diffBlock==0)]
            block.sameBlock <- data.x[index.sameBlock,block.var]
            dt.sameBlock <- rowSums(diffLp[index.sameBlock,coltime.var,drop=FALSE]) ## 0 value in the other strata so summing is like taking the strata-specific value
            if(type=="UN"){
                code.rho[index.sameBlock] <- paste("R",indexLp[index.sameBlock,1],indexLp[index.sameBlock,2],sep=sep[2])
                level.rho[index.sameBlock] <- paste(block.sameBlock,paste("(",data.x[index.sameBlock,time.var],
                                                                          ",",data.y[index.sameBlock,time.var],")",sep=""),
                                                    sep=sep[2])                
            }else if(type=="LAG"){
                code.rho[index.sameBlock] <- paste("R",strataLp[index.sameBlock],block.sameBlock,dt.sameBlock,sep=sep[2])
                level.rho[index.sameBlock] <- paste(block.sameBlock,paste("(dt=",dt.sameBlock,")",sep=""),sep=sep[2])
            }else if(type=="CS"){
                code.rho[index.sameBlock] <- paste("R",strataLp[index.sameBlock],block.sameBlock,sep=sep[2])
                level.rho[index.sameBlock] <- paste0(block.sameBlock) ## convert to character
            }

            ## different blocks (with ordered times and block, e.g. block (A,B) and (B,A) are the same)
            index.diffBlock <- index.unequal[which(test.diffBlock!=0)]
            rev.diffBlock <- data.x[index.diffBlock,block.var] > data.y[index.diffBlock,block.var]
            block.diffBlock <- data.frame(x = data.x[index.diffBlock,block.var], y = data.y[index.diffBlock,block.var])
            block.diffBlock[rev.diffBlock,] <- block.diffBlock[rev.diffBlock,c("y","x")]
            t.diffBlock <- data.frame(x = data.x[index.diffBlock,time.var],
                                      y = data.y[index.diffBlock,time.var])
            t.diffBlock[rev.diffBlock,] <- t.diffBlock[rev.diffBlock,c("y","x")]
            dt.diffBlock <- rowSums(diffLp[index.diffBlock,coltime.var,drop=FALSE])
            if(type=="UN"){
                ## constrain constant correlation when measured at the same time \rho = cor(X(t),Y(t))
                indexLp.diffBlock <- indexLp[index.diffBlock,,drop=FALSE]
                indexLp.diffBlock[data.x[index.diffBlock,time.var] == data.y[index.diffBlock,time.var],] <- 0
                txt.diffBlock <- ifelse(indexLp.diffBlock[,1]==0,"dt=0",paste0(t.diffBlock$x,",",t.diffBlock$y))
            
                code.rho[index.diffBlock] <- paste("D",indexLp.diffBlock[,1],indexLp.diffBlock[,2],sep=sep[2])
                level.rho[index.diffBlock] <- paste("(",block.diffBlock$x,",",block.diffBlock$y,",",txt.diffBlock,")",sep="")
            }else if(type=="LAG"){
                code.rho[index.diffBlock] <- paste("D",strataLp[index.diffBlock],block.diffBlock$x,block.diffBlock$y,abs(dt.diffBlock),sep=sep[2])
                level.rho[index.diffBlock] <- paste("(",block.diffBlock$x,",",block.diffBlock$y,",dt=",abs(dt.diffBlock),")",sep="")
            }else if(type=="CS"){
                code.rho[index.diffBlock] <- paste("D",strataLp[index.diffBlock],block.diffBlock$x,block.diffBlock$y,as.numeric(dt.diffBlock!=0),sep=sep[2])
                level.rho[index.diffBlock] <- paste("(",block.diffBlock$x,",",block.diffBlock$y,",dt=",as.numeric(dt.diffBlock!=0),")",sep="")
            }

        }else{ ## no block: base correlation coefficient on the time difference
            ## 0 value in the other strata so summing is like taking the strata-specific value
            code.rho[index.unequal] <- paste("D",strataLp[index.unequal],rowSums(diffLp[,coltime.var,drop=FALSE]),sep=sep[2])
            level.rho[index.unequal] <- paste("(",rowSums(diffLp[,coltime.var,drop=FALSE]),")",sep="")
        }        
    }
    test.rho <- !duplicated(code.rho)
    code.Urho <- code.rho[test.rho]

    ## ***  name parameters
    level.Urho <-  level.rho[test.rho]
    strata.Urho <- strataLp[test.rho]
    if(n.strata>1){
        level.Urho <- paste0(level.Urho,sep[2],strata.Urho)
    }

    ## ***  retrieve k
    if(length(param.k)>0){
        indexObs <- do.call(rbind,lapply(XpairPattern$LpU.strata,attr,"index"))
        
        k.x <- param.k.augmented[indexObs[,1]]
        k.y <- param.k.augmented[indexObs[,2]]
    }else{
        k.x <- NA
        k.y <- NA
    }

    ## ***  collect    
    strataIndex.rho <- match(strata.Urho,U.strata)
    structure.rho <- data.frame(name = paste0("rho",level.Urho),
                                index.strata = strataIndex.rho,
                                type = rep("rho",length=length(level.Urho)),
                                constraint = as.numeric(NA),
                                level = level.Urho,
                                code = code.Urho,
                                lp.x = NA,
                                lp.y = NA,
                                sigma = param.sigma[match(strataIndex.rho,strataIndex.sigma)],
                                k.x = k.x[test.rho],
                                k.y = k.y[test.rho],                                  
                                stringsAsFactors = FALSE)

    ## save information about correlation parameters fixed to 0 (e.g. crossed random effects)
    if(length(index.0)>0){
        Ucode.rho0 <- unique(code.rho[index.0])
        if(length(intersect(Ucode.rho0,code.rho[-index.0]))){
            stop("Something went wrong with the constraint when identifying the correlation parameters. \n")
        }
        structure.rho$constraint[structure.rho$code %in% Ucode.rho0] <- 0
    }
    ## ensure proper ordering
    code.x <- tapply(indexLp[,1], code.rho, base::identity, simplify = FALSE)
    code.y <- tapply(indexLp[,2], code.rho, base::identity, simplify = FALSE)
    structure.rho$lp.x[match(names(code.x),structure.rho$code)] <- code.x
    structure.rho$lp.y[match(names(code.y),structure.rho$code)] <- code.y

    ## ** export
    rownames(structure.rho) <- NULL
    structure$param <- rbind(structure$param, structure.rho)
    return(structure)
}

## * .skeletonRho.UN
.skeletonRho.UN <- function(structure, data, 
                            U.cluster, index.cluster,
                            U.time, index.clusterTime, 
                            U.strata, index.clusterStrata,
                            sep){

    ## ** handle special case (no repetition)
    if(all(lengths(index.cluster)==1)){return(structure)}

    ## ** extract information
    ## parameters
    param.sigma <- structure$param[structure$param$type=="sigma","name"]
    strataIndex.sigma <- structure$param[structure$param$type=="sigma","index.strata"]

    param.k <- structure$param[structure$param$type=="k","name"]

    ## design matrix (reduce sample size to unique replicates)
    XpairPattern <- .pairPatternX(structure$cor,
                                  U.cluster = U.cluster, index.cluster = index.cluster,
                                  U.strata = U.strata, index.clusterStrata = index.clusterStrata)

    ## k parameters
    param.k.augmented <- structure$param$name[match(rownames(structure$var$lp2X)[structure$var$lp],structure$param$code)]
    param.k.augmented[param.k.augmented %in% param.sigma] <- NA

    ## variables
    time.var <- structure$name$cor[[1]]
    strata.var <- structure$name$strata
    n.strata <- length(U.strata)

    ## levels
    M.level <- structure$cor$lp2data

    ## ** identify and name parameters
    ## name parameters
    indexLp <- do.call(rbind,XpairPattern$LpU.strata)
    level.rho <- paste0("(",nlme::collapse(M.level[indexLp[,1],time.var,drop=FALSE], sep = sep[1], as.factor = FALSE),
                        ",",nlme::collapse(M.level[indexLp[,2],time.var,drop=FALSE], sep = sep[1], as.factor = FALSE),
                        ")")
    if(n.strata==1){
        strata.rho <- U.strata
    }else if(n.strata>1){
        strata.rho <- M.level[indexLp[,1],strata.var]
        level.rho <- paste0(level.rho,sep[2],strata.rho)        
    }
    ## generate code
    diffLp <- do.call(rbind,XpairPattern$diffU.strata)
    code.rho <- paste("D",strata.rho,nlme::collapse(diffLp, as.factor=FALSE, sep=sep[1]), sep = sep[2])

    ## retrive k
    if(length(param.k)>0){
        indexObs <- do.call(rbind,lapply(XpairPattern$LpU.strata,attr,"index"))
        k.x <- param.k.augmented[indexObs[,1]]
        k.y <- param.k.augmented[indexObs[,2]]
    }else{
        k.x <- NA
        k.y <- NA
    }

    ## collect
    strataIndex.rho <- match(strata.rho,U.strata)
    structure.rho <- data.frame(name = paste0("rho",level.rho),
                                index.strata = strataIndex.rho,
                                type = rep("rho",length=length(level.rho)),
                                constraint = as.numeric(NA),
                                level = level.rho,
                                code = code.rho,
                                lp.x = NA,
                                lp.y = NA,
                                sigma = param.sigma[match(strataIndex.rho,strataIndex.sigma)],
                                k.x = k.x,
                                k.y = k.y,                                  
                                stringsAsFactors = FALSE)

    structure.rho$lp.x <- as.list(indexLp[,1])
    structure.rho$lp.y <- as.list(indexLp[,2])
    rownames(structure.rho) <- NULL
    structure$param <- rbind(structure$param, structure.rho)

    ## ** export    
    return(structure)
}

## * .skeletonRho.EXP
## .skeletonRho.EXP <- function(structure, data, 
##                            U.cluster, index.cluster,
##                            U.time, index.clusterTime, 
##                            U.strata, index.clusterStrata){

##     ## *** param rho
##     regressor <- colnames(X.cor)[which(attr(X.cor, "assign") == max(attr(X.cor, "assign")))]

##     if(n.strata==1){
##         param.rho <- "lambda"
##         strata.rho <- 1
##         code.rho <- regressor
##         level.rho <- ""
##         if(structure$heterogeneous){
##             param.rho <- c(param.rho,"nugget")
##             strata.rho <- c(strata.rho,1)
##             code.rho <- c(code.rho,NA)
##             level.rho <- c(code.rho,"")
##         }
##     }else{
##         param.rho <- paste0("lambda",U.strata)
##         strata.rho <- 1:n.strata
##         code.rho <- regressor
##         level.rho <- U.strata
##         if(structure$heterogeneous){
##             param.rho <- c(param.rho,paste0("nugget",U.strata))
##             strata.rho <- c(strata.rho,1:n.strata)
##             code.rho <- c(code.rho,rep(NA, n.strata))
##             level.rho <- c(level.rho,U.strata)
##         }
##     }

## }

## * helper
## ** .pairPatternX
##' @description Generate combinations of pairwise differences from linear predictors
##' 
##' @param object [list] output of model.matrix containing at least the design matrix in element \code{"X"}.
##' @param U.cluster [character vector] cluster levels.
##' @param index.cluster [list of integer vector] position of each cluster in the dataset.
##' @param U.strata [character vector] strata levels.
##' @param indexCluster.strata [character vector] strata associated to each cluster.
##' @param sep [character] character used to collapse the covariate values into a linear predictor.
##' 
##' @noRd
##' @return When \code{blockwise} is TRUE, a list with the following elements \itemize{
##' \item lpW [numeric matrix, n.lp x 2] index of the linear predictor (diagonal blocks).
##' \item strataW [numeric vector, n.lp] index of the strata.
##' \item indexW [numeric matrix, n.lp x 2] index of the observations forming the pair in the dataset.
##' \item clusterW [numeric vector, n.lp] index of the cluster.
##' \item lpB [numeric matrix, n.lp x 2] index of the linear predictor (off-diagonal blocks).
##' \item diffXB [numeric matrix, n.lp n.X] difference in design matrix value in the  pair.
##' \item strataB [numeric vector, n.lp] index of the strata.
##' \item indexB [numeric matrix, n.lp x 2] index of the observations forming the pair in the dataset.
##' \item clusterB [numeric vector, n.lp] index of the cluster.
##' }
##' where n.lp is the number of distinct value that the linear predictor can take
##' and n.X the number of columns in the design matrix.
##' Otherwise, a list with the following elements \itemize{
##' \item lp [numeric matrix, n.lp x 2] index of the linear predictor (all blocks).
##' \item diffX [numeric matrix, n.lp] difference in design matrix value in the  pair.
##' \item strata [numeric vector, n.lp] index of the strata.
##' \item index [numeric matrix, n.lp x 2] index of the observations forming the pair in the dataset.
##' \item cluster [numeric vector, n.lp] index of the cluster.
##' }
##' @examples
##' data(gastricbypassL, package = "LMMstar")
##'
##' X.cor <- model.matrix(~1, data = gastricbypassL)
##' .pairPatternX(list(X=X.cor),
##'               U.cluster = unique(gastricbypassL$id),
##'               index.cluster = tapply(1:NROW(gastricbypassL), gastricbypassL$id, identity),
##'               U.strata = "1", index.clusterStrata = rep(1, length(unique(gastricbypassL$id))),
##'               blockwise = TRUE)
##' 
##' X.cor <- model.matrix(~visit, data = gastricbypassL)
##' .pairPatternX(list(X=X.cor),
##'               U.cluster = unique(gastricbypassL$id),
##'               index.cluster = tapply(1:NROW(gastricbypassL), gastricbypassL$id, identity),
##'               U.strata = "1", index.clusterStrata = rep(1, length(unique(gastricbypassL$id))),
##'               blockwise = TRUE)
.pairPatternX <- function(object, 
                          U.cluster, index.cluster,
                          U.strata, index.clusterStrata,
                          sep = c("")){

    ## ** normalize user input
    n.strata <- length(U.strata)
    
    ## ** extract unique pairs of linear predictor 
    lp.pair <- list("within" = NULL, "between" = NULL)
    index <- list("within" = NULL, "between" = NULL)
    cluster <- list("within" = NULL, "between" = NULL)
    strata <- list("within" = NULL, "between" = NULL)

    if("X.cross" %in% names(object)){
        blocks <- c("within","between")
    }else{
        blocks <- "within"
    }
        
    for(iBlock in blocks){

        iSuffix <- switch(iBlock,
                          "within" = "",
                          "between" = ".cross")

        for(iS in 1:n.strata){ ## iS <- 1
            
            ## index of the clusters in the strata
            iIndex.clusterStrata <- which(index.clusterStrata==iS)
            ## index of the clusters with unique linear predictor pattern
            iIndex.clusterU <- iIndex.clusterStrata[.patternINtable12(pattern = object[[paste0("pattern",iSuffix)]][iIndex.clusterStrata],
                                                                      pattern2lp = object[[paste0("pattern2lp",iSuffix)]])
                                                    ]
            ## linear predictor for each cluster
            iLs.lp.clusterU <- object[[paste0("pattern2lp",iSuffix)]][object[[paste0("pattern",iSuffix)]][iIndex.clusterU]]

            ## remove triplicates
            iLs.lpU.clusterU <- lapply(1:length(iIndex.clusterU),function(iC){ ## iC <- 1
                iIndex <- which(!triplicated(iLs.lp.clusterU[[iC]]))
                iIndex2 <- iIndex[order(iLs.lp.clusterU[[iC]][iIndex])]
                iOut <- iLs.lp.clusterU[[iC]][iIndex2]
                attr(iOut,"index") <- index.cluster[[iIndex.clusterU[[iC]]]][iIndex2]
                return(iOut)
            })

            ## form all pairs
            iLS.pairs <- lapply(iLs.lpU.clusterU, unorderedPairs, distinct = TRUE)
            
            iNpair.lpU.clusterU <- sapply(iLS.pairs,NCOL)
            iM.pairs <- base::t(do.call(cbind,iLS.pairs))
            ## combine while remove pairs having the same linear predictors
            iTest.duplicated <- !duplicated(iM.pairs)
            ## retrieve cluster corresponding to the linear predictor
            iCluster <- do.call(c,mapply(x = iIndex.clusterU, y = iNpair.lpU.clusterU, function(x,y){rep(x, times = y)}, SIMPLIFY = FALSE))
            ## retrieve index of the observation corresponding to the linear predictor
            iIndex <- do.call(cbind,lapply(1:length(iIndex.clusterU), function(iC){
                iiIndex <- unorderedPairs(attr(iLs.lpU.clusterU[[iC]],"index"), distinct = TRUE)
                iiTest <- iTest.duplicated[cumsum(c(1,iNpair.lpU.clusterU))[iC]:cumsum(iNpair.lpU.clusterU)[iC]]
                return(iiIndex[,iiTest,drop=FALSE])
            }))

            ## update
            lp.pair[[iBlock]] <- rbind(lp.pair[[iBlock]], iM.pairs[iTest.duplicated,,drop=FALSE])
            index[[iBlock]] <- rbind(index[[iBlock]], t(iIndex))
            cluster[[iBlock]] <- c(cluster[[iBlock]], iCluster[iTest.duplicated])
            strata[[iBlock]] <- c(strata[[iBlock]],rep(iS, sum(iTest.duplicated)))
        }
    }

    if(length(blocks)==1){
        lp.pair$between <- lp.pair$within
        index$between <- index$within
        cluster$between <- cluster$within
        strata$between <- strata$within
        object$X.cross <- object$X
    }

    ## ** contrast in linear predictor
    ##  cbind(object$X[index[,1],], object$X[index[,2],]) - lp.pair
    diffX.pair <- list(within = object$X[index$within[,2],,drop=FALSE] - object$X[index$within[,1],,drop=FALSE],
                       between = object$X.cross[index$between[,2],,drop=FALSE] - object$X.cross[index$between[,1],,drop=FALSE])
    rownames(diffX.pair$within) <- NULL
    rownames(diffX.pair$between) <- NULL

    ## ** export    
    index.within <- list(within = which(rowSums(diffX.pair$within!=0)==0),
                         between = which(rowSums(object$X[index$between[,2],,drop=FALSE] - object$X[index$between[,1],,drop=FALSE]!=0)==0))
    
    return(list(lpW = lp.pair$within[index.within$within,,drop=FALSE],
                strataW = strata$within[index.within$within],
                indexW = index$within[index.within$within,,drop=FALSE],
                clusterW = cluster$within[index.within$within],
                lpB = lp.pair$between[-index.within$between,,drop=FALSE],
                diffXB = diffX.pair$between[-index.within$between,,drop=FALSE],
                strataB = strata$between[-index.within$between],
                indexB = index$between[-index.within$between,,drop=FALSE],
                clusterB = cluster$between[-index.within$between]
                ))
}

## ** .patternINtable12
##' @description Find the subset with unique pairs from a table
##'
##' @return column/cluster index for which unique pairs are observed
##' @noRd
##' @examples
##'
##' M1 <- cbind(c(1,0,0,0),c(1,1,1,1),c(1,0,1,0))
##' .patternINtable12(M1)
##' 
##' M2 <- cbind(c(1,0,0,0,0,0),c(1,1,1,0,0,0),c(0,0,0,1,1,1),c(0,1,0,0,1,1))
##' .patternINtable12(M2)
##' 
##' M3 <- cbind(c(1,0,0,0,0,0),c(1,1,1,0,0,0),c(0,0,0,1,1,1),c(0,2,0,0,1,1))
##' .patternINtable12(M3)
##' 
##' M3.bis <- cbind(c(0,0,0,2,1,1),c(1,0,0,0,0,0),c(1,1,1,0,0,0),c(0,0,0,1,1,1))
##' .patternINtable12(M3.bis)
##'
##' M4 <- cbind(c(0,0,0,1,1,1),c(2,0,0,0,0,0),c(1,1,1,0,0,0),c(0,0,0,1,1,1))
##' .patternINtable12(M4)
.patternINtable12 <- function(object, pattern = NULL, pattern2lp = NULL){
    
    ## *** initialization
    if(missing(object) && !is.null(pattern) && !is.null(pattern2lp)){

        index.unique <- which(!duplicated(pattern)) ## only keep one cluster per pattern
        Upattern <- pattern[index.unique] ## find unique patterns
        Upattern2lp <- pattern2lp[Upattern] ## find linear predictor per pattern
        Ulp <- unique(unlist(Upattern2lp)) ## find unique linear predictor

        tablePatternLp <- do.call(cbind,lapply(Upattern2lp, function(iVec){table(factor(iVec, levels = Ulp))})) ## table of linear predictors per patterns
        object <- pmin(tablePatternLp,2) ## observing more than twice a linear predictor does not lead to different parameter
                                         ## but observing twice may lead to \rho(lp1,lp1) instead of \rho(lp1,lp2)
    }else{
        index.unique <- 1:NCOL(object)
    }

    ## *** special case 1: single pattern
    n.pattern <- NCOL(object)
    if(n.pattern==1){
        return(1)
    }

    ## *** special case 2: single linear predictor
    n.lp <- NROW(object)
    if(n.lp==1){
        return(which.max(object[1,]))
    }


    ## *** find unique patterns
    index.current <- 1:n.pattern
    index.keep <- NULL
    ## max.row <- apply(object,1,max)
    ## max.current <- rep(0,n.lp)

    for(iP in 1:n.pattern){ ## iP <- 1
        iNewKeep <- index.current[which.max(colSums(object[,index.current,drop=FALSE]))]
        index.keep <- c(index.keep,iNewKeep)

        ## update remaining patterns
        index.current <- setdiff(index.current,iNewKeep)
        if(length(index.current)>0){
            ## also remove patterns nested in the current one
            iDiff <- sweep(object[,index.current,drop=FALSE], FUN = "-", MARGIN = 1, STATS = object[,iNewKeep])
            index.current <- index.current[which(colSums(iDiff>0)>0)]
        }
        if(length(index.current)==0){ ## no other pattern remaining
            break
        }
    }

    ## *** export
    return(index.unique[index.keep])
}


##----------------------------------------------------------------------
### structure-skeletonRho.R ends here
