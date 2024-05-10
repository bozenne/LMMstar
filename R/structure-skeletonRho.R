### structure-skeletonRho.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 11 2023 (13:27) 
## Version: 
## Last-Updated: maj 10 2024 (11:43) 
##           By: Brice Ozenne
##     Update #: 1083
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
             U.strata, index.clusterStrata) UseMethod(".skeletonRho")

## * .skeletonRho.ID
.skeletonRho.ID <- function(structure, data, 
                            U.cluster, index.cluster,
                            U.time, index.clusterTime, 
                            U.strata, index.clusterStrata){

    ## no correlation parameter
    return(structure)

}

## * .skeletonRho.IND
.skeletonRho.IND <- .skeletonRho.ID

## * .skeletonRho.CS
.skeletonRho.CS <- function(structure, data, 
                            U.cluster, index.cluster,
                            U.time, index.clusterTime, 
                            U.strata, index.clusterStrata){

    ## ** handle special case (no repetition)
    if(all(lengths(index.cluster)==1)){return(structure)}

    ## ** extract information
    ## parameters
    index.sigma <- structure$param[structure$param$type=="sigma","index.level"]
    param.sigma <- structure$param[structure$param$type=="sigma","name"]
    strata.sigma <- structure$param[structure$param$type=="sigma","index.strata"]

    param.k <- structure$param[structure$param$type=="k","name"]
    strata.k <- structure$param[structure$param$type=="k","index.strata"]

    ## design matrix (reduce sample size to unique replicates)
    XpairPattern <- .pairPatternX(structure$cor,
                                  U.cluster = U.cluster, index.cluster = index.cluster,
                                  U.strata = U.strata, index.clusterStrata = index.clusterStrata)

    ## k parameters
    param.k.augmented <- structure$param$name[match(rownames(structure$var$lp2X)[structure$var$lp],structure$param$code)]
    param.k.augmented[param.k.augmented %in% param.sigma] <- NA

    ## variables
    cluster.var <- stats::na.omit(c(attr(structure$name$cluster,"original"),structure$name$cluster))[1]
    if(cluster.var=="XXcluster.indexXX"){
        cluster.varNULL <- NULL
        cluster.varPNULL <- NULL
    }else{
        cluster.varNULL <- cluster.var
        cluster.varPNULL <- paste0("(",cluster.var,")")
    }
    strata.var <- structure$name$strata
    reg.var <- stats::na.omit(structure$name$cor[[1]])
    n.reg <- length(reg.var)
    n.strata <- length(U.strata)

    ## levels
    M.level.cor <- structure$cor$lp2data
    lp2X.cor <- structure$cor$lp2X
    lp2data.cor <- structure$cor$lp2data
    
    ## sep
    sep <- LMMstar.options()$sep[c("rho.name","rho.strata")]
    

    ## ** special case (no covariate, only strata)
    if(length(reg.var)==0){
        vec.strataRho <- which(lengths(XpairPattern$LpU.strata)>0)
        if(n.strata==1){
            level.strataRho <- ""
        }else{
            level.strataRho <- paste0(sep[2],U.strata[vec.strataRho])
        }
        structure.rho <- data.frame(name = paste0("rho",cluster.varPNULL,level.strataRho),
                                    index.strata = vec.strataRho,
                                    type = "rho",
                                    constraint = as.numeric(NA),
                                    level = paste0(cluster.varPNULL,level.strataRho),
                                    code = paste0("R.",vec.strataRho,".",vec.strataRho),
                                    lp.x = NA,
                                    lp.y = NA,
                                    sigma = param.sigma[match(vec.strataRho, strata.sigma)],
                                    k.x = NA,
                                    k.y = NA,                                  
                                    stringsAsFactors = FALSE)        
        structure.rho$lp.x <- as.list(vec.strataRho)
        structure.rho$lp.y <- as.list(vec.strataRho)
        structure$param <- rbind(structure$param, structure.rho)
        rownames(structure$param) <- NULL
        return(structure)
    }

    ## ** identify and name parameters

    ## *** combine linear predictor accross strata
    diffLp <- do.call(rbind,XpairPattern$diffU.strata) ## difference in linear predictor index for pair of observation
    if(n.strata==1){
        diffLp.name <- colnames(diffLp)
    }else{ ## remove mention of the level relative to the strata in the column name to have nicer naming, i.e. rho(school):DK instead of rho(id:countryDK):DK
        diffLp.name <- gsub(paste(paste0(":",strata.var,U.strata,"$"),collapse="|"),"",colnames(diffLp))
    }
    indexLp <- do.call(rbind,XpairPattern$LpU.strata) ## linear predictor index for each observation of each pair
    clusterLp <- do.call(rbind,lapply(XpairPattern$LpU.strata, attr, "cluster"))
    lp.x <- lp2X.cor[indexLp[,1],,drop=FALSE] ## design matrix for one observation of each pair of observations
    data.x <- lp2data.cor[indexLp[,1],,drop=FALSE] ## data for one observation of each pair of observations
    data.y <- lp2data.cor[indexLp[,2],,drop=FALSE] ## data for the other observation of each pair of observations
    strataLp <- index.clusterStrata[clusterLp[,1]] ## strata for pair of observation (would be the same with clusterLp[,2])
    
    n.Lp <- NROW(lp.x)
    
    ## *** identify pairs with equal or non-equal linear predictors
    index.equal <- which(rowSums(diffLp!=0)==0)
    index.unequal <- setdiff(1:n.Lp, index.equal)
    index.0 <- NULL ## cells in the design matrix constrained to be 0

    ## *** generate code
    code.rho <- rep(NA, length = n.Lp)
    level.rho <- rep(NA, length = n.Lp)

    ## pairs with distinct linear predictors
    if(length(index.unequal)>0){
        M.level.cor <- attr(structure$cor$X,"M.level")
        test.active.covariate <- colSums(M.level.cor[colSums(diffLp!=0)>0,,drop=FALSE]!=FALSE)
        reg.var.unequal <- intersect(reg.var,names(test.active.covariate)[test.active.covariate>0])
        if(structure$type == "homogeneous"){
            ## contrast at the design matrix level
            test.n0 <- 1*(data.x[index.unequal,reg.var.unequal,drop=FALSE]!=data.y[index.unequal,reg.var.unequal,drop=FALSE])
            code.rho[index.unequal] <- paste("D",strataLp[index.unequal],nlme::collapse(test.n0,sep=sep[1],as.factor=FALSE),sep=sep[2])
            if(all(as.numeric(factor(code.rho[index.unequal], levels = unique(code.rho[index.unequal]))) == as.numeric(factor(strataLp[index.unequal], levels = unique(strataLp[index.unequal]))))){
                ## name according to cluster variable when a single within subject covariate per strata
                level.rho[index.unequal] <- cluster.varPNULL
            }else{
                ## name according to the constant variables when multiple within covariates per strata
                index.example <- !duplicated(code.rho[index.unequal])
                diffLp.example <- diffLp[index.unequal[index.example],colSums(diffLp!=0)>0,drop=FALSE] ## exclude columns corresponding to between subject covariates
                strataLp.example <- strataLp[index.unequal[index.example]]
                if(n.strata>1){
                    for(iS in 1:n.strata){ ## iS <- 1
                        diffLp.example[strataLp.example==iS,M.level.cor[strata.var]!=U.strata[iS]] <- NA
                    }
                }
                ls.rhovar <- apply(diffLp.example, MARGIN = 1, function(iRow){                    
                    return(paste(c(cluster.varNULL,diffLp.name[which(iRow==0)]),collapse="/"))
                }, simplify = FALSE)
                label.Urho.unequal <- stats::setNames(paste0("(",unlist(ls.rhovar),")"), code.rho[index.unequal][index.example])
                level.rho[index.unequal] <- unname(label.Urho.unequal[code.rho[index.unequal]])
            }
            
        }else if(structure$type == "heterogeneous"){
            code.rho[index.unequal] <- paste("D",rownames(lp.x)[index.unequal],nlme::collapse(diffLp[index.unequal,,drop=FALSE], sep = sep[1], as.factor = FALSE),sep=sep[2])
            level.rho[index.unequal] <- paste0("(",nlme::collapse(data.x[index.unequal,reg.var.unequal,drop=FALSE], sep = sep[1], as.factor = FALSE),
                                               ",",nlme::collapse(data.y[index.unequal,reg.var.unequal,drop=FALSE], sep = sep[1], as.factor = FALSE),
                                               ")")
        }
    }

    ## pairs containing the same linear predictor
    if(length(index.equal)>0){
        ## if(length(reg.var)==0){ ## case with no covariate (only strata) already dealt before
        if(structure$type == "heterogeneous" || length(index.unequal)==0){ ## homogeneous case with only across cluster covariates or heterogeneous case
            code.rho[index.equal] <- paste("R",rownames(lp.x)[index.equal],sep=sep[2])
            level.rho[index.equal] <- paste0("(",nlme::collapse(data.x[index.equal,,drop=FALSE], sep = sep[1], as.factor = FALSE),")")
        }else if(length(reg.var)==length(reg.var.unequal)){ ## homogeneous case where all covariates are within cluster
            code.rho[index.equal] <- paste("R",strataLp[index.equal],sep=sep[2])
            ## unique to handle the case with strata
                level.rho[index.equal] <- paste0("(",paste(c(cluster.varNULL,unique(diffLp.name)),collapse="/"),")")
        }else{ ## homogeneous case where some covariates are within cluster but some are across clusters
            stop("Cannot define the correlation parameters. \n",
                 "The right-hand side of the formula mixes within-cluster and between cluster covariates.\n ",
                 "Consider putting the covariate \"",paste0(setdiff(reg.var,reg.var.unequal),collapse = "\", \""),"\" on the left hand side of the formula. \n")
        }
    }


    ## *** unique correlation coefficients
    test.rho <- !duplicated(code.rho)
    code.Urho <- code.rho[test.rho]
    strata.Urho <- strataLp[test.rho]
    if(n.strata==1){
        level.Urho <-  level.rho[test.rho]
        nameX.intercept <- cluster.var
    }else{
        level.Urho <- paste0(level.rho[test.rho],sep[2],U.strata[strata.Urho])
        nameX.intercept <- paste0(cluster.var,sep[2],structure$name$strata,U.strata)
    }

    ## find which columns in the design matrix are identical between the pairs relative to a given correlation coefficients
    ls.cstCol <- by(diffLp[test.rho,,drop=FALSE], code.rho[test.rho], function(iM){
        apply(iM, MARGIN = 2, function(iCol){all(iCol==0)})
    }, simplify = FALSE)
    cstCol.Urho <- cbind(matrix(TRUE, nrow = length(code.Urho), ncol = n.strata,
                                dimnames = list(NULL, nameX.intercept)),
                         do.call(rbind,ls.cstCol)[code.Urho,,drop=FALSE])

    ## convertion between variables and columns names (useful in presence of strata: time --> id:strataA, id:strataB)
    attr(cstCol.Urho,"col2var") <- lapply(attr(structure$cor$X,"ls.level"), function(iVec){setdiff(names(iVec),strata.var)})
    attr(cstCol.Urho,"var2col") <- c(stats::setNames(list(nameX.intercept),cluster.var),
                                     tapply(names(unlist(attr(cstCol.Urho,"col2var"))),unlist(attr(cstCol.Urho,"col2var")), base::identity, simplify = FALSE)
                                     )

    ## constrain independence in the case of different groups of correlation parameters
    ## i.e. remove some of the correlation parameters
    if(length(unique(structure$group.type))>1){

        var2col <- attr(cstCol.Urho,"var2col")
        ## parameters with non-identical variables accross different groups of correlation parameters should be removed
        ls.index.0 <- tapply(names(structure$group.type),structure$group.type, function(iVar){ ## iVar <- names(structure$group.type)[1]
            iCol.group <- unlist(var2col[names(var2col) %in% iVar])
            iCol.othergroup <- unlist(var2col[names(var2col) %in% setdiff(names(structure$group.type),iVar)])
            iRho.group <- intersect(which(rowSums(cstCol.Urho[,iCol.group,drop=FALSE]==TRUE)>0), ## at least one non-identical variable
                                    which(rowSums(cstCol.Urho[,iCol.othergroup,drop=FALSE]==TRUE)>0)) ## some variables corresponding to other groups are not identical
            return(iRho.group)
        }, simplify = FALSE)

        
        index.0 <- sort(unique(c(unlist(ls.index.0), which(rowSums(cstCol.Urho[,-1,drop=FALSE]==TRUE)==0)))) ## add parameter corresponding to pairs that are all different
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
    structure.rho <- data.frame(name = paste0("rho",level.Urho),
                                index.strata = strata.Urho,
                                type = rep("rho",length=length(level.Urho)),
                                constraint = as.numeric(NA),
                                level = level.Urho,
                                code = code.Urho,
                                lp.x = NA,
                                lp.y = NA,
                                sigma = param.sigma[match(strata.Urho,strata.sigma)],
                                k.x = k.x[test.rho],
                                k.y = k.y[test.rho],                                  
                                stringsAsFactors = FALSE)

    ## save information about correlation parameters fixed to 0 (e.g. crossed random effects)
    if(length(index.0)>0){
        structure.rho$constraint[structure.rho$code %in% code.Urho[index.0]] <- 0
    }

    ## ensure proper ordering
    code.x <- tapply(indexLp[,1], code.rho, base::identity, simplify = FALSE)
    code.y <- tapply(indexLp[,2], code.rho, base::identity, simplify = FALSE)
    structure.rho$lp.x[match(names(code.x),structure.rho$code)] <- code.x
    structure.rho$lp.y[match(names(code.y),structure.rho$code)] <- code.y

    ## ** export
    rownames(structure.rho) <- NULL
    structure$param <- rbind(structure$param, structure.rho)
    structure$rho$param2var <- cstCol.Urho
    return(structure)
}

## * .skeletonRho.RE
.skeletonRho.RE <- function(structure, data, 
                            U.cluster, index.cluster,
                            U.time, index.clusterTime, 
                            U.strata, index.clusterStrata){

    ## ** handle special case (no repetition)
    if(all(lengths(index.cluster)==1)){return(structure)}

    ## ** identify variance and correlation parameters
    structure <- .skeletonRho.CS(structure = structure, data = data, 
                                 U.cluster = U.cluster, index.cluster = index.cluster,
                                 U.time = U.time, index.clusterTime = index.clusterTime, 
                                 U.strata = U.strata, index.clusterStrata = index.clusterStrata)

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

        ## *** check crossed random effects structure
        ## it is expected to have perfectly crossed factors, e.g. (a,1) (a,2) (a,3) ... (a,n) (b,1) (b,2) ... (b,n)
        ##                                                 , e.g. (A,a,1) (A,b,2) (A,c,3) ... (B,a,2)
        ## but not (a,1) (a,1) or (A,a,1) (A,a,2)
        var2col <- attr(structure$rho$param2var, "var2col")
        if(length(unique(structure$ranef$index.terms))>1){
            param.constraint <- structure$param[!is.na(structure$param$constraint) & structure$param$type == "rho",c("name","index.strata","level","code"),]
            param2var.contraint <- structure$rho$param2var[param.constraint$code,,drop=FALSE]

            if(any(param2var.contraint[,-1])){ ## first column corresponds to cluster
                warning("Crossed random effects have been specified but the associated factors are not perfectly crossed. \n",
                        "Corresponding correlation parameter have been constrain to 0 but may lead to a rank-deficient variance-covariance matrix. \n")
            }
        }
        

        ## *** identify the variable constant outside the hierarchy and varying inside the hierarchy
        structure.activeRho <- structure$param[which(structure$param$type == "rho" & is.na(structure$param$constrain)),c("name","index.strata","code")]
        activeRho2var <- structure$rho$param2var[structure.activeRho$code,,drop=FALSE]
        code2name <- stats::setNames(structure.activeRho$name, structure.activeRho$code)

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
                                  U.strata, index.clusterStrata){

    ## ** handle special case (no repetition)
    if(all(lengths(index.cluster)==1)){return(structure)}

    ## ** extract information
    ## parameters
    index.sigma <- structure$param[structure$param$type=="sigma","index.level"]
    param.sigma <- structure$param[structure$param$type=="sigma","name"]
    strata.sigma <- structure$param[structure$param$type=="sigma","index.strata"]

    param.k <- structure$param[structure$param$type=="k","name"]
    strata.k <- structure$param[structure$param$type=="k","index.strata"]

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
    
    ## sep
    sep <- LMMstar.options()$sep[c("rho.name","rho.strata")]

    ## type
    type <- structure$type
    
    ## ** identify and name parameters

    ## *** combine linear predictor accross strata
    diffLp <- do.call(rbind,XpairPattern$diffU.strata) ## difference in linear predictor index for pair of observation
    indexLp <- do.call(rbind,XpairPattern$LpU.strata) ## linear predictor index for each observation of each pair
    lp.x <- lp2X.cor[indexLp[,1],,drop=FALSE] ## design matrox for one observation of each pair of observations
    data.x <- lp2data.cor[indexLp[,1],,drop=FALSE] ## data for one observation of each pair of observations
    data.y <- lp2data.cor[indexLp[,2],,drop=FALSE] ## data for the other observation of each pair of observations
    strataLp <- index.clusterStrata[attr(index.cluster,"vectorwise")[indexLp[,1]]] ## strata for pair of observation (would be the same with indexLp[,2])
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
        if(structure$bloc){ ## block
            time.var <- reg.var[2]
            block.var <- reg.var[1]
            test.diffBlock <- rowSums(diffLp[,block.var,drop=FALSE]!=0)

            ## same block
            index.sameBlock <- index.unequal[which(test.diffBlock==0)]
            block.sameBlock <- data.x[index.sameBlock,block.var]
            dt.sameBlock <- diffLp[index.sameBlock,time.var]
            if(type=="UN"){
                code.rho[index.sameBlock] <- paste("R",indexLp[index.sameBlock,1],indexLp[index.sameBlock,2],sep=sep[2])
                level.rho[index.sameBlock] <- paste(block.sameBlock,paste("(",data.x[index.sameBlock,time.var],",",data.y[index.sameBlock,time.var],")",sep=""),sep=sep[2])                
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
            t.diffBlock <- data.frame(x = data.x[index.diffBlock,time.var], y = data.y[index.diffBlock,time.var])
            t.diffBlock[rev.diffBlock,] <- t.diffBlock[rev.diffBlock,c("y","x")]
            dt.diffBlock <- diffLp[index.diffBlock,time.var]
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
            time.var <- reg.var[1]
            code.rho[index.unequal] <- paste("D",strataLp[index.unequal],diffLp[,time.var],sep=sep[2])
            level.rho[index.unequal] <- paste("(",diffLp[,time.var],")",sep="")
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
    structure.rho <- data.frame(name = paste0("rho",level.Urho),
                                index.strata = strata.Urho,
                                type = rep("rho",length=length(level.Urho)),
                                constraint = as.numeric(NA),
                                level = level.Urho,
                                code = code.Urho,
                                lp.x = NA,
                                lp.y = NA,
                                sigma = param.sigma[match(strata.Urho,strata.sigma)],
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
                            U.strata, index.clusterStrata){

    ## ** handle special case (no repetition)
    if(all(lengths(index.cluster)==1)){return(structure)}

    ## ** extract information
    ## parameters
    index.sigma <- structure$param[structure$param$type=="sigma","index.level"]
    param.sigma <- structure$param[structure$param$type=="sigma","name"]
    strata.sigma <- structure$param[structure$param$type=="sigma","index.strata"]

    param.k <- structure$param[structure$param$type=="k","index.level"]
    param.k <- structure$param[structure$param$type=="k","name"]
    strata.k <- structure$param[structure$param$type=="k","index.strata"]

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

    ## sep
    sep <- LMMstar.options()$sep[c("rho.name","rho.strata")]

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
    structure.rho <- data.frame(name = paste0("rho",level.rho),
                                index.strata = strata.rho,
                                type = rep("rho",length=length(level.rho)),
                                constraint = as.numeric(NA),
                                level = level.rho,
                                code = code.rho,
                                lp.x = NA,
                                lp.y = NA,
                                sigma = param.sigma[match(strata.rho,strata.sigma)],
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

## ## * .skeletonRho.EXP
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
##' @description Generate unique combinations of pairwise differences from linear predictors
##' 
##' @noRd
##' @return A list with the following elements \itemize{
##' \item LpU.strata [list of matrices] linear predictor (and position of the observations in the attribute index) used for the pairwise differences.
##' \item diffU.strata [list of matrices] pairwise difference between the linear predictor index
##' }
##' @examples
##' data(gastricbypassL, package = "LMMstar")
##'
##' X.cor <- model.matrix(~1, data = gastricbypassL)
##' .pairPatternX(X.cor,
##'               U.cluster = unique(gastricbypassL$id),
##'               index.cluster = tapply(1:NROW(gastricbypassL), gastricbypassL$id, identity),
##'               U.strata = "1", index.clusterStrata = rep(1, length(unique(gastricbypassL$id))))
##' 
##' X.cor <- model.matrix(~visit, data = gastricbypassL)
##' .pairPatternX(list(X=X.cor),
##'               U.cluster = unique(gastricbypassL$id),
##'               index.cluster = tapply(1:NROW(gastricbypassL), gastricbypassL$id, identity),
##'               U.strata = "1", index.clusterStrata = rep(1, length(unique(gastricbypassL$id))))
.pairPatternX <- function(object, 
                          U.cluster, index.cluster,
                          U.strata, index.clusterStrata,
                          sep = c("")){

    ## ** normalize user input
    object.X <- object$X
    object.lp <- object$lp
    object.pattern <- object$pattern
    object.pattern2lp <- object$pattern2lp
    n.strata <- length(U.strata)

    ## code for linear predictor for each observation
    if(is.null(object.lp)){
        lp.obs <- nlme::collapse(object.X, sep = "", as.factor = TRUE)
    }else{
        lp.obs <- object.lp
    }

    ## same but grouped by cluster
    if(is.null(object.pattern)){
        pattern <- lapply(U.cluster, function(iC){ lp.obs[index.cluster[[iC]]] })
    }else{
        pattern <- object.pattern
    }

    ## ** extract pattern
    ## list containing for each strata the unique pairs of linear predictors
    ## (and the position in the dataset of these linear predictors as an attribute)
    lpU.clusterStrata <- lapply(1:n.strata, function(iS){ ## iS <- 1

        ## index of the clusters in the strata
        iIndex.clusterStrata <- which(index.clusterStrata==iS)
        ## index of the clusters with unique linear predictor pattern
        iIndex.clusterU <- iIndex.clusterStrata[.patternINtable12(pattern = pattern[iIndex.clusterStrata], pattern2lp = object.pattern2lp)]
        ## linear predictor for each cluster
        iLs.lp.clusterU <- object.pattern2lp[pattern[iIndex.clusterU]]

        ## remove triplicates
        iLs.lpU.clusterU <- lapply(1:length(iIndex.clusterU),function(iC){ ## iC <- 1
            iIndex <- which(!triplicated(iLs.lp.clusterU[[iC]]))
            iIndex2 <- iIndex[order(iLs.lp.clusterU[[iC]][iIndex])]
            iOut <- iLs.lp.clusterU[[iC]][iIndex2]
            attr(iOut,"index") <- index.cluster[[iIndex.clusterU[[iC]]]][iIndex2]
            attr(iOut,"cluster") <- rep(iIndex.clusterU[[iC]],length(iIndex2))
            return(iOut)
        })
        iVec.lpU.clusterU <- do.call(base::c,iLs.lpU.clusterU)

        ## form all pairs
        iM <- base::t(do.call(cbind,lapply(iLs.lpU.clusterU, unorderedPairs, distinct = TRUE)))
        iMatch.lpU.clusterU <- match(iM[],iVec.lpU.clusterU)

        attr(iM,"index") <- matrix(do.call(c,lapply(iLs.lpU.clusterU,attr,"index"))[iMatch.lpU.clusterU],
                                   nrow = NROW(iM), ncol = NCOL(iM))
        attr(iM,"cluster") <- matrix(do.call(c,lapply(iLs.lpU.clusterU,attr,"cluster"))[iMatch.lpU.clusterU],
                                     nrow = NROW(iM), ncol = NCOL(iM))

        ## remove duplicates
        iIndex.Upairs <- which(!duplicated(iM))
        iOut <- iM[iIndex.Upairs,,drop=FALSE]
        attr(iOut,"index") <- attr(iM,"index")[iIndex.Upairs,,drop=FALSE]
        attr(iOut,"cluster") <- attr(iM,"cluster")[iIndex.Upairs,,drop=FALSE]
        return(iOut)
    })

    ## list containing for each strata the unique difference in design matrix between pairs
    diffU.strata <- lapply(lpU.clusterStrata, function(iS){
        iM <- object.X[attr(iS,"index")[,2],,drop=FALSE] - object.X[attr(iS,"index")[,1],,drop=FALSE]
        rownames(iM) <- NULL
        return(iM)
    })

    ## ** export    
    out <- list(LpU.strata = stats::setNames(lpU.clusterStrata, U.strata),
                diffU.strata = stats::setNames(diffU.strata, U.strata))
    return(out)
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
