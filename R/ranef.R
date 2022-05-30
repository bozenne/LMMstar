### ranef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 26 2022 (11:18) 
## Version: 
## Last-Updated: May 30 2022 (22:45) 
##           By: Brice Ozenne
##     Update #: 165
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .ranef
##' @description estimate random effect in a given strata
##' @param object a \code{lmm} object.
##' @param p [numeric vector] value of the model coefficients to be used. Only relevant if differs from the fitted values.
##' @param nestingStructure [list] output of the \code{.nestingRanef} function.
##' @noRd
##'
##' @details Consider the following mixed model:
##' \deqn{Y = X\beta + \epsilon}
##' where \eqn{\Sigma_{\epsilon}}, the variance of \eqn{\epsilon}, has a (possibly stratified) compound symmetry structure.
##' Denoting by \eqn{I} the identity matirx, this mean that \eqn{\Sigma_{\epsilon} = \sigma^2 I + Z \Sigma_{\eta} Z^T}
##' where \(\Sigma_{\eta}\) is the covariance relative to the design matrix \eqn{Z} (e.g. same student or school). So implicitely we have:
##' \deqn{Y = X\beta + Z \eta + \varepsilon}
##' where \eqn{\varepsilon \sim \mathcal{N}(0, \sigma^2 I)}. So we can estimate the random effets via:
##' \deqn{E[Y|\eta] = Z \Sigma_{\eta} \Omega^{-1} (Y-X\beta)}
.ranef <- function(object, p = NULL, nestingStructure = NULL){


    ## ** extract from object
    param.name <- object$design$param$name
    param.type <- stats::setNames(object$design$param$type,param.name)
    param.rho <- param.name[param.type=="rho"]

    cluster.var <- object$cluster$var
    U.cluster <- object$design$cluster$levels
    attr(cluster.var,"original") <- NULL
    X.cor <- object$design$vcov$X$cor

    Xpattern.cor <- object$design$vcov$X$Xpattern.cor
    index.cluster <- object$design$index.cluster
    Vindex.cluster <- attr(index.cluster, "vectorwise")
    pattern.cluster <- object$design$vcov$X$pattern.cluster$pattern
    Upattern <- object$design$vcov$X$Upattern
    index.na <- object$index.na
    
    ## ** normalize user input
    if(!is.null(p)){
        if(any(duplicated(names(p)))){
            stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
        }
        if(any(param.name %in% names(p) == FALSE)){
            stop("Incorrect argument \'p\': missing parameter(s) \"",paste(param.name[param.name %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
        }
        p <- p[param.name]
    }else{
        p <- object$param
    }
    cumtau <- coef(object, p = p, transform.rho = "cov", transform.names = FALSE)
    
    ## ** converting correlation parameters into random effect variance
    if(missing(nestingStructure)){
        nestingStructure <- .nestingRanef(object)
    }
    strata.var <- attr(nestingStructure,"strata.var")
    nesting.var <- attr(nestingStructure,"nesting.var")
    index.clusterStrata <- as.character(attr(nestingStructure,"index.clusterStrata"))

    ls.tau <- lapply(nestingStructure, function(iVec){ ## iVec <- nestingStructure[[1]]
        iTau <- rev(c(utils::tail(cumtau[iVec],1),diff(rev(cumtau[iVec]))))
        if(any(iTau<0)){
            stop("Variance for the random effects is found to be negative - cannot estimate the random effects. \n")
        }
        return(iTau)
    })
    tau <- unlist(unname(ls.tau))

    ## ** extract raw residuals
    df.epsilon <- stats::residuals(object, p = p, keep.data = TRUE, type = "response", format = "long")
    if(length(index.na)>0){
        df.epsilon <- df.epsilon[-index.na,,drop=FALSE]
    }
    ## ** extract inverse residual variance-covariance matrix
    U.indexcluster <- sort(unique(df.epsilon$XXcluster.indexXX))
    ls.OmegaM1 <- stats::sigma(object, p = p, cluster = U.indexcluster, inverse = TRUE, simplifies = FALSE)
    ls.epsilon <- base::tapply(df.epsilon$r.response,df.epsilon$XXcluster.indexXX,list)## split residuals by id

    ## ** design matrix
    OmegaM1epsilon <- mapply(x = ls.OmegaM1, y = ls.epsilon, FUN = `%*%`, SIMPLIFY = FALSE)

    ## ** estimate first random effect
    ls.tau1 <- lapply(ls.tau,function(iTau){utils::tail(iTau,1)})
    G <- unlist(ls.tau1[index.clusterStrata])
    out <- cbind(G * sapply(OmegaM1epsilon, colSums))
    colnames(out) <- cluster.var
    rownames(out) <- U.cluster[U.indexcluster]
        
    ## ** estimate following random effects
    if(length(nesting.var)>0){
        ls.rho2variable <- attr(nestingStructure, "rho2variable")
        rho2variable <- do.call(rbind, ls.rho2variable)
        rho2variable <- rho2variable[rho2variable$param %in% sapply(ls.tau1, names) == FALSE,,drop=FALSE]

        ## flatten residuals
        vec.OmegaM1epsilon <- unlist(OmegaM1epsilon)
        Vindex.cluster.sorted <- Vindex.cluster[unlist(index.cluster)]

        ## align design matrix with residuals
        X.cor.sorted <- X.cor[unlist(index.cluster),rho2variable$variable,drop=FALSE]
        index.nesting.var <- which(colnames(X.cor.sorted) %in% nesting.var)
        colnames(X.cor.sorted) <- stats::setNames(rho2variable$variable2, rho2variable$variable)[colnames(X.cor.sorted)]

        ## for each level of nesting
        ls.ranef <- lapply(unique(rho2variable$assign), function(iV){ ## iV <- 1
            
            iTable <- rho2variable[rho2variable$assign == iV,,drop=FALSE]
            ## combine X across strata
            iVec <- rowSums(X.cor.sorted[,iTable$variable2,drop=FALSE])
            ## convert to factor and rename
            iDf <- stats::setNames(list(as.factor(iVec)), rho2variable$term.labels2[iV])
            ## expand design matrix relative to each factor level
            iX <- model.matrix(stats::as.formula(paste("~ 0+",names(iDf))), iDf)[,paste0(names(iDf),setdiff(sort(unique(iVec)),0))]
            ## find variance parameters (one for each strata)
            iParam <- stats::setNames(tau[iTable$param], iTable$strata)
            ## expand variance parameters across strata
            iG <- iParam[index.clusterStrata]
            ## compute random effects

            iLs.zranef <- by(sweep(iX, MARGIN = 1, FUN = "*", STATS = vec.OmegaM1epsilon), INDICES = Vindex.cluster.sorted, FUN = colSums, simplify = FALSE)
            iRanef <- sweep(do.call(rbind,iLs.zranef), MARGIN = 1, FUN = "*", STATS = iG)
            colnames(iRanef) <- gsub("_X_XX_X_",":",colnames(iRanef))
            rownames(iRanef) <- U.cluster[U.indexcluster]
            return(iRanef)

        })
        out <- cbind(out,do.call(cbind,ls.ranef))
    }

    ## ** export
    return(out)

}

## * .nestingRanef
##' @description identify nesting of the random effects
##' @noRd
.nestingRanef <- function(object){

    ## ** extract from object
    n.strata <- object$strata$n
    table.rho <- object$design$vcov$param[object$design$vcov$param$type == "rho",]
    name.rho <- table.rho$name
    n.rho <- NROW(table.rho)
    n.cluster <- object$design$cluster$n
    index.cluster <- object$design$index.cluster
    Vindex.cluster <- attr(index.cluster, "vectorwise")
    pattern.cluster <- object$design$vcov$X$pattern.cluster
    X.cor <- object$design$vcov$X$cor
    Xpattern.cor <- object$design$vcov$X$Xpattern.cor    
    name.X.cor <- colnames(X.cor)

    if((object$design$vcov$type != "CS") || object$design$vcov$heterogeneous){
        stop("Identification of random effect structure only implemented for \"CS\" structure with argument heterogenous=FALSE. \n")
    }
    if(any(grepl("_X_XX_X_",colnames(X.cor)))){
        stop("Cannot identify random effects when the name some variables contains \"_X_XX_X_\". \n")
    }

    ## ** Deal with no correlation parameter
    if(NROW(table.rho)==0){
        return(NULL)
    }

    ## ** Find strata (explicit or implicit)
    M.testUnique <- do.call(rbind,by(X.cor,attr(index.cluster, "vectorwise"), function(iX){apply(iX,2,function(iiX){sum(!duplicated(iiX))})}, simplify = FALSE))
    col.strata <- colnames(M.testUnique)[which(colSums(M.testUnique!=1)==0)]
    col.within <- setdiff(colnames(M.testUnique), col.strata)
       
    ## ** deal with nested cases
    ## all clusters
    index.clusterStrata2 <- interaction(as.data.frame(X.cor[sapply(index.cluster,"[",1),col.strata]), drop = TRUE)

    ## relevant subset of clusters
    U.clusterR <- sapply(object$design$vcov$X$Upattern$index.cluster,"[",1)
    strata.clusterR <- index.clusterStrata2[U.clusterR]

    ## *** pattern of parameters for each representative individual
    ls.cluster.design <- lapply(U.clusterR, FUN = function(iC){ ## iC <- 2

        iPattern <- as.character(pattern.cluster[iC,"cor"])
        iIndex.pair <- attr(Xpattern.cor[[iPattern]],"index.pair")
        if(length(col.within)==0){
            iiDF <- data.frame(matrix(0, nrow = 1, ncol = n.rho), NA, NA, NA)
            names(iiDF) <- c(name.rho, "index", "Z", "value")
            iiDF[1,unique(iIndex.pair$param)] <- 1
            return(iiDF)
        }
        iNtime <- NROW(Xpattern.cor[[iPattern]])
        iIndex.obs <- index.cluster[[iC]]
        iZ <- X.cor[iIndex.obs,col.within,drop=FALSE]

        iCol.within <- names(which(colSums(iZ!=0)>0))
        iGrid <- do.call(rbind,lapply(iCol.within, function(iCol){data.frame(col = iCol, value = c(NA,unique(iZ[,iCol.within])))}))

        ls.iG <- lapply(1:NROW(iGrid), function(iiRow){  ##  iiRow <- 1
            iiCol <- iGrid[iiRow,"col"]
            iiValue <- iGrid[iiRow,"value"]
            if(is.na(iiValue)){
                iiZindex.obs <- which(iZ[,iiCol]!=0)
            }else{
                iiZindex.obs <- which(iZ[,iiCol]==iiValue)
            }
            iiZindex.pair <- iIndex.pair[iIndex.pair[,"row"] %in% iiZindex.obs & iIndex.pair[,"col"] %in% iiZindex.obs,]
            iiZindex.pairR <- iiZindex.pair[iiZindex.pair[,"row"]<iiZindex.pair[,"col"],,drop=FALSE]

            iiDF <- data.frame(matrix(0, nrow = 1, ncol = n.rho), NA, iiCol, iiValue)
            names(iiDF) <- c(name.rho, "index", "Z", "value")
            iiDF[1,unique(iiZindex.pairR$param)] <- 1
            iiDF$index <- list(unique(c(iiZindex.pairR[,1],iiZindex.pairR[,2])))
            return(iiDF)
        })

        return(do.call(rbind,ls.iG))

    })

    ## *** aggregate at strata level
    out <- tapply(1:length(strata.clusterR), strata.clusterR, function(iIndex){ ## iIndex <- 1
        iDesign <- do.call(rbind,ls.cluster.design[iIndex])
        iOccurence <- colSums(iDesign[,name.rho,drop=FALSE])
        iOccurence <- iOccurence[iOccurence!=0]
        if(any(duplicated(iOccurence))){
            warning("Non strict nesting of the correlation components - identify of the random effect may be unreliable. \n")
        }
        iNesting <- names(iOccurence)[order(iOccurence, decreasing = TRUE)]
        return(iNesting)         
    }, simplify = FALSE)

    ## *** associate columns of the design matrix to correlation parameters
    df.param <- do.call(rbind,lapply(Xpattern.cor, function(iPattern){
        iIndex.pair <- attr(iPattern,"index.pair")
        iGrid <- do.call(rbind,lapply(c(col.strata,col.within), function(iiVar){ ## iVar <- col.strata[1]
            if(all(iPattern[,iiVar]==0)){
                return(NULL)
            }else{
                iUvalue <- unique(iPattern[,iiVar])
                iTable <- table(iPattern[,iiVar])
                return(data.frame(variable = iiVar, value = iUvalue, rep = as.double(iTable[as.character(iUvalue)])))
            }
        }))
        ## only consider levels with at least two observations, i.e. at least one pair
        iIndex.grid <- which(iGrid$rep>1)
        iGrid.param <- do.call(rbind,lapply(iIndex.grid, function(iiG){ ## iiG <- 1
            iiVar <- iGrid[iiG,"variable"]
            iiValue <- iGrid[iiG,"value"]
            iiPair <- which(iPattern[,iiVar]==iiValue)
            iiVec.param <- iIndex.pair[(iIndex.pair[,"row"] %in% iiPair) & (iIndex.pair[,"col"] %in% iiPair),"param"]
            return(table(factor(iiVec.param, levels = name.rho)))            
        }))
        return(cbind(strata = index.clusterStrata2[attr(iPattern,"index.cluster")[1]], iGrid[iIndex.grid,,drop=FALSE], iGrid.param>0))
    }))
    rownames(df.param) <- NULL

    X.cor.assign <- stats::setNames(attr(X.cor,"assign"), colnames(X.cor))
    X.cor.terms <- stats::setNames(attr(X.cor,"term.labels"), colnames(X.cor))
    
    param2variable <- by(df.param, INDICES = df.param$strata, FUN = function(iDF){ ## iDF <- df.param

        ## rownames(iDF) <- NULL
        iCol <- colSums(iDF[,name.rho,drop=FALSE])
        iCol2 <- iCol[iCol!=0]
        iRow <- stats::setNames(rowSums(iDF[,name.rho,drop=FALSE]), iDF$variable)
        iRow2 <- tapply(iRow, names(iRow), max)

        if(length(iCol2)!=length(iRow2)){
            stop("Something went wrong when associating the correlation parameters to the random effects. \n")
        }
    
        iOut <- data.frame(strata = iDF$strata[1],
                           param = names(sort(iCol2)), ## which parameter is the least there, i.e. is more global as it requires more pairs
                           variable = names(sort(iRow2, decreasing = TRUE)), ## which factor is contains the most parameters, i.e. when TRUE there are the most pairs
                           assign = X.cor.assign[names(iRow2)],
                           term.labels = X.cor.terms[names(iRow2)])
        iOut$variable2 <- gsub(":","_X_XX_X_",iOut$variable, fixed = TRUE) ## as otherwise model.matrix will try to form an interaction instead of taking the variable as it is
        iOut$term.labels2 <- gsub(":","_X_XX_X_",iOut$term.labels)
        return(iOut)
    })
    class(param2variable) <- "list"
    attr(param2variable,"call") <- NULL

    ## ** export
    attr(out,"nested") <- length(col.within)>0
    attr(out,"strata.var") <- col.strata
    attr(out,"nesting.var") <- col.within
    attr(out,"index.clusterStrata") <- index.clusterStrata2
    attr(out,"rho2variable") <- param2variable
    return(out)
}


##----------------------------------------------------------------------
### ranef.R ends here
