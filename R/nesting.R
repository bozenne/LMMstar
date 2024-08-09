### nesting.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 15 2024 (11:57) 
## Version: 
## Last-Updated: aug  8 2024 (11:24) 
##           By: Brice Ozenne
##     Update #: 149
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * checkNesting
##' @title Check Nesting between Linear Mixed Models
##' @description Check whether one linear mixed model has its mean, variance, and correlation structure nested in the other.
##'
##' @param objectH0 [lmm]
##' @param objectH1 [lmm]
##' @param sep [character] character used to combine covariate values.
##' @param tol [numeric] tolerance when comparing mean and variance structure.
##' 
##' @noRd
.checkNesting <- function(objectH0, objectH1, sep = ".", tol = 1e-8){
   
    ## ** number of observations
    nobsH0 <- stats::nobs(objectH0)
    nobsH1 <- stats::nobs(objectH1)
    if(any(nobsH0 != nobsH1)){
        if(nobsH0["missing"]!=nobsH1["missing"]){
            stop("Mismatch between the number of observations between the two models - could be due to missing data. \n",
                 "H0: ",paste(paste(names(nobsH0),"=",nobsH0), collapse = ", "),".\n",
                 "H1: ",paste(paste(names(nobsH1),"=",nobsH0), collapse = ", "),".\n")
        }else{
            stop("Mismatch between the number of observations between the two models. \n",
                 "H0: ",paste(paste(names(nobsH0),"=",nobsH0), collapse = ", "),".\n",
                 "H1: ",paste(paste(names(nobsH1),"=",nobsH0), collapse = ", "),".\n")
        }
    }

    ## ** check outcome
    if(any(abs(objectH1$design$Y-objectH0$design$Y)>tol)){
        stop("Mismatch in outcome between the two models. \n")
    }

    ## ** extract parameters
    table.paramH0 <- stats::model.tables(objectH0, effects = "param")
    table.paramH1 <- stats::model.tables(objectH1, effects = "param")
    
    name.paramH0 <- table.paramH0$trans.name
    name.paramH1 <- table.paramH1$trans.name
    type.paramH0 <- table.paramH0$type
    type.paramH1 <- table.paramH1$type
    if(any(table(factor(type.paramH1, levels = c("mu","sigma","k","rho"))) < table(factor(type.paramH0, levels = c("mu","sigma","k","rho"))))){
        ## table check
        tableType.paramH0 <- table(factor(type.paramH0, levels = c("mu","sigma","k","rho")))
        tableType.paramH1 <- table(factor(type.paramH1, levels = c("mu","sigma","k","rho")))
        stop("Cannot perform a likelihood ratio test when the model are not nested. \n",
             "The model with the largest log-likelihood has fewer parameters of type \"",paste(names(tableType.paramH1)[tableType.paramH1<tableType.paramH0],collapse ="\", \""),"\". \n")
    }
    mismatchH0 <- stats::setNames(name.paramH0 %in% name.paramH1 == FALSE, name.paramH0)
    mismatchH1 <- stats::setNames(name.paramH1 %in% name.paramH0 == FALSE, name.paramH1)
    rhs <- stats::setNames(rep("0", sum(mismatchH1)), names(which(mismatchH1)))
    current.mismatchH0 <- mismatchH0
    current.mismatchH1 <- mismatchH1

    ## ** mean structure
    rhs <- NULL
    mu.paramH0 <- name.paramH0[type.paramH0=="mu"]
    mu.paramH1 <- name.paramH1[type.paramH1=="mu"]
    if(length(mu.paramH0)==length(mu.paramH1)){
        if(all(sort(mu.paramH0)==sort(mu.paramH1))){
            ## same parameter names
            equal.mean <- TRUE
        }else{
            ## same mean fit
            meanfitH0 <- stats::lm.fit(x = objectH0$design$mean, y = objectH0$design$Y)$fitted.values
            meanfitH1 <- stats::lm.fit(x = objectH1$design$mean, y = objectH1$design$Y)$fitted.values
            if(all(abs(meanfitH1-meanfitH0)<tol)){
                equal.mean <- TRUE
            }else{
                stop("Cannot perform a likelihood ratio test when the model are not nested. \n",
                     "The two models have the same number of mean parameters but do not seem to induce a similar mean structure (discrepancy ",max(abs(meanfitH1-meanfitH0)),"). \n")
            }
        }
    }else if(length(mu.paramH0)<length(mu.paramH1)){
        if(all(mu.paramH0 %in% mu.paramH1)){
            ## common parameter names + extra
            equal.mean <- FALSE
            rhs <- c(rhs, stats::setNames(rep(0, length(mu.paramH1)-length(mu.paramH0)), setdiff(mu.paramH1,mu.paramH0)))
        }else{
             ## at least as good mean fit            
            meanfitH0 <- stats::lm.fit(x = objectH0$design$mean, y = objectH0$design$Y)$fitted.values
            ## Since space spanned by X.H0 is a subspace of the one spanned by X.H1  should be able to perfectly fit 
            meanresH1 <- stats::lm.fit(x = objectH1$design$mean, y = meanfitH0)$residuals
            rhs <- c(rhs, stats::setNames("?[mean]?","?"))
            if(all(abs(meanresH1)<tol)){
                equal.mean <- FALSE
            }else{
                stop("Cannot perform a likelihood ratio test when the model are not nested. \n",
                     "The mean structure of the model with the larger likelihood does not include the mean model of the other lmm (discrepancy ",max(abs(meanresH1)),"). \n")
            }
            
        }
    } ## third case (length(mu.paramH0)>length(mu.paramH1)) already dealt with in the table check


    ## ** variance-covariance structure
    structure.H0 <- objectH0$design$vcov$class
    structure.H1 <- objectH0$design$vcov$class
    test.sameStructure <- (structure.H0 == structure.H1) || (structure.H0 %in% c("ID","IND","CS","UN") && structure.H1 %in% c("ID","IND","CS","UN"))

    strata.H0 <- stats::na.omit(objectH0$design$vcov$name$strata)
    strata.H1 <- stats::na.omit(objectH1$design$vcov$name$strata)

    ## *** variance structure
    var.H0 <- stats::na.omit(unlist(objectH0$design$vcov$name$var))
    var.H1 <- stats::na.omit(unlist(objectH1$design$vcov$name$var))
    sigmak.paramH0 <- name.paramH0[type.paramH0 %in% c("sigma","k")]
    sigmak.paramH1 <- name.paramH1[type.paramH1 %in% c("sigma","k")]

    if(length(sigmak.paramH0)==length(sigmak.paramH1)){
        if(all(sort(sigmak.paramH0)==sort(sigmak.paramH1))){
            ## same parameter names
            equal.var <- TRUE
        }else if(test.sameStructure){
            ## same mean fit so probably same var fit (TRUE for ID,IND,CS,UN likely for CUSTOM)
            varfitH0 <- stats::lm.fit(x = objectH0$design$vcov$var$X, y = objectH0$design$Y)$fitted.values
            varfitH1 <- stats::lm.fit(x = objectH1$design$vcov$var$X, y = objectH1$design$Y)$fitted.values
            if(all(abs(varfitH1-varfitH0)<tol)){
                equal.var <- TRUE
            }else{
                stop("Cannot perform a likelihood ratio test when the model are not nested. \n",
                     "The two models have the same number of variance parameters but do not seem to induce a similar variance structure (discrepancy ",max(abs(meanfitH1-meanfitH0)),"). \n")
            }
        }else{
            stop("Cannot decided about the nesting between the two models when a ",setdiff(union(structure.H0,structure.H1),c("ID","IND","CS","UN"))[1]," covariance structure is used. \n")
        }
    }else if(length(sigmak.paramH0)<length(sigmak.paramH1)){

        if(all(sigmak.paramH0 %in% sigmak.paramH1) && test.sameStructure){
            ## common parameter names + extra
            equal.var <- FALSE
            rhs.var <- stats::setNames(rep(1, length(sigmak.paramH1)-length(sigmak.paramH0)), setdiff(sigmak.paramH1,sigmak.paramH0))
            rhs.var[type.paramH1[setdiff(sigmak.paramH1,sigmak.paramH0)]=="sigma"] <- NA
            rhs <- c(rhs, rhs.var)
        }else if(test.sameStructure && length(setdiff(var.H0, var.H1))==0 && length(setdiff(strata.H0, strata.H1))==0){
            ## common variable names + extra
            equal.var <- FALSE
            extra.var <- c(setdiff(var.H1, var.H0),setdiff(strata.H1, strata.H0))

            if(length(var.H0) == 0 && length(strata.H0) == 0){
                ## no covariate for the reference model, only sigma parameter
                rhs.var <- NULL
                if(sum(type.paramH1[sigmak.paramH1]=="sigma")>1){
                    rhs.var <- c(rhs.var,
                                 stats::setNames(sigmak.paramH1[type.paramH1[sigmak.paramH1]=="sigma"][1], paste(sigmak.paramH1[type.paramH1[sigmak.paramH1]=="sigma"][-1], collapse = "=="))
                                 )
                }
                if(sum(type.paramH1[sigmak.paramH1]=="k")>0){                
                    rhs.var <- c(rhs.var,
                                 stats::setNames(rep(1, sum(type.paramH1[sigmak.paramH1]=="k")),sigmak.paramH1[type.paramH1[sigmak.paramH1]=="k"])
                                 )
                }

            }else if(identical(sort(var.H1),sort(var.H0))){
                ## stratified structure
                n.strata <- objectH1$strata$n
                sigmak.strataH1 <- stats::setNames(objectH1$design$param$index.strata,objectH1$design$param$name)[sigmak.paramH1]
                rhs.var <- stats::setNames(rep(NA, sum(sigmak.strataH1!=1)),sigmak.paramH1[sigmak.strataH1!=1])
                if(sum(sigmak.strataH1==1)*(n.strata-1) == sum(sigmak.strataH1!=1)){
                    rhs.var[] <- rep(sigmak.paramH1[sigmak.strataH1==1], n.strata-1)
                }                
            }else if(length(extra.var)>0){

                ## find levels not matching the reference
                Mvar.level <- attr(objectH1$design$vcov$var$X,"M.level")
                ls.testref <- lapply(extra.var, function(iVar){Mvar.level[[iVar]]!=objectH1$xfactor$var[[iVar]][1]})
                name.vartest <- colnames(objectH1$design$vcov$var$X)[Reduce("+",ls.testref)>0]
                rhs.var <- stats::setNames(rep(NA, length(name.vartest)), name.vartest)

                ## find matching coefficient among common covariates that are not being tested
                Mvar.level2 <- nlme::collapse(Mvar.level[,setdiff(colnames(Mvar.level),extra.var),drop=FALSE], sep = sep[1])
                rhs.var[] <- sapply(match(names(rhs.var),colnames(objectH1$design$vcov$var$X)), function(iIndexH1){
                    iIndexH0 <- setdiff(which(Mvar.level2[iIndexH1]==Mvar.level2),iIndexH1)
                    if(length(iIndexH0)==1){
                        return(colnames(objectH1$design$vcov$var$X)[iIndexH0])
                    }else{
                        return(NA)
                    }
                })
            }else{
                rhs.var <- stats::setNames("?","?[var]?")
            }
            rhs <- c(rhs, rhs.var)

        }else {
            ## at least as good mean fit so probably at least as god var fit (TRUE for ID,IND,CS,UN likely for CUSTOM)
            varfitH0 <- stats::lm.fit(x = objectH0$design$vcov$var$X, y = objectH0$design$Y)$fitted.values
            ## Since space spanned by X.H0 is a subspace of the one spanned by X.H1  should be able to perfectly fit 
            varresH1 <- stats::lm.fit(x = objectH1$design$vcov$var$X, y = varfitH0)$residuals
            
            rhs <- c(rhs, stats::setNames("?","?[var]?"))

            if(all(abs(varresH1)<tol)){
                equal.var <- FALSE
            }else{
                stop("Cannot perform a likelihood ratio test when the model are not nested. \n",
                     "The variance structure of the model with the larger likelihood does not include the variance model of the other lmm (discrepancy ",max(abs(varresH1)),"). \n")
            }
            
        }
    } ## third case (length(sigmak.paramH0)>length(sigmak.paramH1)) already dealt with in the table check

    ## *** correlation structure
    if(structure.H0 != structure.H1 && (structure.H0 %in% c("TOEPLITZ","CUSTOM") || structure.H0 %in% c("TOEPLITZ","CUSTOM"))){
        stop("Cannot decided about the nesting between the two models when a ",setdiff(union(structure.H0,structure.H1),c("ID","IND","CS","UN"))[1]," covariance structure is used. \n")
    }

    cor.H0 <- stats::na.omit(unlist(objectH0$design$vcov$name$cor))
    cor.H1 <- stats::na.omit(unlist(objectH1$design$vcov$name$cor))
    rho.paramH0 <- name.paramH0[type.paramH0 %in% "rho"]
    rho.paramH1 <- name.paramH1[type.paramH1 %in% "rho"]
        
    if(length(rho.paramH0)==length(rho.paramH1)){
        if(all(sort(rho.paramH0)==sort(rho.paramH1))){
            ## same parameter names
            equal.cor <- TRUE
        }else if(test.sameStructure){

            table.cor <- Reduce("+",lapply(1:objectH0$cluster$n, function(iC){ ## iC <- 1
                iPatternH0 <- objectH0$design$vcov$cor$pattern[[iC]]
                iPatternH1 <- objectH1$design$vcov$cor$pattern[[iC]]
                table(factor(attr(objectH0$design$vcov$cor$Xpattern[[iPatternH0]],"index.pair")$param, levels = rho.paramH0),
                      factor(attr(objectH1$design$vcov$cor$Xpattern[[iPatternH1]],"index.pair")$param, levels = rho.paramH1))
            }))

            if(all(rowSums(table.cor!=0)==1) && all(colSums(table.cor!=0)==1)){
                equal.cor <- TRUE
            }else{
                stop("Cannot perform a likelihood ratio test when the model are not nested. \n",
                     "The two models have the same number of corelation parameters but do not seem to induce a similar correlation structure. \n")
            }
        } ## third case (!test.sameStructure) already delta with just below *** corelation structure
    }else if(length(rho.paramH0)<length(rho.paramH1)){

        if(all(rho.paramH0 %in% rho.paramH1) && test.sameStructure){
            ## common parameter names + extra
            equal.cor <- FALSE
            rhs <- c(rhs, stats::setNames(rep(0, length(rho.paramH1)-length(rho.paramH0)), setdiff(rho.paramH1,rho.paramH0)))

        }else if(test.sameStructure && length(setdiff(cor.H0, cor.H1))==0 && length(setdiff(strata.H0, strata.H1))==0){
            ## common covariate/strata names + extra
            equal.cor <- FALSE

            extra.var <- c(setdiff(cor.H1, cor.H0),setdiff(strata.H1, strata.H0))
            if(identical(sort(cor.H1),sort(cor.H0))){
                ## stratified structure
                n.strata <- objectH1$strata$n
                rho.strataH1 <- stats::setNames(objectH1$design$param$index.strata,objectH1$design$param$name)[rho.paramH1]
                rhs.rho <- stats::setNames(rep(NA, times = sum(rho.strataH1!=1)),rho.paramH1[rho.strataH1!=1])
                if(sum(rho.strataH1==1)*(n.strata-1) == sum(rho.strataH1!=1)){
                    rhs.rho[] <- rep(rho.paramH1[rho.strataH1==1], n.strata-1)
                }
            }else if(length(extra.var)>0){                
                ## additional covariate possibly also strata

                ## find levels not matching the reference
                ls.testref <- lapply(extra.var, function(iVar){attr(objectH1$design$vcov$cor$X,"M.level")[[iVar]]!=objectH1$xfactor$cor[[iVar]][1]})
                name.cortest <- colnames(objectH1$design$vcov$cor$X)[Reduce("+",ls.testref)>0]

                ## find 'extra' correlation parameters as those who correspond to not the reference level
                ls.rhotest <- lapply(objectH1$design$vcov$cor$Xpattern, function(iX){
                    iX.cortest <- iX[,name.cortest,drop=FALSE]
                    if(sum(abs(iX.cortest))==0){
                        return(NULL)
                    }else{
                        iIndex <- which(rowSums(abs(iX.cortest))>0)
                        iRowRho <- attr(iX,"index.pair")[attr(iX,"index.pair")$row %in% iIndex,"param"]
                        iColRho <- attr(iX,"index.pair")[attr(iX,"index.pair")$col %in% iIndex,"param"]
                        return(union(iRowRho,iColRho))                        
                    }
                })
                rhs.rho <- unique(unlist(ls.rhotest))

                ## find the 'orginal' correlation parameters (i.e. those from H0)
                if(length(setdiff(rho.paramH1, rhs.rho))==1){
                    ## all vs. 1
                    names(rhs.rho) <- stats::setNames(rep(setdiff(rho.paramH1, rhs.rho),length(rhs.rho)),rhs.rho)
                }else if(length(setdiff(rho.paramH1,rhs.rho))==0){
                    ## all equal (within strata)
                    strata.rho <- stats::setNames(table.paramH1$index.strata,table.paramH1$name)[rhs.rho]
                    rhs.rho <- unlist(unname(tapply(names(strata.rho),strata.rho,function(iVec){ ## iVec <- names(strata.rho)[strata.rho==1]
                        stats::setNames(iVec[1],paste(iVec[-1],collapse="=="))
                    }, simplify = FALSE)))
                }else{
                    ## specific equality (rho:male = rho:female)
                    Mcor.H02H1 <- .corH02H1(objectH0, objectH1)
                    names(rhs.rho) <- sapply(rhs.rho, function(iName){
                        paste(names(which(rowSums(Mcor.H02H1[setdiff(rho.paramH1,rhs.rho),which(Mcor.H02H1[iName,]),drop=FALSE])>0)), collapse = ", ")
                    })
                    
                }
                
            }else{
                rhs.rho <- stats::setNames("?", "?[cor]?")
            }
            rhs <- c(rhs, rhs.rho)

        }else{            
            Mcor.H02H1 <- .corH02H1(objectH0, objectH1)
            name.rhoref <- apply(Mcor.H02H1, MARGIN = 2, FUN = function(iCol){names(which(iCol))[1]})

            if(all(rowSums(Mcor.H02H1)==1) && all(!duplicated(name.rhoref))){
                equal.cor <- FALSE

                rhs.rho <- setdiff(rho.paramH1,name.rhoref)
                names(rhs.rho) <- sapply(rhs.rho, function(iName){
                    paste(names(which(rowSums(Mcor.H02H1[name.rhoref,which(Mcor.H02H1[iName,]),drop=FALSE])>0)), collapse = ", ")
                })
                rhs <- c(rhs, rhs.rho)                
            }else{
                stop("Cannot perform a likelihood ratio test when the model are not nested. \n",
                     "The correlation structure of the model with the larger likelihood does not include the correlation model of the lmm. \n")
            }
            
        }
    } ## third case (length(rho.paramH0)>length(rho.paramH1)) already dealt with in the table check

    ## ** export
    out <- c(mean = equal.mean, var = equal.var, cor = equal.cor)
    attr(out,"rhs") <- rhs
    return(out)
}

## * .corH02H1
##' @title Map correlation from H0 to H1
##' @description Map correlation parameters from H0 to H1
##'
##' @return a matrix with TRUE/FALSE for which H1 correlation parameters correspond to which H0 correlation parameters.
##' @noRd
.corH02H1  <- function(objectH0, objectH1){

    ## ** extract from object
    table.paramH0 <- stats::model.tables(objectH0, effects = "param")
    table.paramH1 <- stats::model.tables(objectH1, effects = "param")
    
    name.paramH0 <- table.paramH0$trans.name
    name.paramH1 <- table.paramH1$trans.name
    type.paramH0 <- table.paramH0$type
    type.paramH1 <- table.paramH1$type
    rho.paramH0 <- name.paramH0[type.paramH0 %in% "rho"]
    rho.paramH1 <- name.paramH1[type.paramH1 %in% "rho"]

    ## ** table H0-H1
    table.cor <- Reduce("+",lapply(1:objectH0$cluster$n, function(iC){ ## iC <- 1
        iPatternH0 <- objectH0$design$vcov$cor$pattern[[iC]]
        iPatternH1 <- objectH1$design$vcov$cor$pattern[[iC]]
        table(factor(attr(objectH0$design$vcov$cor$Xpattern[[iPatternH0]],"index.pair")$param, levels = rho.paramH0),
              factor(attr(objectH1$design$vcov$cor$Xpattern[[iPatternH1]],"index.pair")$param, levels = rho.paramH1))
    }))

    ## ** matrix H0-H1
    M.out <- do.call(cbind,apply(table.cor, MARGIN = 1, FUN = function(iRow){
        rho.paramH1 %in% names(which(iRow!=0))
    }, simplify = FALSE))
    rownames(M.out) <- rho.paramH1

    ## ** export
    return(M.out)
}

##----------------------------------------------------------------------
### nesting.R ends here

