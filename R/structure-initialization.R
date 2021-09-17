### structure-initialization.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2021 (13:20) 
## Version: 
## Last-Updated: sep 17 2021 (15:14) 
##           By: Brice Ozenne
##     Update #: 50
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * initialization
##' @title Initialize Variance-Covariance Structure
##' @description Initialize the parameters of the variance-covariance structure using residual variance and correlations.
##'
##' @param structure [structure]
##' @param residuals [vector] vector of residuals.
##'
##' @keywords internal
##' 
##' @examples
##' \dontrun{
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$gender <- c("M","F")[as.numeric(gastricbypassL$id) %% 2+1]
##' dd <- gastricbypassL[!duplicated(gastricbypassL[,c("time","gender")]),]
##'
##' eDD.lm <- lm(weight ~ visit, data = dd)
##' eGas.lm <- lm(weight ~ visit*gender, data = gastricbypassL)
##' 
##' ## independence
##' Sid1 <- .skeleton(IND(~1, var.time = "time"), data = dd)
##' Sid4 <- .skeleton(IND(~1|id, var.time = "time"), data = dd)
##' Sdiag1 <- .skeleton(IND(~visit), data = dd)
##' Sdiag4 <- .skeleton(IND(~visit|id), data = dd)
##' Sdiag24 <- .skeleton(IND(~visit+gender|id, var.time = "time"), data = gastricbypassL)
##'
##' .initialize(Sid1, residuals = residuals(eDD.lm))
##' ## sd(residuals(eDD.lm))
##' .initialize(Sdiag4, residuals = residuals(eDD.lm))
##' ## tapply(residuals(eDD.lm),dd$visit,sd)
##' .initialize(Sdiag24, residuals = residuals(eGas.lm))
##' ## tapply(residuals(eGas.lm),interaction(gastricbypassL[,c("visit","gender")]),sd)
##' 
##' ## compound symmetry
##' Scs4 <- .skeleton(CS(~1|id, var.time = "time"), data = gastricbypassL)
##' Scs24 <- .skeleton(CS(gender~time|id), data = gastricbypassL)
##' 
##' .initialize(Scs4, residuals = residuals(eGas.lm))
##' 
##' ## unstructured
##' Sun4 <- skeleton(UN(~visit|id), data = gastricbypassL)
##' param4 <- setNames(c(1,1.1,1.2,1.3,0.5,0.45,0.55,0.7,0.1,0.2),Sun4$param$name)
##' Sun24 <- skeleton(UN(gender~visit|id), data = gastricbypassL)
##' param24 <- setNames(c(param4,param4*1.1),Sun24$param$name)
##' 
##' .calc_dOmega(Sun4, param = param4)
##' .calc_dOmega(Sun24, param = param24)
##' }
`.initialize` <-
    function(object, residuals, p, ssc) UseMethod(".initialize")

## * initialization.IND
.initialize.IND <- function(object, residuals, p = 1, ssc = TRUE){

    param.type <- setNames(object$param$type,object$param$name)
    param.strata <- setNames(object$param$strata,object$param$name)
    Upattern.name <- object$X$Upattern$name
    cluster.pattern <- object$X$cluster.pattern

    ## combine all residuals and all design matrices
    M.res <- do.call(rbind,lapply(object$X$var, function(iPattern){ ## iPattern <- object$X$var[[1]]
        X.iPattern <- iPattern
        cluster.iPattern <- attr(iPattern,"index.cluster")
        obs.iPattern <- unlist(attr(iPattern,"index.obs"))
        attr(X.iPattern,"indicator.param") <- NULL
        attr(X.iPattern,"index.cluster") <- NULL
        attr(X.iPattern,"index.obs") <- NULL
        iOut <- cbind(residuals = residuals[obs.iPattern], do.call(rbind,rep(list(X.iPattern),length(cluster.iPattern))))
        return(iOut)
    }))

    ## extract information
    epsilon2 <- M.res[,1]^2
    X <- M.res[,-1,drop=FALSE]
    paramVar.type <- param.type[colnames(X)]
    paramVar.strata <- param.strata[colnames(X)]
    n.strata <- length(unique(paramVar.strata))
    n.obs <- NROW(X)

    ## small sample correction (inflate residuals)
    if(ssc){
        UX <- interaction(as.data.frame(X), drop = TRUE)
        n.UX <- table(UX)
        epsilon2.ssc <- epsilon2 * (n.UX/(n.UX-p))[UX]
    }else{
        epsilon2.ssc <- epsilon2
    }

    ## fit
    e.res <- lm.fit(y=epsilon2.ssc,x=X) 
    if(all(paramVar.type=="sigma")){
        out <- sqrt(e.res$coef)
    }else{
        ## try to move from additive to full interaction model
        ls.Z <- lapply(1:n.strata, function(iStrata){
            iParamVar.type <- paramVar.type[paramVar.strata == iStrata]
            iX <- X[,paramVar.strata==iStrata,drop=FALSE]
            if(any("k" %in% iParamVar.type)){
                iX[,iParamVar.type=="sigma"] <- iX[,iParamVar.type=="sigma"] - rowSums(iX[,iParamVar.type!="sigma"])
            }
            return(iX)
        })
        
        Z <- do.call(cbind,ls.Z)[,colnames(X)]
        eTest.res <- lm.fit(y=epsilon2.ssc,x=Z)

        if(all(abs(e.res$fitted.value-eTest.res$fitted.value)<1e-8)){
            ls.out <- lapply(1:n.strata, function(iStrata){
                iParamVar.type <- paramVar.type[paramVar.strata == iStrata]
                iOut <- sqrt(eTest.res$coef[names(iParamVar.type)])
                if(any("k" %in% iParamVar.type)){
                    iOut[iParamVar.type=="k"] <- iOut[iParamVar.type=="k"]/iOut[iParamVar.type=="sigma"]
                }
                return(iOut)
            })
            out <- unlist(ls.out)[names(paramVar.type)]
        }else{ ## failure of the full interaction model. Use a log transform
            e.res <- lm.fit(y=log(epsilon2*n.obs/(n.obs-p)),x=X)
            out <- exp(0.5*e.res$coef)
        }
    }

    ##  standardize residuals
    if(identical(attr(residuals,"studentized"),TRUE)){
        attr(residuals,"studentized") <- NULL
        attr(out,"studentized") <- residuals/exp(X %*% log(out))
    }

    ## export
    return(out)
}

## * initialize.CS
.initialize.CS <- function(object, residuals, p = 1, ssc = TRUE){
    ## extract information
    param.type <- setNames(object$param$type,object$param$name)
    param.strata <- setNames(object$param$strata,object$param$name)
    Upattern.name <- object$X$Upattern$name
    cluster.pattern <- object$X$cluster.pattern

    ## estimate variance and standardize residuals
    attr(residuals,"studentized") <- TRUE ## to return studentized residuals
    sigma <- .initialize.IND(object = object, residuals = residuals, p = p, ssc = ssc)
    residuals.studentized <- attr(sigma, "studentized")
    attr(sigma, "studentized") <- NULL

    ## combine all residuals and all design matrices
    M.res <- do.call(rbind,lapply(1:length(object$X$cor), function(iPattern){ ## iPattern <- 1
        browser()
        X.iPattern <- object$X$cor[[iPattern]]
        cluster.iPattern <- attr(X.iPattern,"index.cluster")
        obs.iPattern <- do.call(rbind,attr(X.iPattern,"index.obs"))
        iIndex.pairtime <- attr(X.iPattern,"index.pairtime")

        obs.iPattern[,iIndex.pairtime[1,]]
        pairobs.iPattern <- lapply(obs.iPattern, function(iC){
            cbind(obs.iPattern[[iC]][],obs.iPattern[iIndex.pairtime[2,]])
        })
        
        attr(X.iPattern,"index.vec2matrix") <- NULL
        attr(X.iPattern,"indicator.param") <- NULL
        attr(X.iPattern,"index.strata") <- NULL
        attr(X.iPattern,"index.pairtime") <- NULL
        attr(X.iPattern,"index.Utime") <- NULL
        attr(X.iPattern,"index.cluster") <- NULL
        attr(X.iPattern,"index.obs") <- NULL
        iOut <- cbind(residuals1 = residuals[obs.iPattern], residuals1 = residuals[obs.iPattern],
                      do.call(rbind,rep(list(X.iPattern),length(cluster.iPattern))))
        return(iOut)
    }))

}

## * initialize.CS
.initialize.UN <- function(object, residuals){
}



##----------------------------------------------------------------------
### structure-initialization.R ends here
