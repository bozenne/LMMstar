### ranef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 26 2022 (11:18) 
## Version: 
## Last-Updated: jul  6 2023 (17:42) 
##           By: Brice Ozenne
##     Update #: 233
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * ranef.lmm (documentation)
##' @title Estimate Random Effect From a Linear Mixed Model
##' @description Recover the random effects from the variance-covariance parameter of a linear mixed model.
##' @param object a \code{lmm} object.
##' @param p [numeric vector] value of the model coefficients to be used. Only relevant if differs from the fitted values.
##' @param nestingStructure [list] output of the \code{.nestingRanef} function.
##'
##' @details Consider the following mixed model:
##' \deqn{Y = X\beta + \epsilon}
##' where \eqn{\Sigma_{\epsilon}}, the variance of \eqn{\epsilon}, has a (possibly stratified) compound symmetry structure.
##' Denoting by \eqn{I} the identity matirx, this mean that \eqn{\Sigma_{\epsilon} = \sigma^2 I + Z \Sigma_{\eta} Z^T}
##' where \eqn{\Sigma_{\eta}} is the covariance relative to the design matrix \eqn{Z} (e.g. same student or school). So implicitely we have:
##' \deqn{Y = X\beta + Z \eta + \varepsilon}
##' where \eqn{\varepsilon \sim \mathcal{N}(0, \sigma^2 I)}. So we can estimate the random effets via:
##' \deqn{E[Y|\eta] = Z \Sigma_{\eta} \Omega^{-1} (Y-X\beta)}
##'
##' @keywords methods
##' 
##' @examples
##' data(gastricbypassL, package = "LMMstar")
##' 
##' ## random intercept
##' e.RI <- lmm(weight ~ time + (1|id), data = gastricbypassL)
##' ranef(e.RI)
##' 
##' ## nested random effects
##'
##' ## crossed random effects
##' 

## * ranef.lmm (code)
##' @export
ranef.lmm <- function(object, p = NULL){

    ## ** extract from object
    if(!inherits(object$design$vcov,"RE")){
        stop("Cannot estimate random effects linear mixed models defined by covariance structure (argument \'structure\'). \n",
             "Consider adding random effects in the argument \'formula\' instead. \n")
    }

    ## param
    param.name <- object$design$param$name
    param.type <- stats::setNames(object$design$param$type,param.name)
    param.rho <- param.name[param.type=="rho"]
    param.strata <- unlist(object$design$param[param.type=="rho","index.strata"])

    ## cluster 
    cluster.var <- object$cluster$var
    U.cluster <- object$design$cluster$levels
    attr(cluster.var,"original") <- NULL
    index.cluster <- object$design$index.cluster
    Vindex.cluster <- attr(index.cluster, "vectorwise")
    
    ## strata
    strata.var <- object$design$vcov$name$strata
    indexCluster.strata <- tapply(object$design$vcov$X$Upattern$index.cluster,object$design$vcov$X$Upattern$index.strata,unlist, simplify = FALSE)

    ## design
    X.cor <- object$design$vcov$X$cor
    Xpattern.cor <- object$design$vcov$X$Xpattern.cor
    pattern.cluster <- object$design$vcov$X$pattern.cluster$pattern
    Upattern <- object$design$vcov$X$Upattern
    infoRanef <- object$design$vcov$ranef

    ## missing values
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
    
    ## ** converting correlation parameters into random effect variance
    cumtau <- coef(object, p = p, effects = "correlation", transform.rho = "cov", transform.names = FALSE)
    cumtau.strata <- tapply(cumtau,param.strata,identity, simplify = FALSE)

    n.hierarchy <- length(infoRanef$param)
    index.hierarchy <- unlist(lapply(1:n.hierarchy, function(iH){rep(iH, length(infoRanef$param[[iH]]))}))
    name.RE <- unname(unlist(lapply(infoRanef$param, rownames)))
    n.RE <- length(name.RE)
    varRE <- stats::setNames(rep(NA, n.RE), name.RE)

    for(iH in 1:n.hierarchy){ ## iH <- 1
        iHierarchy <- infoRanef$param[[iH]]
        iCumTau <- apply(iHierarchy, MARGIN = 2, FUN = function(iName){
            cumtau[iName] - head(c(0, cumsum(cumtau[iName])), -1)
        })
        browser()
        varRE[] <- iCumTau
    }
    

    if(infoRanef$type[,"nested"] && infoRanef$type[,"crossed"]){
        stop("Cannot handle nested and crossed random effects. \n")
    }else if(infoRanef$type[,"nested"]){
        rho.strata <- cumrho.strata
        tau.strata <- cumtau.strata
    }else if(infoRanef$type[,"crossed"]){
        rho.strata <- cumrho.strata
        tau.strata <- cumtau.strata
    }else{
        rho.strata <- cumrho.strata
        tau.strata <- cumtau.strata
    }

ranefStructure <- object$design$vcov$ranef
    
    
    browser()

    


    nestingStructure <- .nestingRanef(object)
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



##----------------------------------------------------------------------
### ranef.R ends here
