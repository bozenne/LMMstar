### ranef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 26 2022 (11:18) 
## Version: 
## Last-Updated: jul  7 2023 (18:16) 
##           By: Brice Ozenne
##     Update #: 264
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
    var.cluster <- object$cluster$var
    n.cluster <- object$design$cluster$n
    index.cluster <- object$design$index.cluster
    Vindex.cluster <- attr(index.cluster, "vectorwise")
    
    ## strata
    var.strata <- object$design$vcov$name$strata
    index.clusterStrata <- object$design$index.clusterStrata
    U.strata <- object$strata$levels
    n.strata <- length(U.strata)

    ## design
    X.cor <- object$design$vcov$X$cor
    Xpattern.cor <- object$design$vcov$X$Xpattern.cor
    pattern.cluster <- object$design$vcov$X$pattern.cluster$pattern
    Upattern <- object$design$vcov$X$Upattern
    infoRanef <- object$design$vcov$ranef
    name.hierarchy <- unlist(infoRanef$hierarchy, use.names = FALSE)

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
    varRE <- matrix(NA, nrow = n.RE, ncol = n.strata,
                    dimnames =  list(name.RE, U.strata))
    for(iH in 1:n.hierarchy){ ## iH <- 1
        iHierarchy <- infoRanef$param[[iH]]
        varRE[rownames(iHierarchy),] <- apply(iHierarchy, MARGIN = 2, FUN = function(iName){
            cumtau[iName] - c(0,utils::head(cumtau[iName],-1))
        })
    }
    
    if(any(varRE<=0)){
        stop("Variance for the random effects is found to be negative - cannot estimate the random effects. \n")
    }
    
    ## ** extract raw residuals
    df.epsilon <- stats::residuals(object, p = p, keep.data = TRUE, type = "response", format = "long")
    if(length(index.na)>0){
        df.epsilon <- df.epsilon[-index.na,,drop=FALSE]
    }
    df.epsilon.order <- df.epsilon[order(df.epsilon$XXcluster.indexXX),,drop=FALSE]
    
    ## ** extract inverse residual variance-covariance matrix
    ls.OmegaM1 <- stats::sigma(object, p = p, cluster = unique(df.epsilon$XXcluster.indexXX), inverse = TRUE, simplifies = FALSE)
    ls.epsilon <- base::tapply(df.epsilon$r.response,df.epsilon$XXcluster.indexXX,list)## split residuals by id

    ## ** design matrix
    OmegaM1epsilon <- mapply(x = ls.OmegaM1 , y = ls.epsilon, FUN = `%*%`, SIMPLIFY = FALSE)

    ## ** estimate random effects
        
    ## flatten residuals
    vec.OmegaM1epsilon <- unlist(OmegaM1epsilon)
    out <- stats::setNames(vector(mode = "list", length = n.hierarchy),
                           sapply(infoRanef$hierarchy,paste,collapse=":"))

    if(infoRanef$crossed == FALSE && infoRanef$nested == FALSE){
        name.hierarchy <- infoRanef$hierarchy[[1]]
        out[[1]] <- tapply(vec.OmegaM1epsilon, df.epsilon.order[[name.hierarchy]], sum) * varRE[,tapply(df.epsilon.order$XXstrataXX, df.epsilon.order[[name.hierarchy]], "[", 1)]
        out[[1]] <- out[[1]][unique(df.epsilon[[name.hierarchy]])] ## re-order according to original data
        
    }else{
        X.corA <- cbind(attr(index.cluster,"vectorwise"), X.cor)
        colnames(X.corA)[1] <- var.cluster

        if(infoRanef$nested == FALSE){
            ## variance-covariance of the random effects per strata
            ls.Tau <- apply(varRE, MARGIN = 2, diag, nrow = n.hierarchy, simplify = FALSE)        
        }else{
        }
    
    browser()
    lapply(1:n.cluster, function(iC){ ## iC <- 1
        iTau <- ls.Tau[[index.clusterStrata[iC]]]
        iX <- X.corA[index.cluster[[iC]], name.hierarchy,drop=FALSE]
        iZ1 <- model.matrix(~0+plate, data = data.frame(plate = as.factor(iX[,"plate"])))
        iZ2 <- model.matrix(~0+sample, data = data.frame(sample = as.factor(iX[,"sample"])))

         t(iZ1) %*% vec.OmegaM1epsilon[index.cluster[[iC]]] / GS$plate
         t(iZ2) %*% vec.OmegaM1epsilon[index.cluster[[iC]]] / GS$sample
    })
    ## tapply(vec.OmegaM1epsilon, df.epsilon.order$XXcluster.indexXX, function(iEp))

    for(iH in 1:n.hierarchy){ ## iH <- 1
        iVarH <- 
        if(length(iVarH)==1){
            browser()
            

            tapply(unlist(ls.epsilon), df.epsilon.order[[iVarH]], mean) / coef(object, effects = "variance")
            
            
        }else if(length(iVarH)>1){
            browser()
        }
    }
    }

    ## ** export
    return(out)

}



##----------------------------------------------------------------------
### ranef.R ends here
