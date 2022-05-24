### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: maj 24 2022 (21:03) 
##           By: Brice Ozenne
##     Update #: 488
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * coef.lmm (documentation)
##' @title Extract Coefficients From a Linear Mixed Model
##' @description Extract coefficients from a linear mixed model.
##' @name coef
##'
##' @param object a \code{lmm} object.
##' @param effects [character] Should all coefficients be output (\code{"all"}),
##' or only coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only coefficients relative to the variance structure (\code{"variance"}),
##' or only coefficients relative to the correlation structure (\code{"correlation"}).
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param p [numeric vector] value of the model coefficients to be used. Only relevant if differs from the fitted values.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##' 
##'
##' @details \bold{transform.sigma}: \cr
##' \itemize{
##' \item \code{"none"} ouput residual standard error.
##' \item \code{"log"} ouput log-transformed residual standard error.
##' \item \code{"square"} ouput residual variance.
##' \item \code{"logsquare"} ouput log-transformed residual variance.
##' }
##'
##'  \bold{transform.k}: \cr
##' \itemize{
##' \item \code{"none"} ouput ratio between the residual standard error of the current level and the reference level.
##' \item \code{"log"} ouput log-transformed ratio between the residual standard errors.
##' \item \code{"square"} ouput ratio between the residual variances.
##' \item \code{"logsquare"} ouput log-transformed ratio between the residual variances.
##' \item \code{"sd"} ouput residual standard error of the current level.
##' \item \code{"logsd"} ouput residual log-transformed standard error of the current level.
##' \item \code{"var"} ouput residual variance of the current level.
##' \item \code{"logvar"} ouput residual log-transformed variance of the current level.
##' }
##' 
##'  \bold{transform.rho}: \cr
##' \itemize{
##' \item \code{"none"} ouput correlation coefficient.
##' \item \code{"atanh"} ouput correlation coefficient after tangent hyperbolic transformation.
##' \item \code{"cov"} ouput covariance coefficient.
##' }
##'
##' When using a (pure) compound symmetry covariance structure (\code{structure = "CS"}),
##' estimated random effects can be extracted by setting argument \code{effects} to \code{"ranef"}.
##'
##' @return A vector with the value of the model coefficients.
##' 
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit linear mixed model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##'
##' ## output coefficients
##' coef(eUN.lmm)
##' coef(eUN.lmm, effects = "mean")
##' coef(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none")

## * coef.lmm (code)
##' @rdname coef
##' @export
coef.lmm <- function(object, effects = NULL, p = NULL,
                     transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE, ...){

    ## ** extract from object
    param.name <- object$design$param$name
    param.type <- stats::setNames(object$design$param$type,param.name)
    param.level <- stats::setNames(object$design$param$level,param.name)
    param.sigma <- stats::setNames(object$design$param$sigma,param.name)
    param.strata <- stats::setNames(object$design$param$strata,param.name)
    param.k.x <- stats::setNames(object$design$param$k.x,param.name)
    param.k.y <- stats::setNames(object$design$param$k.y,param.name)

    object.reparametrize.name <- names(object$reparametrize$p)
    object.reparametrize.value <- object$reparametrize$p
    object.reparametrize.newname <- object$reparametrize$newname

    index.na <- object$index.na
    type.pattern <- object$design$vcov$type
    
    U.strata <- object$strata$levels
    strata.var <- object$strata$var
    n.strata <- object$strata$n
    U.cluster.original <- object$design$cluster$levels.original
    cluster.var <- object$cluster$var
    n.cluster <- object$cluster$n
    X.cor <- object$design$vcov$X$cor
    Xpattern.cor <- object$design$vcov$X$Xpattern.cor
    index.cluster <- object$design$index.cluster
    pattern.cluster <- object$design$vcov$X$pattern.cluster$pattern
    Upattern <- object$design$vcov$X$Upattern

    ## ** normalize user imput
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    
    if(is.null(effects)){
        if(transform.sigma == "none" && transform.k == "none" && transform.rho == "none"){
            effects <- options$effects
        }else{
            effects <- c("mean","variance","correlation")
        }
    }else if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","fixed","variance","correlation","ranef"), several.ok = TRUE)
    effects[effects== "fixed"] <- "mean"
    
    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform
    
    effects2 <- effects
    if("ranef" %in% effects){
        if(type.pattern!="CS"){
            stop("Can only extract random effects for \"CS\" structure. \n")
        }
        if(object$design$vcov$heterogeneous){
            stop("Can only extract random effects for \"CS\" structure with homogeneous structure. \n")
        }
        if(length(effects)>1){
            stop("Argument \'effects\' should be of length 1 when it contains \"ranef\". \n")
        }
    }
    if(transform.rho == "cov"){
        if(all("correlation" %in% effects == FALSE)){
            stop("Cannot use the argument \'transform.rho\' set to \"cov\" when \"correlation\" is not in argument \'effect\'. \n")
        }
        if(all("variance" %in% effects == FALSE)){
            effects2 <- c("variance",effects2)
        }
    }
    if(!is.null(p)){
        if(any(duplicated(names(p)))){
            stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
        }
        if(any(param.name %in% names(p) == FALSE)){
            stop("Incorrect argument \'p\': missing parameter(s) \"",paste(param.name[param.name %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
        }
        p <- p[param.name]
        if(object$reparametrize$transform){
            reparametrize.p <- .reparametrize(p = p[object.reparametrize.name],  
                                              type = param.type[object.reparametrize.name],
                                              sigma = param.sigma[object.reparametrize.name],
                                              k.x = param.k.x[object.reparametrize.name],
                                              k.y = param.k.y[object.reparametrize.name],
                                              level = param.level[object.reparametrize.name],                                              
                                              Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                              transform.sigma = transform.sigma,
                                              transform.k = transform.k,
                                              transform.rho = transform.rho,
                                              transform.names = FALSE)$p
        }else{
            reparametrize.p <- p[object.reparametrize.name]
        }
    }else{
        p <- object$param
        reparametrize.p <- object.reparametrize.value
    }

    ## ** extract
    out <- NULL
    if("mean" %in% effects2){
        out <- c(out, p[param.type=="mu"])
    }

    if("ranef" %in% effects2){
        ## raw residuals
        df.epsilon <- stats::residuals(object, p = p, keep.data = TRUE, type = "response")
        if(length(index.na)>0){
            df.epsilon <- df.epsilon[-index.na,,drop=FALSE]
        }
        ## inverse residual variance-covariance matrix
        U.cluster <- unique(df.epsilon$XXcluster.indexXX)
        ls.OmegaM1 <- stats::sigma(object, p = p, cluster = U.cluster, inverse = TRUE, simplifies = FALSE) 
        ls.epsilon <- base::tapply(df.epsilon$r.response,df.epsilon$XXcluster.indexXX,list)[names(ls.OmegaM1)]## split residuals by id
        cluster.strata <- U.strata[unlist(Upattern$index.strata)[match(attr(ls.OmegaM1,"pattern"), Upattern$name)]]
        ## covariance parameter(s)
        tau <- stats::coef(object, effects = "correlation", transform.rho = "cov", transform.names = FALSE)
        tau.strata <- U.strata[unlist(param.strata[match(names(tau),param.name)])]

        ## design matrix
        df.Z <- as.data.frame(X.cor) ## e.g. school=1,id=1,2,3,....
        if(n.strata==1){
            Z.strata <- stats::setNames(rep(U.strata, NCOL(df.Z)),names(df.Z))
        }else{
            Z.strata <- sapply(attr(X.cor,"ls.level"), function(iRow){iRow[,strata.var]})
        }
        M.ranef <- do.call(rbind,lapply(U.strata, function(iStrata){ ## iStrata <- U.strata[1]
            iCluster <- U.cluster[cluster.strata==iStrata]
            iIndex.cluster <- index.cluster[iCluster]
            iUpattern <- Upattern[U.strata[unlist(Upattern$index.strata)]==iStrata,,drop=FALSE]
            iZ <- df.Z[unlist(iIndex.cluster),names(which(Z.strata==iStrata)),drop=FALSE]
            iNobs <- sapply(iIndex.cluster,length)
            iIndex.cluster2 <- mapply(x = cumsum(iNobs), y = c(1,cumsum(iNobs[-length(iNobs)])+1), function(x,y){y:x}, SIMPLIFY = FALSE)

            iRanef <- .ranef(design = iZ,
                             tau = tau[tau.strata == iStrata],
                             OmegaM1 = ls.OmegaM1[iCluster],
                             epsilon = ls.epsilon[iCluster],
                             cluster.var = cluster.var,
                             index.cluster = iIndex.cluster2,
                             Upattern = iUpattern,
                             Xpattern.cor = Xpattern.cor[unique(iUpattern$cor)])
            rownames(iRanef) <- as.character(U.cluster.original[iCluster])        
            return(iRanef)
        }))
        return(M.ranef[as.character(U.cluster.original),,drop=FALSE])
    }
    if(any(c("variance","correlation") %in% effects2)){
        pVar <- NULL
        if("variance" %in% effects2){
            if(test.notransform){
                index.sigmak <- names(param.type)[param.type %in% c("sigma","k")]
                if(transform.names && !is.null(object.reparametrize.newname)){
                    pVar <- c(pVar, stats::setNames(reparametrize.p[index.sigmak],object.reparametrize.newname[match(index.sigmak,names(reparametrize.p))]))
                }else{
                    pVar <- c(pVar, reparametrize.p[index.sigmak])
                }                    
            }else{
                pVar <- c(pVar, p[param.name[param.type %in% c("sigma","k")]])
            }
        }
        if("correlation" %in% effects2){
            if(test.notransform){
                index.rho <- names(param.type)[param.type %in% c("rho")]
                if(transform.names && !is.null(object.reparametrize.newname)){
                    pVar <- c(pVar, stats::setNames(reparametrize.p[index.rho],object.reparametrize.newname[match(index.rho,names(reparametrize.p))]))
                }else{
                    pVar <- c(pVar, reparametrize.p[index.rho])
                }                    
            }else{
                pVar <- c(pVar, p[param.name[param.type %in% c("rho")]])
            }
        }
        if(!test.notransform){
            ls.reparam <- .reparametrize(p = pVar,  
                                         type = param.type[names(pVar)], 
                                         sigma = param.sigma[names(pVar)], 
                                         k.x = param.k.x[names(pVar)], 
                                         k.y = param.k.y[names(pVar)], 
                                         level = param.level[names(pVar)], 
                                         Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
            outVar <- ls.reparam$p
            if(ls.reparam$transform){
                newname <- stats::setNames(ls.reparam$newname,names(pVar))
            }else{
                newname <- NULL
            }            
        }else{
            outVar <- pVar
            newname <- NULL
        }
        out <- c(out,outVar)

        }else{
            newname <- NULL
        }
    
    ## ** post process
    if(("variance" %in% effects2) && ("variance" %in% effects == FALSE)){
        index.rm <- which(names(newname) %in% param.name[param.type %in% c("sigma","k")])
        newname <- newname[-index.rm]
        out <- out[-index.rm]
    }
    if(length(newname)>0){
        ## rename
        names(out)[match(names(newname),names(out))] <- as.character(newname)
    }
    return(out)
}


## * .ranef
##' @description estimate random effect in a given strata
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
.ranef <- function(design, tau, OmegaM1, epsilon,
                   cluster.var, index.cluster, Upattern, Xpattern.cor){

    n.cluster <- length(index.cluster)

    ## ** build design matrix
    if(NCOL(design)==1 && all(design==1)){ 
        ls.designA <- stats::setNames(list(design), cluster.var)        
    }else{
        ls.designA <- c(stats::setNames(list(rep(1,NROW(design))), cluster.var),
                        lapply(design,as.factor))
    }
    df.designA <- as.data.frame(ls.designA)
    names(df.designA) <- names(ls.designA)

    ls.Z.tempo <- lapply(names(ls.designA), function(iZ){stats::model.matrix(stats::as.formula(paste0("~0+",iZ)),df.designA)})
    X.Z <- do.call(cbind,ls.Z.tempo) ## convert factor to dummy variables
    ls.Z <- lapply(index.cluster, function(iIndex){t(X.Z[iIndex,,drop=FALSE])})

    ## ** get variance of the random effects
    if(length(tau)==1){
        G <- diag(tau, nrow = length(tau)) ## variance-covariance matrix of the random effects
    }else{
        index.cluster2 <- rep(NA, NROW(X.Z))
        pattern.cluster <- rep(NA, n.cluster)

        for(iCluster in 1:length(index.cluster)){
            index.cluster2[index.cluster[[iCluster]]] <- iCluster
        }
        U.indexcluster <- unique(index.cluster2)
        for(iPattern in 1:NROW(Upattern)){
            pattern.cluster[Upattern$index.cluster[[iPattern]]] <- Upattern$name[iPattern]
        }

        ls.corparam <- lapply(colnames(X.Z), FUN = function(iName){ ## iName <- colnames(X.Z)[2]
browser()
                iCluster <- min(index.cluster2[X.Z[,iName]==1]) ## find cluster relevant for the type of random effect
                iPattern <- which(Upattern$name == pattern.cluster[iCluster]) ## covariance pattern corresponding to the cluster
                iPattern.index <- attr(Xpattern.cor[[iPattern]],"index.pairtime") ## 2 x T matrix: all pairs of observations
                index.obs <- which(ls.Z[[iCluster]][iName,]==1) ## identify observations relevant for the type of random effect in the cluster
                index.pair <- which((iPattern.index[1,] %in% index.obs) * (iPattern.index[2,] %in% index.obs)==1) ## identify all pairs affected by the random effect
                param.obs <- Xpattern.cor[[iPattern]][index.pair,,drop=FALSE] ## identify correlation parameter(s) corresponding to the random effect
                return(names(which(colSums(param.obs)!=0)))
                
            })
            if(all(sapply(ls.corparam, length)==1)){
                G <- diag(tau[unlist(ls.corparam)], nrow = length(ls.corparam)) ## variance-covariance matrix of the random effects
            }else{
                ## attempt to fix the case where the is confusion between two levels of nesting
                ## e.g. ID/Day/Scan: At baseline (i.e. Day 1) only a single scan was performed
                ##                   so by default it would be associated to Scan instead of Day
                grp.corparam <- mapply(cumsum(c(0,sapply(ls.Z.tempo[1:(length(ls.Z.tempo)-1)],NCOL)))+1,cumsum(sapply(ls.Z.tempo,NCOL)),
                                       FUN = function(x,y){x:y})
                for(iG in 1:length(grp.corparam)){ ## iG <- 2
                    if(sum(!duplicated(ls.corparam[grp.corparam[[iG]]]))>1){
                        warning("Possible issue when identifying the correlation parameters associated to the random effects. \n")
                        ls.corparam[grp.corparam[[iG]]] <- lapply(1:length(grp.corparam[[iG]]), function(iL){unique(unlist(ls.corparam[grp.corparam[[iG]]]))})
                    }
                }

                nametau.order <- names(sort(table(unlist(ls.corparam)), decreasing = FALSE)) ## identify nesting
                tau.order <- tau[nametau.order] - c(0,tau[nametau.order[-length(nametau.order)]]) ## re-order and reparametrize covariance parameters
                G <- diag(sapply(ls.corparam, function(iRho){ ## xxx <- ls.corparam[[1]]
                    tau.order[sort(factor(iRho, levels = nametau.order))[1]]
                }))
            }
           
        }
browser()
    ## ** compute random effects
    M.ranef <- t(do.call(cbind,lapply(1:n.cluster, function(iC){ ##  iC <- 2 
        G %*% ls.Z[[iC]] %*% OmegaM1[[iC]] %*% epsilon[[iC]]
    })))
    rownames(M.ranef) <- names(OmegaM1)
    colnames(M.ranef) <- colnames(X.Z)

    return(M.ranef)
}

##----------------------------------------------------------------------
### coef.R ends here
