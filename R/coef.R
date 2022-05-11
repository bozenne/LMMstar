### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: maj  9 2022 (19:11) 
##           By: Brice Ozenne
##     Update #: 354
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
##' @param strata [character vector] When not \code{NULL}, only output coefficient relative to specific levels of the variable used to stratify the mean and covariance structure.
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
coef.lmm <- function(object, effects = NULL, strata = NULL, p = NULL,
                     transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE, ...){

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
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }
    
    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform

    effects2 <- effects
    if("ranef" %in% effects){
        if(object$design$vcov$type!="CS"){
            stop("Can only extract random effects for \"CS\" structure. \n")
        }
        if(length(all.vars(object$design$vcov$formula$var))>0){
            stop("Cannot extract random effects when the variance is not constant. \n")
        }
        if(length(effects)>1){
            stop("Argument \'effects\' should be of length 1 when it contains \"ranef\". \n")
        }
    }
    ## if(all("variance" %in% effects == FALSE) && transform.sigma != "none"){
    ##     stop("Cannot use the argument \'transform.sigma\' when \"variance\" is not in argument \'effect\'. \n")
    ## }
    ## if(all("variance" %in% effects == FALSE) && transform.k != "none"){
    ##     stop("Cannot use the argument \'transform.k\' when \"variance\" is not in argument \'effect\'. \n")
    ## }
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
        if(any(names(object$param$type) %in% names(p) == FALSE)){
            stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(object$param$type)[names(object$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
        }
        p <- p[names(object$param$type)]
        if(object$reparametrize$transform){
            reparametrize.p <- .reparametrize(p = p[names(object$reparametrize$p)],  
                                              type = object$param$type[names(object$reparametrize$p)], strata = object$param$strata[names(object$reparametrize$p)],
                                              time.k = object$design$param$time.k, time.rho = object$design$param$time.rho,
                                              name2sd = stats::setNames(object$design$vcov$param$name2,object$design$vcov$param$name),
                                              Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                              transform.sigma = object$reparametrize$transform.sigma,
                                              transform.k = object$reparametrize$transform.k,
                                              transform.rho = object$reparametrize$transform.rho,
                                              transform.names = FALSE)$p
        }else{
            reparametrize.p <- p[names(object$reparametrize$p)]
        }
    }else{
        p <- object$param$value
        reparametrize.p <- object$reparametrize$p
    }

    ## ** extract
    out <- NULL
    if("mean" %in% effects2){
        out <- c(out, p[object$param$type=="mu"])
    }

    if("ranef" %in% effects2){
        ## raw residuals
        df.epsilon <- stats::residuals(object, p = p, keep.data = TRUE, type = "response")
        if(length(object$index.na)>0){
            df.epsilon <- df.epsilon[-object$index.na,,drop=FALSE]
        }
        ## inverse residual variance-covariance matrix
        ls.OmegaM1 <- stats::sigma(object, p = p, cluster = unique(df.epsilon$XXclusterXX), inverse = TRUE) 
        ls.epsilon <- base::tapply(df.epsilon$r.response,df.epsilon$XXclusterXX,list)[names(ls.OmegaM1)]## split residuals by id
        ## covariance parameter(s)
        tau <- stats::coef(object, effects = "correlation", transform.rho = "cov", transform.names = FALSE)
        ## design matrix
        df.Z <- as.data.frame(object$design$vcov$X$cor) ## e.g. school=1,id=1,2,3,....
        df.Z[-1] <- lapply(df.Z[-1],as.factor)
        names(df.Z)[1] <- object$cluster$var
        ls.Z.tempo <- lapply(names(df.Z), function(iZ){stats::model.matrix(stats::as.formula(paste0("~0+",iZ)),df.Z)})
        X.Z <- do.call(cbind,ls.Z.tempo) ## convert factor to dummy variables
        ls.Z <- base::by(X.Z,object$design$index.cluster, function(iDF){t(iDF)}, simplify = FALSE)

        if(length(tau)==1){
            G <- diag(tau, nrow = length(tau)) ## variance-covariance matrix of the random effects
        }else{
            ls.corparam <- lapply(colnames(X.Z), FUN = function(iName){ ## iName <- colnames(X.Z)[2]

                iCluster <- min(object$design$index.cluster[X.Z[,iName]==1]) ## find cluster relevant for the type of random effect
                iPattern <- which(object$design$vcov$X$Upattern$name == object$design$vcov$X$pattern.cluster[iCluster]) ## covariance pattern corresponding to the cluster
                iPattern.index <- attr(object$design$vcov$X$cor.pairwise[[iPattern]],"index.pairtime") ## 2 x T matrix: all pairs of observations
                index.obs <- which(ls.Z[[iCluster]][iName,]==1) ## identify observations relevant for the type of random effect in the cluster
                index.pair <- which((iPattern.index[1,] %in% index.obs) * (iPattern.index[2,] %in% index.obs)==1) ## identify all pairs affected by the random effect
                param.obs <- object$design$vcov$X$cor.pairwise[[iPattern]][index.pair,,drop=FALSE] ## identify correlation parameter(s) corresponding to the random effect
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
        M.ranef <- t(do.call(cbind,lapply(1:object$cluster$n, function(iC){ ##  iC <- 5
            G %*% ls.Z[[iC]][,,drop=FALSE] %*% ls.OmegaM1[[iC]] %*% ls.epsilon[[iC]]
        })))
        rownames(M.ranef) <- names(ls.OmegaM1)
        colnames(M.ranef) <- colnames(X.Z)
        return(M.ranef)
    }

    if(any(c("variance","correlation") %in% effects2)){
        pVar <- NULL
        if("variance" %in% effects2){
            if(test.notransform){
                index.sigmak <- names(object$param$type)[object$param$type %in% c("sigma","k")]
                if(transform.names && !is.null(object$reparametrize$newname)){
                    pVar <- c(pVar, stats::setNames(reparametrize.p[index.sigmak],object$reparametrize$newname[match(index.sigmak,names(reparametrize.p))]))
                }else{
                    pVar <- c(pVar, reparametrize.p[index.sigmak])
                }                    
            }else{
                pVar <- c(pVar, p[object$param$type %in% c("sigma","k")])
            }
        }
        if("correlation" %in% effects2){
            if(test.notransform){
                index.rho <- names(object$param$type)[object$param$type %in% c("rho")]
                if(transform.names && !is.null(object$reparametrize$newname)){
                    pVar <- c(pVar, stats::setNames(reparametrize.p[index.rho],object$reparametrize$newname[match(index.rho,names(reparametrize.p))]))
                }else{
                    pVar <- c(pVar, reparametrize.p[index.rho])
                }                    
            }else{
                pVar <- c(pVar, p[object$param$type %in% c("rho")])
            }
        }
        if(!test.notransform){
            ls.reparam <- .reparametrize(p = pVar,  
                                         type = object$param$type[names(pVar)], strata = object$param$strata[names(pVar)],
                                         time.k = object$design$param$time.k, time.rho = object$design$param$time.rho,
                                         name2sd = stats::setNames(object$design$vcov$param$name2,object$design$vcov$param$name),
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
    if(!identical(effects,effects2)){
        ## remove variance parameters
        newname <- newname[setdiff(names(newname),names(p[object$param$type %in% c("sigma","k")]))]
        out <- out[setdiff(names(out),names(p[object$param$type %in% c("sigma","k")]))]
    }else if("variance" %in% effects && transform.k %in% c("sd","var","logsd","logvar") && object$strata$n>1 && transform.names){
        ## re-order values when converting to sd with strata (avoid sd0:0 sd0:1 sd1:0 sd1:1 sd2:0 sd2:1 ...)
        ## out.strata <- object$param$strata[names(pVar)]
        ## out.type <- object$param$type[names(pVar)]
        ## index.sd <- which(out.type %in% c("sigma","k"))
        ## savenames <- names(out)
        ## savenames[savenames %in% names(out.type)[index.sd]] <- names(out.type)[index.sd][order(out.strata[index.sd])]
        ## out <- out[savenames]        
    }
    if(!is.null(strata)){
        ## only keep parameters relative to certain strata
        out <- out[object$param$strata[names(out)] %in% strata]
    }
    if(length(newname)>0){
        ## rename
        names(out)[match(names(newname),names(out))] <- as.character(newname)
    }
    return(out)
}
##----------------------------------------------------------------------
### coef.R ends here
