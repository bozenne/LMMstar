### ranef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 26 2022 (11:18) 
## Version: 
## Last-Updated: maj 10 2024 (19:09) 
##           By: Brice Ozenne
##     Update #: 703
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
##' @param effects [character] should the estimated random effects (\code{"mean"}) or the estimated variance/standard deviation of the random effects (\code{"variance"},\code{"std"}) be output?
##' @param scale [character] should the total variance, variance relative to each random effect, and residual variance be output (\code{"absolute"}).
##' Or the ratio of these variances relative to the total variance (\code{"relative"}).
##' @param se [logical] should standard error and confidence intervals be evaluated using a delta method?
##' Will slow down the execution of the function.
##' @param df [logical] Should degrees of freedom, computed using Satterthwaite approximation, be output.
##' @param p [numeric vector] value of the model coefficients to be used. Only relevant if differs from the fitted values.
##' @param newdata [data.frame] dataset relative to which the random effects should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param format [character] should each type of random effect be output in a data.frame (\code{format="long"})
##' @param transform [logical] should confidence intervals for the variance estimates (resp. relative variance estimates) be evaluated using a log-transform (resp. atanh transformation)?
##' @param simplify [logical] when relevant will convert list with a single element to vectors and omit unessential output.
##' @param ... for internal use.
##'
##' @details Consider the following mixed model:
##' \deqn{Y = X\beta + \epsilon = X\beta + Z\eta + \xi}
##' where the variance of \eqn{\epsilon} is denoted \eqn{\Omega},
##' the variance of \eqn{\eta} is denoted \eqn{\Omega_{\eta}},
##' and the variance of \eqn{\xi} is \eqn{\sigma^2 I} with \eqn{I} is the identity matrix. \cr
##' The random effets are estimating according to:
##' \deqn{E[Y|\eta] = \Omega_{\eta} Z^{t} \Omega^{-1} (Y-X\beta)}
##' 
##' @keywords methods
##' 
##' @return A data.frame or a list depending on the argument \code{format}.
##' 
##' @examples
##' if(require(nlme)){
##' data(gastricbypassL, package = "LMMstar")
##' 
##' ## random intercept
##' e.RI <- lmm(weight ~ time + (1|id), data = gastricbypassL)
##' ranef(e.RI, effects = "mean")
##' ranef(e.RI, effects = "mean", se = TRUE)
##' 
##' ranef(e.RI, effects = "variance")
##' ranef(e.RI, effects = "variance", format = "wide")
##'
##' }

## * ranef.lmm (code)
##' @export
ranef.lmm <- function(object, effects = "mean", scale = "absolute", se = FALSE, df = NULL, transform = (effects %in% c("std","variance")),
                      p = NULL, newdata = NULL, format = "long", simplify = TRUE, ...){

    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    mycall <- match.call()
    if(!inherits(object$design$vcov,"RE")){
        stop("Cannot estimate random effects linear mixed models defined by covariance structure (argument \'structure\'). \n",
             "Consider adding random effects in the argument \'formula\' instead. \n")
    }
    effects <- match.arg(effects, c("mean","std","variance"))
    if(length(scale)==1 && scale == "all"){
        scale <- c("absolute","relative")
    }
    scale <- match.arg(scale, c("absolute","relative"), several.ok = TRUE)
    format <- match.arg(format, c("wide","long"))
    if(!is.null(se)){        
        if(length(se)!=1 || (!is.numeric(se) && !is.logical(se))){
            stop("Argument \'se\' must be a logical value. \n")
        }
        if(se == TRUE && format == "wide"){
            se <- FALSE
            message("Argument \'se\' ignored when using the wide format. \n",
                    "Considering setting argument \'format\' to \"long\". \n")
        }
        if(se>0 && !is.null(p)){
            stop("Cannot evaluate the uncertainty when the argument \'p\' is specified")
        }
    }
    if(is.null(df)){
        df <- se & (effects == "mean")
    }else{
        if(length(df)!=1 || (!is.numeric(df) && !is.logical(df))){
            stop("Argument \'df\' must be a logical value. \n")
        }
        if(se == FALSE && df>0){
            df <- FALSE
            message("Argument \'df\' ignored when the argument \'se\' is FALSE. \n",
                    "Considering setting argument \'se\' to TRUE. \n")
        }
    }
    if(length(transform)!=1 || (!is.numeric(transform) && !is.logical(transform))){
        stop("Argument \'transform\' must be a logical value. \n")
    }
    if(!is.null(mycall$scale) && effects == "mean"){
        message("Argument \'scale\' ignored when argument \'effects\' is set to \"mean\". \n")
    }

    table.param <- object$design$param
    param.name <- table.param$name

    ## p
    if(!is.null(p)){
        init <- .init_transform(p = p, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, 
                                x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho,
                                table.param = object$design$param)
        theta <- init$p
        cumTau <- coef(object, p = theta, effects = c("variance","correlation"), transform.sigma = "square", transform.rho = "cov", transform.names = FALSE)
    }else{
        theta <- object$param
        cumTau <- coef(object, effects = c("variance","correlation"), transform.sigma = "square", transform.rho = "cov", transform.names = FALSE)
    }

    if(transform && se && effects == "mean"){
        transform <- FALSE
        message("Argument \'transform\' is ignore when evaluating confidence intervals for the random effects. \n")
    }

    if(any(object$design$vcov$name$cor[[1]] %in% c("total","relative"))){
        stop("The ranef function cannot be used when the variable(s) w.r.t. which the random are defined are called \"total\" or \"relative\". \n")
    }
    
    ## ** extract from object    
    ## param
    param.type <- stats::setNames(table.param$type,param.name)
    param.rho <- param.name[param.type=="rho"]
    param.sigma <- param.name[param.type=="sigma"]

    ## cluster 
    var.cluster <- object$cluster$var
    n.cluster <- object$design$cluster$n
    index.cluster <- object$design$index.cluster
    Vindex.cluster <- attr(index.cluster, "vectorwise")
    
    ## strata
    var.strata <- object$strata$var
    index.clusterStrata <- object$design$index.clusterStrata
    U.strata <- object$strata$levels
    n.strata <- length(U.strata)

    ## design
    X.cor <- object$design$vcov$cor$X
    Xpattern.cor <- object$design$vcov$cor$Xpattern
    pattern.cluster <- object$design$vcov$pattern
    Upattern <- object$design$vcov$Upattern
    infoRanef <- object$design$vcov$ranef
    
    ## hierarchy
    n.hierarchy <- length(infoRanef$hierarchy)
    df.hierarchy <- data.frame(hierarchy = unlist(lapply(1:n.hierarchy, function(iH){rep(iH, length(infoRanef$hierarchy[[iH]]))})),
                               level = unlist(lapply(1:n.hierarchy, function(iH){1:length(infoRanef$hierarchy[[iH]])})),
                               cluster = unlist(lapply(1:n.hierarchy, function(iH){rep(infoRanef$hierarchy[[iH]][1], length(infoRanef$hierarchy[[iH]]))})),
                               variable = unname(unlist(infoRanef$hierarchy)))
    df.hierarchy$param <- unname(unlist(lapply(infoRanef$param, function(iDF){apply(iDF, MARGIN = 1, identity, simplify = FALSE)}), recursive = FALSE))
    n.RE <- NROW(df.hierarchy)
    
    ## missing values
    index.na <- object$index.na

    ## all vars
    var.all <- unique(lava::manifest(object))

    ## ** converting correlation parameters into random effect variance

    ## *** prepare contrast matrix
    contrastRE <- array(0, dim = c(n.RE+2, length(cumTau), n.strata),
                        dimnames = list(c("total",df.hierarchy$variable,"residual"), names(cumTau), U.strata))

    ## total variance
    tableSigma.param <- table.param[match(param.sigma,table.param$name),,drop=FALSE]
    if(n.strata==1){
        contrastRE["total",param.sigma,U.strata] <- 1
    }else{
        diag(contrastRE["total",tableSigma.param[order(tableSigma.param$index.strata),"name"],U.strata]) <- 1
    }
    ## change in variance within hierarchy
    for(iH in 1:n.hierarchy){ ## iH <- 1
        iHierarchy <- infoRanef$param[[iH]]
        iContrast <- diag(1, nrow = NROW(iHierarchy), ncol = NROW(iHierarchy))
        iContrast[sdiag(iContrast,-1)] <- -1
        for(iS in 1:n.strata){ ## iS <- 1
            contrastRE[rownames(iHierarchy),iHierarchy[,U.strata[iS]],U.strata[iS]] <- iContrast
        }        
    }
    ## residual variance
    if(length(df.hierarchy$variable)==1){
        contrastRE["residual",,] <- contrastRE["total",,] - contrastRE[df.hierarchy$variable,,]
    }else{
        for(iS in 1:n.strata){ ## iS <- 1
            contrastRE["residual",,iS] <- contrastRE["total",,iS] - colSums(contrastRE[df.hierarchy$variable,,iS])
        }
    }

    ## *** variance decomposition
    varDecomp <- apply(contrastRE, MARGIN = 3, FUN = `%*%`, cumTau)
    rownames(varDecomp) <- c("total",df.hierarchy$variable,"residual")
    if(any(varDecomp<=0)){
        stop("Variance for the random effects is found to be negative - cannot estimate the random effects. \n")
    }

    if(effects %in% c("std","variance")){

        ## absolute
        if(se){
            cumTau.var <- vcov(object, p = p, effects = c("variance","correlation"), transform.sigma = "square", transform.rho = "cov", transform.names = FALSE)
            varDecomp.vcov <- apply(contrastRE, MARGIN = 3, FUN = function(iM){iM %*% cumTau.var %*% t(iM)}, simplify = FALSE)
            varDecomp.var <- do.call(cbind,lapply(varDecomp.vcov,diag))
        }else{
            varDecomp.var <- stats::setNames(vector(mode =  "numeric", length = n.strata), U.strata)
        }

        ## relative
        RvarDecomp <- sweep(varDecomp, FUN = "/", MARGIN = 2, STATS = varDecomp["total",])
        if(se){
            RvarDecomp.var <- do.call(cbind,lapply(1:n.strata, function(iS){ ## iS <- 1
                A <- varDecomp[-1,iS]
                B <- varDecomp["total",iS]
                c(total = 0,diag(varDecomp.vcov[[iS]])[-1]/B^2 + varDecomp.vcov[[iS]]["total","total"]*A^2/B^4 - 2*varDecomp.vcov[[iS]][-1,"total"]*A/B^3)
            }))
            colnames(RvarDecomp.var) <- U.strata
        }else{
            RvarDecomp.var <- stats::setNames(vector(mode =  "numeric", length = n.strata), U.strata)
        }

        ## take square root
        if(effects=="std"){
            varDecomp <- sqrt(varDecomp)
            RvarDecomp <- sqrt(RvarDecomp)            

            if(se){
                varDecomp.var <- varDecomp.var/(2*varDecomp)^2
                RvarDecomp.var <- RvarDecomp.var/(2*RvarDecomp)^2
            }
        }

        ## gather results into a single data.frame
        out <- do.call(rbind,mapply(u = as.data.frame(varDecomp), v = as.data.frame(varDecomp.var),
                                    x = as.data.frame(RvarDecomp), y = as.data.frame(RvarDecomp.var), s = colnames(varDecomp),
                                    FUN = function(u,v,x,y,s){ ## x <-
                                        iDF <- data.frame(strata = s,
                                                          type = rep(c("total",df.hierarchy$variable,"residual"),2),
                                                          scale = c(rep("absolute",n.RE+2),rep("relative",n.RE+2)),
                                                          estimate = c(u,x))
                                        if(se){iDF$se <- sqrt(c(v,y))}
                                        return(iDF)
                                    }, SIMPLIFY = FALSE))
        if(se & transform){
            out$lower = exp(log(out$estimate) + stats::qnorm(0.025)*out$se/out$estimate)
            out$upper = exp(log(out$estimate) + stats::qnorm(0.975)*out$se/out$estimate)
        }else if(se){
            out$lower <- out$estimate + stats::qnorm(0.025)*out$se
            out$upper <- out$estimate + stats::qnorm(0.975)*out$se
        }
 
        if(format=="long"){
            out <- out[out$scale %in% scale,,drop=FALSE]
            rownames(out) <- NULL
            if(simplify){
                if(length(scale)==1){
                    out$scale <- NULL
                }
                if(n.strata==1){
                    out$strata <- NULL
                }
                if(se == FALSE && n.strata==1 && length(scale)==1){
                    out <- stats::setNames(out$estimate,out$type)
                }
            }
        }else if(format=="wide"){
            if(!simplify || n.strata>1){
                out$scale <- paste(out$scale,out$strata,sep=".")
            }
            out$strata <- NULL
            out <- stats::reshape(out, direction  = "wide",
                                  idvar = "type", timevar = "scale",
                                  varying = unique(out$scale))
            rownames(out) <- NULL
        
        }
        return(out)

    }else if(effects == "mean"){

        ## extract normalized residuals
        df.epsilon <- stats::residuals(object, newdata = newdata, p = p, keep.data = TRUE, type = "normalized2", format = "long",
                                       simplify = (se==FALSE), fitted.ci = se) ## obtain gradient
        if (n.strata == 1) {
            df.epsilon$XXstrata.indexXX <- U.strata
        }
        ## derivatives
        if(se>0){
            ## normalized residuals 
            grad.epsilon <- attr(df.epsilon,"grad")[,,"r.normalized2"]

            ##  random effect variance (with respect to original parametrisation instead of the cov parametrisation)
            table.param.vcov <- table.param[match(c(param.sigma,param.rho),table.param$name),,drop=FALSE]
            jacobian.none2cov <- .reparametrize(p = theta[table.param.vcov$name], type = table.param.vcov$type, level = table.param.vcov$level, 
                                                sigma = table.param.vcov$sigma, k.x = table.param.vcov$k.x, k.y = table.param.vcov$k.y,
                                                Jacobian = TRUE, dJacobian = FALSE, inverse = FALSE, 
                                                transform.sigma = "none",
                                                transform.k = "none",
                                                transform.rho = "cov",
                                                transform.names = FALSE)$Jacobian
            jacobian.none2default <- .reparametrize(p = theta[table.param.vcov$name], type = table.param.vcov$type, level = table.param.vcov$level, 
                                                  sigma = table.param.vcov$sigma, k.x = table.param.vcov$k.x, k.y = table.param.vcov$k.y,
                                                  Jacobian = TRUE, dJacobian = FALSE, inverse = FALSE, 
                                                  transform.sigma = object$reparametrize$transform.sigma,
                                                  transform.k = object$reparametrize$transform.k,
                                                  transform.rho = object$reparametrize$transform.rho,
                                                  transform.names = FALSE)$Jacobian
            jacobian.cov2default <- solve(jacobian.none2cov) %*% jacobian.none2default
            ls.contrastRE2 <-  apply(contrastRE, MARGIN = 3, function(iC){iC %*% jacobian.cov2default}, simplify = FALSE)
            contrastRE2 <- array(unlist(ls.contrastRE2), dim = dim(contrastRE), dimnames = dimnames(contrastRE))
        }

        ## evaluate random effects
        ls.out <- vector(mode = "list", length = n.RE)
        keep.datacol <- stats::na.omit(unique(c(attr(var.strata,"original"), df.hierarchy$variable)))
        if(se>0){
            options <- LMMstar.options()
            vcov.theta <- vcov(object, effects = "all", df = 2*df, transform.names = FALSE)                
        }

        for(iRE in 1:n.RE){ ## iRE <- 1
            iHierarchy <- df.hierarchy[iRE,"hierarchy"]
            iLevel <- df.hierarchy[iRE,"level"]
            iVar <- df.hierarchy[iRE,"variable"]
            iHvar <- df.hierarchy[df.hierarchy$hierarchy == iHierarchy & df.hierarchy$level <= iLevel,"variable"]

            iSplit <- nlme::collapse(as.data.frame(df.epsilon[iHvar]), as.factor = FALSE)
            if(length(unique(iSplit))==1){
                iSplit.Mdummy <- model.matrix(~1, data.frame(split = iSplit))
            }else{
                iSplit.Mdummy <- model.matrix(~0+split, data.frame(split = iSplit))
            }
            iData <- as.data.frame(matrix(NA, nrow = NCOL(iSplit.Mdummy), ncol = length(keep.datacol),
                                          dimnames = list(NULL,keep.datacol)))
            
            iData[iHvar] <- df.epsilon[apply(iSplit.Mdummy>0, 2, function(iVec){which(iVec)[1]}),iHvar,drop=FALSE]
            iIndex.NNA <- which(!is.na(df.epsilon$r.normalized2))
            if(length(iIndex.NNA)==0){
                iEpsilon.normalized2.sum <- as.double(rbind(df.epsilon$r.normalized2) %*% iSplit.Mdummy)
            }else{
                iEpsilon.normalized2.sum <- as.double(rbind(df.epsilon$r.normalized2[iIndex.NNA]) %*% iSplit.Mdummy[iIndex.NNA,,drop=FALSE])
                iEpsilon.normalized2.sum[colSums(iSplit.Mdummy[iIndex.NNA,,drop=FALSE])==0] <- NA
            }

            if(n.strata==1){
                iStrata <- df.epsilon[[var.strata]][1]
                iTau.strata <- varDecomp[iVar,]
            }else{
                iStrata <- tapply(df.epsilon[[var.strata]], INDEX = iSplit, FUN = unique, simplify = TRUE)
                iStrata.Mdummy <- model.matrix(~0+strata, data.frame(strata = as.factor(iStrata)))
                iTau.strata <- as.double(iStrata.Mdummy %*% cbind(varDecomp[iVar,]))
                iData[attr(var.strata,"original")] <- iStrata ## attr(U.strata,"original")[iStrata,,drop=FALSE]
            }
            ls.out[[iRE]] <- data.frame(hierarchy = iHierarchy,
                                        level = iLevel,
                                        strata = iStrata,
                                        iData,
                                        estimate = iTau.strata * iEpsilon.normalized2.sum)
            
            if(se>0){
                if(length(iIndex.NNA)==0){
                    iGrad <- (t(iSplit.Mdummy) %*% grad.epsilon) * iTau.strata
                }else{
                    iGrad <- (t(iSplit.Mdummy[iIndex.NNA,,drop=FALSE]) %*% grad.epsilon[iIndex.NNA,,drop=FALSE]) * iTau.strata
                    iGrad[colSums(iSplit.Mdummy[iIndex.NNA,,drop=FALSE])==0,] <- NA
                }
                if(n.strata==1){
                    iGrad[,dimnames(contrastRE2)[[2]]] <- iGrad[,dimnames(contrastRE2)[[2]]] + tcrossprod(iEpsilon.normalized2.sum, contrastRE2[iVar,,])
                }else{
                    iGrad[,dimnames(contrastRE2)[[2]]] <- iGrad[,dimnames(contrastRE2)[[2]]] + (iStrata.Mdummy %*% t(contrastRE2[iVar,,])) * iEpsilon.normalized2.sum
                }
                attr(ls.out[[iRE]], "grad") <- iGrad

                ls.out[[iRE]]$se <- sqrt(rowSums(iGrad %*% vcov.theta * iGrad))
                if(df){
                    ls.out[[iRE]]$df <- pmax(.dfX(X.beta = iGrad, vcov.param = vcov.theta, dVcov.param = attr(vcov.theta,"dVcov")), options$min.df)
                    ls.out[[iRE]]$df[is.na(ls.out[[iRE]]$estimate) | is.na(ls.out[[iRE]]$se)] <- NA        
                }else{
                    ls.out[[iRE]]$df <- Inf
                }
                ls.out[[iRE]]$lower <- ls.out[[iRE]]$estimate + stats::qt(0.025, df = ls.out[[iRE]]$df) * ls.out[[iRE]]$se
                ls.out[[iRE]]$upper <- ls.out[[iRE]]$estimate + stats::qt(0.975, df = ls.out[[iRE]]$df) * ls.out[[iRE]]$se
            }
        }
    }

    ## ** export
    out <- do.call(rbind,ls.out)
    out <- out[do.call(order, out[df.hierarchy$variable]),,drop=FALSE]
    rownames(out) <- NULL
    if(se>0){
        if(simplify<=0){
            attr(out,"grad") <- stats::setNames(lapply(ls.out,attr,"grad"), paste(df.hierarchy$hierarchy,df.hierarchy$level, sep = "."))
        }else{
            attr(out,"grad") <- NULL  
        }
    }

    if(format == "wide"){
        out <- stats::setNames(lapply(1:n.hierarchy, function(iH){ ## iH <- 1
            iVar.cluster <- df.hierarchy[df.hierarchy$hierarchy==iH,"cluster"][1]
            iVar.sub <- setdiff(df.hierarchy[df.hierarchy$hierarchy==iH,"variable"],iVar.cluster)

            iOut <- out[out$hierarchy==iH,,drop=FALSE]
            iOut.order <- iOut[order(iOut$level),,drop=FALSE]
            iOut.order$col <- ifelse(iOut.order$level==1,"",paste0(".",nlme::collapse(iOut.order[iVar.sub],sep=":")))
            
            iOutW <- stats::reshape(iOut.order[c(iVar.cluster,"col","estimate")], direction = "wide", 
                                    idvar = iVar.cluster, timevar = "col", sep="")
            rownames(iOutW) <- NULL
            return(iOutW)
        }), sapply(infoRanef$hierarchy,"[",1))
        if(simplify && n.hierarchy == 1){
            names(out[[1]])[1] <- names(out)
            out <- out[[1]]            
        }
        
    }else if(format == "long"){
        if(simplify){
            if(n.hierarchy==1 && all(lengths(infoRanef$hierarchy)==1) && se == FALSE){
                out <- out$estimate
            }else{
                out$strata <- NULL
                if(n.hierarchy==1){
                    out$hierarchy <- NULL
                }
                if(all(lengths(infoRanef$hierarchy)==1)){
                    out$level <- NULL                
                }
            }
        }
    }
    return(out)
}

## * ranef.mlmm (code)
##' @export
ranef.mlmm <- function(object, ...){

    return(lapply(object$model, ranef, ...))
}

##----------------------------------------------------------------------
### ranef.R ends here
