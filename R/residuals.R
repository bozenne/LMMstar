### residuals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:40) 
## Version: 
## Last-Updated: Dec 19 2021 (00:36) 
##           By: Brice Ozenne
##     Update #: 447
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * residuals.lmm (documentation)
##' @title Extract The Residuals From a Linear Mixed Model
##' @description Extract or compute the residuals of a linear mixed model.
##' @name residuals
##' 
##' @param object a \code{lmm} object.
##' @param type [character] type of residual to output: raw residuals (\code{"response"}), Pearson residuals (\code{"pearson"}), normalized residuals (\code{"normalized"}, scaled residual \code{"scaled"}), or partial residuals (\code{"partial"} or \code{"partial-ref"}). Can also be \code{"all"} to output all except partial residuals. See detail section.
##' @param format [character] Should the residuals be output relative as a vector (\code{"long"}), or as a matrix with in row the clusters and in columns the outcomes (\code{"wide"}).
##' @param data [data.frame] dataset relative to which the residuals should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param var [character vector] name of the variable relative to which the partial residuals should be computed.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the residuals. Only relevant if differs from the fitted values.
##' @param plot [character] Should a qqplot (\code{"qqplot"}), or a heatmap of the correlation between residuals  (\code{"correlation"}, require wide format), or a plot of residuals along the fitted values (\code{"scatterplot"}, require long format) be displayed?
##' @param engine.qqplot [character] Should ggplot2 or qqtest be used to display quantile-quantile plots? Only used when argument \code{plot} is \code{"qqplot"}.
##' @param add.smooth [logical] should a local smoother be used to display the mean of the residual values across the fitted values. Only relevant for \code{plot="scatterplot"}.
##' @param digit.cor [integer, >0] Number of digit used to display the correlation coefficients? No correlation coefficient is displayed when set to 0. Only used when argument \code{plot} is \code{"correlation"}.
##' @param size.text [numeric, >0] Size of the font used to displayed text when using ggplot2.
##' @param keep.data [logical] Should the argument \code{data} be output along side the residuals? Only possible in the long format.
##' @param scales [character] Passed to \code{ggplot2::facet_wrap}.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details The argument \code{type} defines how the residuals are computed:
##' \itemize{
##' \item \code{"fitted"}: fitted value \eqn{X_{ij} \hat{\beta}}.
##' \item \code{"raw"}: observed outcome minus fitted value \eqn{\varepsilon = Y_{ij} - X_{ij} \hat{\beta}}.
##' \item \code{"pearson"}: each raw residual is divided by its modeled standard deviation \eqn{\varepsilon = \frac{Y_{ij} - X_{ij} \hat{\beta}}{\sqrt{\hat{\omega}_{ij}}}}.
##' \item \code{"studentized"}: same as \code{"pearson"} but excluding the contribution of the cluster in the modeled standard deviation  \eqn{\varepsilon = \frac{Y_{ij} - X_{ij} \hat{\beta}}{\sqrt{\hat{\omega}_{ij}-\hat{q}_{ij}}}}.
##' \item \code{"normalized"}: raw residuals are multiplied, within clusters, by the inverse of the (lower) Cholesky factor of the modeled residual variance covariance matrix \eqn{\varepsilon = ( Y_{i} - X_{i} \hat{\beta} )\hat{C}^{-1}}.
##' \item \code{"normalized2"}: same as \code{"normalized"} but excluding the contribution of the cluster in the modeled residual variance covariance matrix \eqn{\varepsilon = ( Y_{i} - X_{i} \hat{\beta} ) \hat{D}_i^{-1}}.
##' \item \code{"scaled"}: corresponds to the scaled residuals of PROC MIXED in SAS.
##' \item \code{"partial"} or \code{"partial-ref"}: the partial residual are computed as the raw residuals plus the effect of the covariates in argument \code{var}.
##' \code{"partial"} uses \eqn{\hat{\beta} X  + \hat{\varepsilon}} where \eqn{X} is centered while \code{"partial-ref"} uses \eqn{\hat{\beta} X + \hat{\gamma} Z  + \hat{\varepsilon}} where the Z value are the same for all observations, i.e. uses a reference level.
##' }
##' where
##' \itemize{
##' \item \eqn{X} the design matrix (default) or the design matrix restricted to the variable(s) in argument \code{var} (partial residuals).
##' \item \eqn{Y} the outcome
##' \item \eqn{Z} not defined (default) or the design matrix restricted to the variable(s) not in argument \code{var} (partial residuals).
##' \item \eqn{\hat{\beta}} the estimated mean coefficients relative to \eqn{X}
##' \item \eqn{\hat{\gamma}} the estimated mean coefficients relative to \eqn{Z}
##' \item \eqn{\hat{\Omega}} the modeled variance-covariance of the residuals and \eqn{\hat{\omega}} its diagonal elements
##' \item \eqn{\hat{C}} the lower Cholesky factor of \eqn{\hat{\Omega}}, i.e. \eqn{\hat{C} \hat{C}^{t} = \hat{\Omega}}
##' \item \eqn{\hat{Q}_i= X_i (X^{t}\hat{\Omega}X)^{-1}X_i^{t}} a cluster specific correction factor, approximating the contribution of cluster i to \eqn{\hat{\Omega}}. Its diagonal elements are denoted \eqn{\hat{q}_i}.
##' \item \eqn{\hat{D}_i} the lower Cholesky factor of \eqn{\hat{\Omega}-\hat{Q}_i}
##' }
##'
##' @return
##' When argument format is \code{"long"} and type.oobject is \code{"lmm"}, a vector containing the value of the residual realtive to each observation.
##' It is a matrix if the argument \code{type} contains several values.
##' When argument format is \code{"wide"} and type.oobject is \code{"lmm"}, a data.frame with the value of the residual relative to each cluster (in rows) at each timepoint (in columns).
##' 
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit Linear Model
##' e.lm <- lmm(Y ~ visit + X1 + X2 + X5, data = dL)
##' residuals(e.lm, type = "partial", var = "X1")
##' residuals(e.lm, type = "partial", var = "X1", keep.data = TRUE)
##' 
##' ## fit Linear Mixed Model
##' eUN.lmm <- lmm(Y ~ visit + X1 + X2 + X5,
##'                repetition = ~visit|id, structure = "UN", data = dL)
##'
##' ## residuals
##' residuals(eUN.lmm, format = "long", type = c("normalized","pearson"))
##' residuals(eUN.lmm, format = "wide", plot = "correlation")
##' residuals(eUN.lmm, format = "wide", type = "normalized")
##' residuals(eUN.lmm, format = "wide", type = "scaled")
##'
##' ## residuals and predicted values
##' residuals(eUN.lmm, type = "all")
##' residuals(eUN.lmm, type = "all", keep.data = TRUE)


## * residuals.lmm (code)
##' @rdname residuals
##' @export
residuals.lmm <- function(object, type = "response", format = "long",
                          data = NULL, p = NULL, keep.data = FALSE, var = NULL,
                          plot = "none", engine.qqplot = "ggplot2", add.smooth = TRUE, digit.cor = 2, size.text = 16, scales = "free", ...){

    options <- LMMstar.options()
    type.residual <- type

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    ## check format
    format <- match.arg(format, c("wide","long"))
    if(keep.data && format == "wide"){
        stop("Argument \'keep.data\' must be \"FALSE\" when using the wide format. \n")
    }
    ## check type.residuals
    if(identical("all",tolower(type.residual))){
        type.residual <- c("fitted","response","studentized","pearson","normalized","normalized2","scaled")
    }
    type.residual <- match.arg(type.residual, c("fitted","response","studentized","pearson","normalized","normalized2","scaled","partial","partial-ref"), several.ok = (format=="long"))
    if(any(type.residual %in% c("studentized","pearson","normalized","normalized2","scaled"))){
        effects <- c("mean","variance")
    }else{
        effects <- "mean"
    }    
    name.residual <- paste0("r.",gsub("-ref","",type.residual,fixed = TRUE))
    if("fitted" %in% type.residual){
        name.residual[type.residual=="fitted"] <- "fitted"
    }
    ## special checks for partial residuals
    if("partial" %in% type.residual || "partial-ref" %in% type.residual){
        if((length(type.residual)>2) || (length(type.residual) == 2 && "response" %in% type.residual == FALSE)){
            stop("Argument \'type.residual\' should have length 1 when it contains \"partial\" or  \"partial-ref\". \n",
                 "It can also have length 2 but then the second element should be \"response\". \n")
        }
        if(is.null(var)){
            stop("Argument \'var\' should indicate the covariate effects to preserve when computing the partial residuals. \n")
        }
        name.X <- attr(object$design$mean,"variable")
        keep.intercept <- "(Intercept)" %in% var
        if(keep.intercept && "(Intercept)" %in%  colnames(object$design$mean) == FALSE){
            stop("Argument \'var\' cannot contain \"(Intercept)\" when the model does no include an intercept. \n")
        }
        var <- setdiff(var,"(Intercept)")
        if(any(var %in% name.X == FALSE)){
            stop("Argument \'var\' should refer to covariate(s) of the mean structure. \n",
                 "Valid covariates: \"",paste(name.X, collapse = "\" \""),"\". \n",
                 "Invalid covariates: \"",paste(var[var %in% name.X == FALSE],collapse="\" \""),"\". \n")
        }
    }else{
        if(!is.null(var)){
            warning("Argument \'var\' ignored when computing residuals other than partial residuals. \n")
        }
        keep.intercept <- TRUE
    }
    ## check plot
    plot <- match.arg(plot, c("none","qqplot","correlation","scatterplot","scatterplot2"))
    if(plot == "correlation" && object$time$n == 1){
        stop("Cannot display the residual correlation over time when there is only a single timepoint. \n")
    }
    if(length(add.smooth)==1){
        add.smooth <- rep(add.smooth,2)
    }
    if(length(type.residual)>1 && plot != "none"){
        stop("Argument \'plot\' must be \"none\" when exporting several types of residuals. \n")
    }
    if(plot!="none"){
        label.residual <- switch(type.residual,
                                 "fitted" = "Fitted values",
                                 "response" = "Raw residuals",
                                 "studentized" = "Studentized residuals",
                                 "pearson" = "Pearson residuals",
                                 "normalized" = "Normalized residuals",
                                 "normalized2" = "Pearson Normalized residuals",
                                 "scaled" = "Scaled residuals",
                                 "partial" = paste("Partial residuals for ",paste(var,collapse=", ")),
                                 "partial-ref" = paste("Partial residuals for ",paste(var,collapse=", ")))
    }
    ## check agreement plot, format, type.residual
    if(length(type.residual)>1 && format == "wide"){
        stop("Argument \'format\' must be \"long\" when exporting several types of residuals. \n")
    }
    if(plot == "correlation" && format == "long"){
        stop("Argument \'format\' must be \"wide\" to display the correlation between residuals. \n")
    }
    ## check data
    if(keep.data && !is.null(data) && any(colnames(data) %in% name.residual)){
        stop("Argument \'data\' should not contain a column named \"",paste(name.residual[name.residual %in% colnames(data)], collapse = "\" \""),"\". \n",
             "This name is used to export the residuals. \n")
    }

    ## ** update design
    if(any(c("partial","partial-ref") %in% type.residual)){
        ## extract data and design matrix
        if(is.null(data)){
            data <- stats::model.frame(object)
            data <- data[,setdiff(colnames(data),c("XXindexXX", "XXstrata.indexXX", "XXvisit.indexXX")),drop=FALSE]
            design <- model.matrix(object, effects = effects, simplifies = FALSE)
        }else{
            design <- model.matrix(object, data = data, effects = effects, simplifies = FALSE)
        }

        ## design matrix relative to the data where the effect of variables no in var has been removed
        reference <- NULL
        centering <- NULL
        if("partial-ref" %in% type.residual){ 
            ## set the dataset at the reference value for all variables not in var
            if(is.null(attr(type.residual,"reference"))){
                name.Xfac <- names(object$xfactor$mean)
                name.Xnum <- setdiff(name.X,name.Xfac)
                if(length(name.Xfac)>0){
                    reference <- c(reference, lapply(object$xfactor$mean,function(iX){factor(iX[1], levels = iX)}))
                }
                if(length(name.Xnum)>0){
                    reference <- c(reference, as.list(stats::setNames(rep(0, length(name.Xnum)), name.Xnum)))
                }
                reference <- data.frame(reference, stringsAsFactors = FALSE)
            }else{
                reference <- attr(type.residual,"reference")
            }
            resdata <- data
            if(length(setdiff(name.X,var))>0){
                for(iVar in setdiff(name.X,var)){
                    data[[iVar]] <- reference[[iVar]]
                }
            }

            ## build design matrix
            design2 <- model.matrix(object, data = data, effects = effects, simplifies = FALSE)

            ## handle intercept term
            if(keep.intercept==FALSE && "(Intercept)" %in% colnames(design2$mean)){
                design2$mean[,"(Intercept)"] <- 0
            }
        }else if("partial" %in% type.residual){
            centering <- colMeans(object$design$mean)
            design2 <- design
            design2$mean <- sweep(design$mean, FUN = "-", MARGIN = 2, STATS = centering)
        }
        
    }else{
        design <- model.matrix(object, data = data, effects = effects, simplifies = FALSE)
        if(keep.data){
            data <- stats::model.frame(object)
            data <- data[,setdiff(colnames(data),c("XXindexXX", "XXstrata.indexXX", "XXvisit.indexXX")),drop=FALSE]
        }
    }
    Y <- design$Y
    X <- design$mean
    structure <- design$vcov
    n.cluster <- design$cluster$n
    precompute.XX <- design$precompute.XX
    cluster.level <- design$cluster$levels
    index.cluster <- design$index.cluster
    index.time <- design$index.time

    index.variance <- structure$X$pattern.cluster
    n.pattern <-  NROW(structure$X$Upattern)

    ## ** update Omega
    if(!is.null(p)){
        if(any(duplicated(names(p)))){
            stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
        }
        if(any(names(object$param$type) %in% names(p) == FALSE)){
            stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(object$param$type)[names(object$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
        }
        beta <- p[names(object$param$type=="mu")]
        if(any(type.residual %in% c("studentized","pearson","normalized","normalized2","scaled"))){
            Omega <- .calc_Omega(object = structure, param = p)
            precision <- lapply(Omega, solve)
        }
    }else{
        beta <- object$param$value[object$param$type=="mu"]
        if(any(type.residual %in% c("studentized","pearson","normalized","normalized2","scaled"))){
            Omega <- object$Omega
            precision <- object$OmegaM1
        }
    }

    ## ** pre-compute
    sqrtPrecision <- list()
    if("pearson" %in% type.residual){
        sqrtPrecision$pearson <- lapply(Omega,function(iM){1/sqrt(diag(iM))})
    }
    if("studentized" %in% type.residual || "normalized2" %in% type.residual){
        tX.precision.X <- matrix(0, nrow = NCOL(X), ncol = NCOL(X), dimnames = list(colnames(X),colnames(X)))

            if(!is.null(precompute.XX)){

                for (iPattern in 1:n.pattern) { ## iPattern <- 1
                    iOmega <- precision[[iPattern]]
                    iTime <- NCOL(iOmega)
                    iTime2 <- length(iOmega)

                    iX <- matrix(unlist(precompute.XX$pattern[[iPattern]]), nrow = iTime2, ncol = dim(precompute.XX$pattern[[iPattern]])[3], byrow = FALSE)
                    tX.precision.X <- tX.precision.X + (as.double(iOmega) %*% iX)[as.double(precompute.XX$key)]
                }
            }else{
                for(iId in 1:n.cluster){
                    iIndex <- which(index.cluster==iId)
                    iOrder <- order(index.time[iIndex])
                    iX <- X[iIndex[iOrder],,drop=FALSE]
                    tX.precision.X <- tX.precision.X + t(iX) %*% precision[[index.variance[iId]]] %*% iX
                }
            }
        ## equal if proper ordering
        ## t(X) %*% bdiag(lapply(1:n.cluster,function(x){precision[[1]]})) %*% X - tX.precision.X
        ## equal in large samples or when using observed information
        ## information(object, effects = "mean") - tX.precision.X
        tX.precision.X.M1 <- solve(tX.precision.X)
        ## vcov(object, effects = "mean") - tX.precision.X.M1            
    }
    if("normalized" %in% type.residual){
        ## from SAS documentation
        ## If Var[Y]=V and C'C=V, then C'-1 Y has uniform dispersion and its elements are uncorrelated
        ## B=chol(M) gives B'B = M 
        sqrtPrecision$normalized <- lapply(precision,function(iP){t(chol(iP))})
    }

    ## ** raw residuals
    if(!is.null(data) || !is.null(p) || "partial" %in% type.residual || "partial-ref" %in% type.residual){
        fitted <- X %*% beta
        res <-  as.vector(Y - fitted)
        M.res <- matrix(NA, nrow = NROW(X), ncol = length(type.residual), dimnames = list(NULL, name.residual))
    }else{
        fitted <- stats::fitted(object)
        res <- as.vector(object$residuals)
        M.res <- matrix(NA, nrow = NROW(object$residuals), ncol = length(type.residual), dimnames = list(NULL, name.residual))
    }

    ## ** normalization
    if ("fitted" %in% type.residual) {
        M.res[,"fitted"] <- fitted
    }
    if ("response" %in% type.residual) {
        M.res[,"r.response"] <- res
    }
    if ("partial" %in% type.residual){
        index.var <- which(attr(object$design$mean,"variable") %in% var)
        index.col <- which(attr(object$design$mean,"assign") %in% index.var)
        M.res[,"r.partial"] <- design2$mean[,index.col,drop=FALSE] %*% beta[index.col] + res
        attr(M.res,"centering") <- centering
    }
    if ("partial-ref" %in% type.residual){
        M.res[,"r.partial"] <- design2$mean %*% beta + res
        attr(M.res,"reference") <- reference
    }
    if (any(type.residual %in% c("studentized", "pearson", "normalized", "normalized2", "scaled"))) {
        for(iId in 1:n.cluster){ ## iId <- 7
            iIndex <- which(index.cluster==iId)
            iOrder <- order(index.time[iIndex])
            iResidual <- res[iIndex[iOrder]]
            iN.time <- length(iIndex)
                
            for(iType in type.residual){
                if("pearson" %in% type.residual){
                    resnorm <- iResidual * sqrtPrecision$pearson[[index.variance[iId]]]
                    M.res[iIndex,"r.pearson"] <- resnorm[order(iOrder)]
                }
                if("studentized" %in% type.residual){
                    iX <- X[iIndex[iOrder],,drop=FALSE]
                    iQ <- iX %*% tX.precision.X.M1 %*% t(iX)
                    resnorm <- iResidual / sqrt(diag(Omega[[index.variance[iId]]] - iQ))
                    M.res[iIndex,"r.studentized"] <- resnorm[order(iOrder)]
                }
                if("normalized" %in% type.residual){
                    ## resnorm <- as.double(sqrtPrecision$normalized[[index.variance[iId]]] %*% iResidual)
                    resnorm <- as.double(iResidual %*% sqrtPrecision$normalized[[index.variance[iId]]])
                    M.res[iIndex,"r.normalized"] <- resnorm[order(iOrder)]
                }
                if("normalized2" %in% type.residual){
                    iX <- X[iIndex[iOrder],,drop=FALSE]
                    iQ <- iX %*% tX.precision.X.M1 %*% t(iX)
                    resnorm <- as.double(iResidual %*% t(chol(solve(Omega[[index.variance[iId]]] - iQ))))
                    M.res[iIndex,"r.normalized2"] <- resnorm[order(iOrder)]
                }
                if("scaled" %in% type.residual){
                    M.res[iIndex[iOrder][1],"r.scaled"] <- iResidual[1]/attr(Omega[[index.variance[iId]]],"sd")[1]
                    if(iN.time>1){
                        for(iTime in 2:iN.time){ ## iTime <- 2
                            iVar <- Omega[[index.variance[iId]]][iTime,iTime]
                            iPrecision_kk <- solve(Omega[[index.variance[iId]]][1:(iTime-1),1:(iTime-1),drop=FALSE])
                            iOmega_lk <- Omega[[index.variance[iId]]][iTime,1:(iTime-1),drop=FALSE]
                            iOmega_kl <- Omega[[index.variance[iId]]][1:(iTime-1),iTime,drop=FALSE]
                                
                            num <- iResidual[iTime] - iOmega_lk %*% as.double(iPrecision_kk %*% iResidual[1:(iTime-1)])
                            denom <- iVar - as.double(iOmega_lk %*% iPrecision_kk %*% iOmega_kl)
                            M.res[iIndex[iOrder][iTime],"r.scaled"] <- num/sqrt(denom) ## issue in term of dimension
                        }
                    }
                }
            }
        }
    }

    ## ** add NA
    if(is.null(match.call()$data) && length(object$index.na)>0){
        inflateNA <-  .addNA(index.na = object$index.na, design = design, time = object$time)
        Msave.res <- M.res
        M.res <- matrix(NA, nrow = inflateNA$n.allobs, ncol = length(type.residual), dimnames = list(NULL, name.residual))
        M.res[-object$index.na,] <- Msave.res

        level.cluster <- factor(inflateNA$level.cluster, levels = cluster.level)
        level.time <- factor(inflateNA$level.time, object$time$levels)
    }else{
        level.cluster <- factor(cluster.level[index.cluster], levels = cluster.level)
        level.time <- factor(object$time$levels[index.time], object$time$levels)
    }

    ## plot
    ##
    if(format=="wide"){

        dfL.res <- data.frame(residuals = as.vector(M.res), cluster = level.cluster, time = level.time, stringsAsFactors = FALSE)
        MW.res <- reshape2::dcast(data = dfL.res,
                                  formula = cluster~time, value.var = "residuals")
        attr(MW.res,"reference") <- attr(M.res,"reference")
        attr(MW.res,"centering") <- attr(M.res,"centering")
        if(plot=="qqplot"){
            if(engine.qqplot=="ggplot2"){
                dfL.res$time <- paste0(object$time$var,": ", dfL.res$time)
                attr(MW.res,"plot") <- ggplot2::ggplot(dfL.res, ggplot2::aes(sample = residuals)) + ggplot2::stat_qq() + ggplot2::stat_qq_line() + ggplot2::facet_wrap(~time, scales = scales) + ggplot2::ggtitle(label.residual) + ggplot2::theme(text = ggplot2::element_text(size=size.text))
                print(attr(MW.res,"plot"))
            }else if(engine.qqplot=="qqtest"){
                requireNamespace("qqtest")
                m <- NCOL(MW.res)-1
                sqrtm.round <- ceiling(sqrt(m))
                sqrtm.round2 <- ceiling(m/sqrtm.round)

                    oldpar <- graphics::par(no.readonly = TRUE)   
                    on.exit(graphics::par(oldpar))            
                    graphics::par(mfrow = c(sqrtm.round,sqrtm.round2))

                    lapply(1:m,function(iCol){qqtest::qqtest(stats::na.omit(MW.res[,iCol+1]), main = paste0(object$time$var,": ",colnames(MW.res)[iCol+1]," (",label.residual,")"))})
                }
            }else if(plot %in% c("scatterplot","scatterplot2")){
                dfL.res$time <- paste0(object$time$var,": ", dfL.res$time)
                dfL.res$fitted <- fitted
                attr(MW.res,"plot") <- ggplot2::ggplot(dfL.res) + ggplot2::xlab("Fitted values") + ggplot2::theme(text = ggplot2::element_text(size=size.text))
                if(plot=="scatterplot"){
                    attr(MW.res,"plot") <- attr(MW.res,"plot") + ggplot2::geom_abline(slope=0,intercept=0,color ="red") + ggplot2::geom_point(ggplot2::aes(x = fitted, y = residuals)) + ggplot2::ylab(label.residual)
                    if(add.smooth[1]){
                        attr(MW.res,"plot") <- attr(MW.res,"plot") + ggplot2::geom_smooth(ggplot2::aes(x = fitted, y = residuals), se = add.smooth[2])
                    }
                }else if(plot=="scatterplot2"){
                    label.residual2 <- paste0("|",label.residual,"|")
                    attr(MW.res,"plot") <- attr(MW.res,"plot") + ggplot2::geom_point(ggplot2::aes(x = fitted, y = sqrt(abs(residuals)))) + ggplot2::ylab(bquote(sqrt(.(label.residual2))))
                    if(add.smooth[1]){
                        attr(MW.res,"plot") <- attr(MW.res,"plot") + ggplot2::geom_smooth(ggplot2::aes(x = fitted, y = sqrt(abs(residuals))), se = add.smooth[2])
                    }
                }
                attr(MW.res,"plot") <- attr(MW.res,"plot") + ggplot2::facet_wrap(~time, scales = scales)
                print(attr(MW.res,"plot"))
            }else if(plot == "correlation"){
                name.time <- colnames(MW.res[,-1])
                
                M.cor  <- stats::cor(MW.res[,-1], use = "pairwise")
                ind.cor <- !is.na(M.cor)
                arr.ind.cor <- which(ind.cor, arr.ind = TRUE)
                arr.ind.cor[] <- name.time[arr.ind.cor]

                df.gg <- data.frame(correlation = M.cor[ind.cor], arr.ind.cor,stringsAsFactors = FALSE)
                df.gg$col <- factor(df.gg$col, levels = name.time)
                df.gg$row <- factor(df.gg$row, levels = name.time)
                df.gg$row.index <- match(df.gg$row,name.time)
                df.gg$col.index <- match(df.gg$col,name.time)
                dfR.gg <- df.gg[df.gg$col.index>=df.gg$row.index,,drop=FALSE]
                
                attr(MW.res,"plot") <- ggplot2::ggplot(dfR.gg, ggplot2::aes_string(x = "col", y = "row", fill = "correlation")) + ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Correlation") + ggplot2::xlab(object$time$var) + ggplot2::ylab(object$time$var) + ggplot2::ggtitle(label.residual) + ggplot2::theme(text = ggplot2::element_text(size=size.text))
                if(!is.na(digit.cor) && digit.cor>0){
                    correlation <- NULL ## [[:forCRANcheck:]]
                    attr(MW.res,"plot") <- attr(MW.res,"plot") + ggplot2::geom_text(ggplot2::aes(label = round(correlation,digit.cor)))
                }
                print(attr(MW.res,"plot"))
            }

            if(plot == "none"){
                names(MW.res)[-1] <- paste0(name.residual,".",names(MW.res)[-1])
                return(MW.res)
            }else{
                return(invisible(MW.res))
            }
            
        }else if(format == "long"){

            if(plot == "qqplot"){
                if(engine.qqplot=="ggplot2"){
                    df.gg <- data.frame(residuals = M.res[,1],stringsAsFactors = FALSE)
                    attr(M.res,"plot") <- ggplot2::ggplot(df.gg, ggplot2::aes(sample = residuals)) + ggplot2::stat_qq() + ggplot2::stat_qq_line() + ggplot2::ggtitle(label.residual) + ggplot2::theme(text = ggplot2::element_text(size=size.text))
                    print(attr(M.res,"plot"))
                }else if(engine.qqplot=="qqtest"){
                    requireNamespace("qqtest")
                    qqtest::qqtest(stats::na.omit(M.res[,1]))
                }
            }else if(plot %in% c("scatterplot","scatterplot2")){
                df.gg <- data.frame(fitted = fitted, residuals = M.res[,1],stringsAsFactors = FALSE)

                attr(M.res,"plot") <- ggplot2::ggplot(df.gg) + ggplot2::xlab("Fitted values") + ggplot2::theme(text = ggplot2::element_text(size=size.text))
                if(plot=="scatterplot"){
                    attr(M.res,"plot") <- attr(M.res,"plot") + ggplot2::geom_abline(slope=0,intercept=0,color ="red") + ggplot2::geom_point(ggplot2::aes(x = fitted, y = residuals)) + ggplot2::ylab(label.residual) 
                    if(add.smooth[1]){
                        attr(M.res,"plot") <- attr(M.res,"plot") + ggplot2::geom_smooth(ggplot2::aes(x = fitted, y = residuals), se = add.smooth[2])
                    }
                }else if(plot=="scatterplot2"){
                    label.residual2 <- paste0("|",label.residual,"|")
                    attr(M.res,"plot") <- attr(M.res,"plot") + ggplot2::geom_point(ggplot2::aes(x = fitted, y = sqrt(abs(residuals)))) + ggplot2::ylab(bquote(sqrt(.(label.residual2))))
                    if(add.smooth[1]){
                        attr(M.res,"plot") <- attr(M.res,"plot") + ggplot2::geom_smooth(ggplot2::aes(x = fitted, y = sqrt(abs(residuals))), se = add.smooth[2])
                    }
                }
                print(attr(M.res,"plot"))
            }

            if(keep.data){
                out <- cbind(data,M.res)
                attr(out,"plot") <- attr(M.res,"plot")
                attr(out,"reference") <- attr(M.res,"reference")
                attr(out,"centering") <- attr(M.res,"centering")
            }else if(length(type.residual)==1){
                out <- as.vector(M.res)
                attr(out,"plot") <- attr(M.res,"plot")
                attr(out,"reference") <- attr(M.res,"reference")
                attr(out,"centering") <- attr(M.res,"centering")
            }else{
                out <- M.res
            }

            if(plot == "none"){
                return(out)
            }else{
                return(invisible(out))
            }

        }
    }

## * .addNA
.addNA <- function(index.na, design, time){

    attr.cluster <- attr(index.na,"cluster")
    attr.time <- attr(index.na,"time")

    ## ** if no missing value or missing information about cluster or time returns nothing
    if(length(index.na)==0 || any(is.na(attr.cluster)) || any(is.na(attr.time))){
        return(NULL)
    }

    ## ** identify all clusters
    allcluster <- sort(union(design$cluster$levels, attr.cluster))
    n.allcluster <- length(allcluster)

    n.allobs <- length(design$index.cluster)+length(index.na)
    
    ## ** identify missing clusters
    missing.cluster <- setdiff(attr.cluster, design$cluster$levels)

    ## ** create extended vector of observations
    level.cluster <- rep(NA, n.allobs)
    level.cluster[-index.na] <- design$cluster$levels[design$index.cluster]
    level.cluster[index.na] <- attr.cluster

    level.time <- rep(NA, n.allobs)
    level.time[-index.na] <- time$levels[design$index.time]
    level.time[index.na] <- time$levels[attr.time]

    ## ** export
    return(list(allcluster = allcluster,
                missing.cluster = missing.cluster,
                n.allcluster = length(allcluster),
                n.allobs = length(level.cluster),
                level.cluster = level.cluster,
                level.time = level.time))
    
}
##----------------------------------------------------------------------
### residuals.R ends here
