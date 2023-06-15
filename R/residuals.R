### residuals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:40) 
## Version: 
## Last-Updated: jun 15 2023 (17:08) 
##           By: Brice Ozenne
##     Update #: 755
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
##' @param type [character] type of residual to output: raw residuals (\code{"response"}), Pearson residuals (\code{"pearson"}), normalized residuals (\code{"normalized"}, scaled residual \code{"scaled"}), or partial residuals (\code{"partial"} or \code{"partial-center"}). Can also be \code{"all"} to output all except partial residuals. See detail section.
##' @param var [character vector] name of the variable relative to which the partial residuals should be computed.
##' @param data [data.frame] dataset relative to which the residuals should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the residuals. Only relevant if differs from the fitted values.
##' @param format [character] Should the residuals be output relative as a vector (\code{"long"}), or as a matrix with in row the clusters and in columns the outcomes (\code{"wide"}).
##' @param keep.data [logical] Should the argument \code{data} be output along side the residuals? Only possible in the long format.
##' @param keep.by [logical or "only"] How should the variable used to index the split of the dataset/LMM be output:
##' only the original variable(s) (FALSE), only the by variable (\code{"only"}), or both (TRUE) .
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details The argument \code{type} defines how the residuals are computed:
##' \itemize{
##' \item \code{"fitted"}: fitted value \eqn{X_{ij} \hat{\beta}}.
##' \item \code{"response"}: raw residual, i.e. observed outcome minus fitted value \eqn{\varepsilon_{ij} = Y_{ij} - X_{ij} \hat{\beta}}.
##' \item \code{"pearson"}: each raw residual is divided by its modeled standard deviation \eqn{\varepsilon_{ij} = \frac{Y_{ij} - X_{ij} \hat{\beta}}{\sqrt{\hat{\omega}_{ij}}}}.
##' \item \code{"studentized"}: same as \code{"pearson"} but excluding the contribution of the cluster in the modeled standard deviation  \eqn{\varepsilon_{ij} = \frac{Y_{ij} - X_{ij} \hat{\beta}}{\sqrt{\hat{\omega}_{ij}-\hat{q}_{ij}}}}.
##' \item \code{"normalized"}: raw residuals are multiplied, within clusters, by the inverse of the (upper) Cholesky factor of the modeled residual variance covariance matrix \eqn{\varepsilon_{ij} = ( Y_{i} - X_{i} \hat{\beta} )\hat{C}^{-1}}.
##' \item \code{"normalized2"}: same as \code{"normalized"} but excluding the contribution of the cluster in the modeled residual variance covariance matrix \eqn{\varepsilon_{ij} = ( Y_{i} - X_{i} \hat{\beta} ) \hat{D}_i^{-1}}.
##' \item \code{"scaled"}: scaled residuals (see PROC MIXED in SAS). Numerically identical to \code{"normalized"} but computed by sequentially scaling and centering the residuals, to make them conditionally independent of previous residuals from the same cluster at previous repetitions.
##' \item \code{"partial"}: partial residuals (\eqn{\gamma E + \hat{\varepsilon}}). A reference level can be also be specified via the attribute \code{"reference"} to change the absolute level of the partial residuals.
##' \code{"partial-center"}: partial residuals with centered covariates (\eqn{\gamma E + \hat{\varepsilon}} where \eqn{E} has been centered, i.e., has 0-mean)
##' }
##' where
##' \itemize{
##' \item \eqn{X=(E,W)} the design matrix. For partial residuals, it is split according to the variable(s) in argument \code{var} (\eqn{E}) and the rest (\eqn{W}).
##' \item \eqn{Y} the outcome
##' \item \eqn{\hat{\beta}=(\hat{\gamma},\hat{\delta})} the estimated mean coefficients relative to \eqn{X=(E,W)}
##' \item \eqn{\hat{\Omega}} the modeled variance-covariance of the residuals and \eqn{\hat{\omega}} its diagonal elements
##' \item \eqn{\hat{C}} the upper Cholesky factor of \eqn{\hat{\Omega}}, i.e. upper triangular matrix satisfying \eqn{\hat{C}^{t} \hat{C} = \hat{\Omega}}
##' \item \eqn{\hat{Q}_i= X_i (X^{t}\hat{\Omega}X)^{-1}X_i^{t}} a cluster specific correction factor, approximating the contribution of cluster i to \eqn{\hat{\Omega}}. Its diagonal elements are denoted \eqn{\hat{q}_i}.
##' \item \eqn{\hat{D}_i} the upper Cholesky factor of \eqn{\hat{\Omega}-\hat{Q}_i}
##' }
##'
##' @return
##' When argument format is \code{"long"} and type.oobject is \code{"lmm"}, a vector containing the value of the residual realtive to each observation.
##' It is a matrix if the argument \code{type} contains several values.
##' When argument format is \code{"wide"} and type.oobject is \code{"lmm"}, a data.frame with the value of the residual relative to each cluster (in rows) at each timepoint (in columns).
##' 
##' @keywords methods
##' 
##' @examples
##' 
##' #### simulate data in the long format ####
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' #### Linear Model ####
##' e.lm <- lmm(Y ~ visit + X1 + X2 + X6, data = dL)
##'
##' ## partial residuals
##' residuals(e.lm, type = "partial", var = "X6")
##' ## residuals(e.lm) + dL$X6 * coef(e.lm)["X6"]
##' e.reslm <- residuals(e.lm, type = "partial", var = "X6", keep.data = TRUE)
##' plot(e.reslm)
##'
##' ## partial residuals with specific reference
##' type <- "partial"
##' attr(type,"reference") <- data.frame(visit=factor(2,1:3),X2=0,X6=3)
##' residuals(e.lm, type = type, var = "X1")
##' ## residuals(e.lm) + dL$X1 * coef(e.lm)["X1"] + coef(e.lm)["visit2"]
##' 
##' ## partial residuals with centered covariates
##' residuals(e.lm, type = "partial-center", var = "X1")
##' ## residuals(e.lm) + (dL$X1-mean(dL$X1)) * coef(e.lm)["X1"]
##'
##' #### Linear Mixed Model ####
##' eUN.lmm <- lmm(Y ~ visit + X1 + X2 + X5 + X6,
##'                repetition = ~visit|id, structure = "UN", data = dL)
##'
##' ## residuals
##' e.reslmm <- residuals(eUN.lmm, type = "normalized", keep.data = TRUE)
##' plot(e.reslmm)
##' plot(e.reslmm, type = "correlation")
##' plot(e.reslmm, type = "scatterplot", labeller = ggplot2::label_both)
##' residuals(eUN.lmm, format = "wide", type = "normalized")
##'
##' ## residuals and predicted values
##' residuals(eUN.lmm, type = "all")
##' residuals(eUN.lmm, type = "all", keep.data = TRUE)
##' 
##' ## partial residuals
##' residuals(eUN.lmm, type = "partial", var = c("(Intercept)","X6"))
##' residuals(eUN.lmm, type = "partial", var = c("X6"))

## * residuals.lmm (code)
##' @export
##' @rdname residuals
residuals.lmm <- function(object, type = "response", var = NULL, 
                          data = NULL, p = NULL, format = "long", keep.data = FALSE, ...){

    options <- LMMstar.options()
    type.residual <- type
    
    ## ** extract from object
    xfactorMu <- object$xfactor$mean
    variableMu.name <- attr(object$design$mean,"variable")
    variableMuFac.name <- names(xfactorMu)
    variableMuNum.name <- setdiff(variableMu.name,variableMuFac.name)

    param.name <- object$design$param$name
    param.type <- stats::setNames(object$design$param$type,param.name)
    param.value <- object$param

    U.time <- object$time$levels
    name.time <- object$time$var
    n.time <- object$time$n
    index.na <- object$index.na
    X.mean <- object$design$mean
    object.residuals <- object$residuals
    object.Omega <- object$Omega
    object.OmegaM1 <- object$OmegaM1

    cluster.levels.original <- object$design$cluster$levels.original
    ## time.levels.original <- object$design$time$levels.original
    
    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    ## check format
    format <- match.arg(format, c("wide","long"))
    if(format=="wide" && identical("all",tolower(type.residual))){
        message("Move to wide format to output all types of residuals. \n")
        format <- "long"
    }
    if(keep.data && format == "wide"){
        stop("Argument \'keep.data\' must be \"FALSE\" when using the wide format. \n")
    }
    ## check type.residuals
    if(identical("all",tolower(type.residual))){        
        type.residual <- c("response","studentized","pearson","normalized","normalized2","scaled")
    }
    attr.ref <- attr(type.residual,"reference")
    type.residual <- match.arg(type.residual, c("response","studentized","pearson","normalized","normalized2","scaled","partial","partial-center"), several.ok = (format=="long"))
    if(any(type.residual %in% c("studentized","pearson","normalized","normalized2","scaled"))){
        effects <- c("mean","variance")
    }else{
        effects <- "mean"
    }    
    name.residual <- paste0("r.",gsub("-center","",type.residual,fixed = TRUE))
    ## special checks for partial residuals
    if("partial" %in% type.residual || "partial-center" %in% type.residual){
        if((length(type.residual)>2) || (length(type.residual) == 2 && "response" %in% type.residual == FALSE)){
            stop("Argument \'type.residual\' should have length 1 when it contains \"partial\" or  \"partial-center\". \n",
                 "It can also have length 2 but then the second element should be \"response\". \n")
        }
        if(is.null(var)){
            stop("Argument \'var\' should indicate the covariate effects to preserve when computing the partial residuals. \n")
        }
        keep.intercept <- "(Intercept)" %in% var
        if(keep.intercept && "(Intercept)" %in%  param.name == FALSE){
            stop("Argument \'var\' cannot contain \"(Intercept)\" when the model does no include an intercept. \n")
        }
        var <- setdiff(var,"(Intercept)")
        if(any(var %in% variableMu.name == FALSE)){
            stop("Argument \'var\' should refer to covariate(s) of the mean structure. \n",
                 "Valid covariates: \"",paste(variableMu.name, collapse = "\" \""),"\". \n",
                 "Invalid covariates: \"",paste(var[var %in% variableMu.name == FALSE],collapse="\" \""),"\". \n")
        }
    }else{
        if(!is.null(var)){
            warning("Argument \'var\' ignored when computing residuals other than partial residuals. \n")
        }
        keep.intercept <- TRUE
    }

    ## check agreement plot, format, type.residual
    if(length(type.residual)>1 && format == "wide"){
        stop("Argument \'format\' must be \"long\" when exporting several types of residuals. \n")
    }

    ## check data
    if(keep.data && !is.null(data) && any(colnames(data) %in% name.residual)){
        stop("Argument \'data\' should not contain a column named \"",paste(name.residual[name.residual %in% colnames(data)], collapse = "\" \""),"\". \n",
             "This name is used to export the residuals. \n")
    }

    ## ** update design
    if(any(c("partial","partial-center") %in% type.residual)){
        ## extract data and design matrix
        if(is.null(data)){
            data <- stats::model.frame(object)
            data <- data[,setdiff(colnames(data),c("XXindexXX",  "XXstrataXX", "XXstrata.indexXX", "XXtimeXX", "XXtime.indexXX", "XXclusterXX", "XXcluster.indexXX")),drop=FALSE]
            design <- stats::model.matrix(object, effects = effects, simplifies = FALSE)
        }else{
            design <- stats::model.matrix(object, data = data, effects = effects, simplifies = FALSE)
        }

        ## design matrix relative to the data where the effect of variables no in var has been removed
        reference <- NULL
        centering <- NULL
        if("partial" %in% type.residual){ 
            ## set the dataset at the reference value for all variables not in var
            if(is.null(attr.ref)){
                if(length(variableMuFac.name)>0){
                    reference <- c(reference, lapply(xfactorMu,function(iX){factor(iX[1], levels = iX)}))
                }
                if(length(variableMuNum.name)>0){
                    reference <- c(reference, as.list(stats::setNames(rep(0, length(variableMuNum.name)), variableMuNum.name)))
                }
                reference <- data.frame(reference, stringsAsFactors = FALSE)
            }else if(!is.data.frame(attr.ref)){
                stop("Attribute \'reference\' must inherit from data.frame. \n")
            }else{
                reference <- attr.ref
            }
            resdata <- data
            if(length(setdiff(variableMu.name,var))>0){
                for(iVar in setdiff(variableMu.name,var)){
                    if(is.factor(data[[iVar]]) && !is.factor(reference[[iVar]])){
                        stop("The reference value of variable ",iVar," should be a factor. \n")
                    }
                    if(is.factor(data[[iVar]]) && !identical(levels(reference[[iVar]]),levels(data[[iVar]]))){
                        stop("Levels of the reference value of variable \'",iVar,"\' should match those of the original data. \n",
                             "Levels: \"",paste(levels(data[[iVar]]), collapse = "\" \""),"\"\n")
                    }
                    data[[iVar]] <- reference[[iVar]]
                }
            }
            ## build design matrix
            design2 <- stats::model.matrix(object, data = data, effects = effects, simplifies = FALSE)
            if(!is.null(index.na)){ 
                design2$mean <- design2$mean[-index.na,,drop=FALSE]
            }
            ## handle intercept term
            if(keep.intercept==FALSE && "(Intercept)" %in% colnames(design2$mean)){
                design2$mean[,"(Intercept)"] <- 0
            }
        }else if("partial-center" %in% type.residual){
            centering <- colMeans(X.mean)
            design2 <- design
            design2$mean <- sweep(design$mean, FUN = "-", MARGIN = 2, STATS = centering)
        }
        
    }else{
        design <- stats::model.matrix(object, data = data, effects = effects, simplifies = FALSE)
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
    index.cluster <- attr(design$index.cluster,"vectorwise")
    index.time <- attr(design$index.clusterTime,"vectorwise")
    index.variance <- as.character(structure$X$pattern.cluster$pattern)
    n.pattern <-  NROW(structure$X$Upattern)

    ## ** update Omega
    if(!is.null(p)){
        if(any(duplicated(names(p)))){
            stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
        }
        if(any(names(param.type) %in% names(p) == FALSE)){
            stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(param.type)[names(param.type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
        }
        beta <- p[names(which(param.type=="mu"))]
        if(any(type.residual %in% c("studentized","pearson","normalized","normalized2","scaled"))){
            Omega <- .calc_Omega(object = structure, param = p)
            precision <- lapply(Omega, solve)
        }
    }else{
        beta <- param.value[param.type=="mu"]
        if(any(type.residual %in% c("studentized","pearson","normalized","normalized2","scaled"))){
            Omega <- object.Omega
            precision <- object.OmegaM1
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
                    iX <- precompute.XX$pattern[[iPattern]]
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

        ## Var(AX) = A VarX A' = A \Omega A' = A \Omega^{1/2} \Omega^{t/2} A'
        ## so A = \Omega^{-1/2} by first taking Cholesky and then inverting
        ## However later on we use X A i.e. X[1,] %*% A i.e. t(A) %*% t(X[1,])
        ## but chol(B) = A' A in R (instead of A A') so we do have A = solve(chol(\Omega))
        sqrtPrecision$normalized <- lapply(Omega,function(iP){solve(chol(iP))})
    }

    ## ** raw residuals
    if(!is.null(data) || !is.null(p) || "partial" %in% type.residual || "partial-center" %in% type.residual){
        fitted <- (X %*% beta)[,1]
        res <-  as.vector(Y - fitted)
        M.res <- matrix(NA, nrow = NROW(X), ncol = length(type.residual), dimnames = list(NULL, name.residual))
    }else{
        fitted <- stats::fitted(object)
        if(!is.null(index.na)){ ## make sure to be consistent with X %*% beta
            fitted <- fitted[-index.na]
        }
        res <- as.vector(object.residuals)
        M.res <- matrix(NA, nrow = NROW(object.residuals), ncol = length(type.residual), dimnames = list(NULL, name.residual))
    }

    ## ** normalization
    if ("response" %in% type.residual) {
        M.res[,"r.response"] <- res
    }
    if ("partial" %in% type.residual){
        M.res[,"r.partial"] <- design2$mean %*% beta + res
        attr(M.res,"reference") <- reference
    }
    if ("partial-center" %in% type.residual){
        index.var <- which(attr(X.mean,"variable") %in% var)
        index.col <- which(attr(X.mean,"assign") %in% index.var)
        M.res[,"r.partial"] <- design2$mean[,index.col,drop=FALSE] %*% beta[index.col] + res
        attr(M.res,"centering") <- centering
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
                    resnorm <- as.double(iResidual %*% solve(chol(Omega[[index.variance[iId]]] - iQ)))
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

    ## ** restaure NA
    if(is.null(match.call()$data) && length(index.na)>0){
        n.allobs <- NROW(M.res)+length(index.na)

        Msave.res <- M.res
        M.res <- matrix(NA, nrow = n.allobs, ncol = NCOL(Msave.res), dimnames = list(NULL, colnames(Msave.res)))
        M.res[-index.na,] <- Msave.res
        attr(M.res,"centering") <- attr(Msave.res,"centering")
        attr(M.res,"reference") <- attr(Msave.res,"reference")

        Msave.fit <- fitted
        fitted <- matrix(NA, nrow = n.allobs, ncol = NCOL(Msave.fit), dimnames = list(NULL, colnames(Msave.fit)))
        fitted[-index.na,] <- Msave.fit

        if(is.null(cluster.levels.original)){
            level.cluster <- levels(attr(index.na,"cluster"))
        }else{
            level.cluster <- rep(NA, n.allobs)
            level.cluster[index.na] <- factor(attr(index.na,"cluster"), levels = cluster.levels.original)
            level.cluster[-index.na] <- factor(index.cluster, levels = 1:length(cluster.levels.original), labels = cluster.levels.original)
        }

        if(is.null(cluster.levels.original)){
            level.time <- levels(attr(index.na,"time"))
        }else{
            level.time <- rep(NA, n.allobs)
            if(NCOL(attr(index.na,"time"))>1){
                level.time[index.na] <- factor(interaction(attr(index.na,"time")), levels = U.time)
                level.time[-index.na] <- factor(index.time, levels = 1:length(U.time), labels = U.time) ## time.levels.original
            }else{
                level.time[index.na] <- factor(attr(index.na,"time"), levels = U.time)
                level.time[-index.na] <- factor(index.time, levels = 1:length(U.time), labels = U.time) ## time.levels.original
            }
            
        }
    }else{
        level.cluster <- index.cluster
        level.time <- index.time
    }

    ## ** export
    if(format=="wide"){
        dfL.res <- data.frame(residuals = as.vector(M.res), cluster = level.cluster, time = factor(U.time[level.time], U.time), stringsAsFactors = FALSE)
        out <- stats::reshape(data = dfL.res[,c("cluster","time","residuals"),drop=FALSE], 
                              direction = "wide", timevar = "time", idvar = "cluster", v.names = "residuals")
        colnames(out)[-1] <- gsub("^residuals.","",colnames(out)[-1])
        if(n.time==1){
            names(out)[-1] <- name.residual
        }else{
            names(out)[-1] <- paste0(name.residual,".",names(out)[-1])
        }
    }else if(format == "long"){

        if(keep.data){
            out <- cbind(data, fitted = fitted, M.res)
        }else if(length(type.residual)==1){
            out <- as.vector(M.res)
        }else{
            out <- M.res
        }


    }

    ## ** export
    attr(out,"reference") <- attr(M.res,"reference")
    attr(out,"centering") <- attr(M.res,"centering")
    attr(out,"args") <- list(type = type, format = format, keep.data = keep.data, var = var, n.time = n.time, name.time = name.time)
    if(format == "long"){
        attr(out,"args")$name.colres <- colnames(M.res)
    }else if(format == "wide"){
        attr(out,"args")$name.colres <- names(out)[-1]
    }
    class(out) <- append("residuals_lmm",class(out))
    return(out)
}

## NOTE: normalizing the residuals can be done by inverting, taking the cholesky transform, possibly transpose
##       but it is not clear in which order to do that
##       here is a small simulation study indicating the "best" solution
##
## rho <- 0.8
## Rho <- matrix(c(1,rho,rho^2,rho^3,
##                 rho,1,rho,rho^2,
##                 rho^2,rho,1,rho,
##                 rho^3,rho^2,rho,1),4,4)
## Sigma <- tcrossprod(1:4) * Rho
##
## library(mvtnorm)
## X <- rmvnorm(100000, mean = rep(0,4), sigma = Sigma)
## var(X)
##
## quantile(var(X %*% solve(t(chol(Sigma)))) - diag(1,4)) ## -2.00926751 -0.78138946 -0.26772457 -0.03288554  2.86152941 
## quantile(var(X %*% solve(chol(Sigma))) - diag(1,4))    ## -0.0022343367 -0.0009080941  0.0013634135  0.0041065544  0.0080840791 
## quantile(var(X %*% chol(solve(Sigma))) - diag(1,4))    ## -0.5963488 -0.1940092  0.3265089  0.6047354  1.7839809 
## quantile(var(X %*% t(chol(solve(Sigma)))) - diag(1,4)) ## -0.0044414902 -0.0004810946  0.0003719232  0.0014844712  0.0098206876

## tSigmaM12 <- t(chol(solve(Sigma)))
## SigmaM12 <- chol(solve(Sigma))
## (X %*% tSigmaM12)[1,]
## X[1,,drop=FALSE] %*% tSigmaM12
## SigmaM12 %*% X[1,]
## 
## quantile(var(t(SigmaM12 %*% t(X))) - diag(1,4))

## * residuals.clmm
##' @export
##' @rdname residuals
residuals.clmm <- function(object, ...){

    object$residuals <- object$design$Y - stats::predict(object, newdata = object$data.original, se = FALSE)$estimate
    object$Omega <- .calc_Omega(object$design$vcov, param = object$param, keep.interim = FALSE)
    object$OmegaM1 <- lapply(object$Omega, solve)
    out <- residuals.lmm(object, ...)
    return(out)

}

## * residuals.mlmm (code)
##' @export
##' @rdname residuals
residuals.mlmm <- function(object, keep.by = FALSE, ...){
    
    if(identical(keep.by,"only")){
        id.col <- "by"
    }else if(keep.by){
        id.col <- c("by",object$object$by)
    }else{
        id.col <- object$object$by
    }
        
    out <- lapply(names(object$model), function(iBy){ ## iBy <- "A"
        iRes <- residuals(object$model[[iBy]], ...)
        if(is.data.frame(iRes)==FALSE){
            iRes <- data.frame(residuals = iRes)
        }
        if(identical(keep.by,"only") || keep.by){
            iRes$by <- iBy
        }
        if(!identical(keep.by,"only")){
            for(iVar in object$object$by){
                iRes[[iVar]] <- object$univariate[object$univariate$by==iBy,iVar]
            }
        }
        return(iRes[,c(id.col,setdiff(names(iRes),id.col)),drop=FALSE])
    })
    
    return(do.call(rbind,out))
}
##----------------------------------------------------------------------
### residuals.R ends here
