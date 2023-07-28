### residuals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:40) 
## Version: 
## Last-Updated: jul 28 2023 (17:29) 
##           By: Brice Ozenne
##     Update #: 1033
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
##' @param format [character] Should the residuals be output
##' in a matrix format with clusters in row and timepoints in columns (\code{"wide"}),
##' or in a data.frame/vector with as many rows as observations (\code{"long"})
##' @param keep.data [logical] Should the dataset relative to which the residuals are evaluated be output along side the residual values?
##' Only possible in the long format.
##' @param simplify [logical] Simplify the data format (vector instead of data.frame) and column names (no mention of the time variable) when possible.
##' Otherwise, information about the call and reference values used for partial residuals be added as an attribute.
##' 
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details The argument \code{type} defines how the residuals are computed:
##' \itemize{
##' \item \code{"fitted"}: fitted value \eqn{X_{ij} \hat{\beta}}.
##' \item \code{"response"}: raw residual, i.e. observed outcome minus fitted value \eqn{\varepsilon_{ij} = Y_{ij} - X_{ij} \hat{\beta}}.
##' \item \code{"pearson"}: each raw residual is divided by its modeled standard deviation \eqn{\varepsilon_{ij} = \frac{Y_{ij} - X_{ij} \hat{\beta}}{\sqrt{\hat{\omega}_{ij}}}}.
##' \item \code{"studentized"}: same as \code{"pearson"} but excluding the contribution of the cluster in the modeled standard deviation  \eqn{\varepsilon_{ij} = \frac{Y_{ij} - X_{ij} \hat{\beta}}{\sqrt{\hat{\omega}_{ij}-\hat{q}_{ij}}}}.
##' \item \code{"normalized"}: raw residuals are multiplied, within clusters, by the inverse of the (upper) Cholesky factor of the modeled residual variance covariance matrix \eqn{\varepsilon_{ij} = ( Y_{i} - X_{i} \hat{\beta} )\hat{C}^{-1}}.
##' \item \code{"normalized2"}: raw residuals are multiplied, within clusters, by the inverse of the modeled residual variance covariance matrix \eqn{\varepsilon_{ij} = ( Y_{i} - X_{i} \hat{\beta} )\hat{Omega}^{-1}}.
##' \item \code{"scaled"}: scaled residuals (see PROC MIXED in SAS). Numerically identical to \code{"normalized"} but computed by sequentially scaling and centering the residuals, to make them conditionally independent of previous residuals from the same cluster at previous repetitions.
##' \item \code{"partial"}: partial residuals (\eqn{\gamma E + \hat{\varepsilon}}). A reference level can be also be specified via the attribute \code{"reference"} to change the absolute level of the partial residuals.
##' \code{"partial-center"}: partial residuals with centered continuous covariates (\eqn{\gamma E + \hat{\varepsilon}} where \eqn{E} has been centered, i.e., has 0-mean)
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
##' \bold{lmm}: a vector or a data.frame when \code{format="long"} (one line per observation, one column per type of residual),
##' a matrix when \code{format="wide"}  (one line per cluster, one column per timepoint).
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
                          data = NULL, p = NULL, format = "long", keep.data = FALSE, simplify = TRUE, ...){

    mycall <- match.call()
    options <- LMMstar.options()
    type.residual <- type
    sep <- options$sep["residuals"]
    
    ## ** extract from object
    xfactorMu <- object$xfactor$mean
    variableMu.name <- attr(object$design$mean,"variable")
    variableMuFac.name <- names(xfactorMu)
    variableMuNum.name <- setdiff(variableMu.name,variableMuFac.name)

    param.name <- object$design$param$name
    param.type <- stats::setNames(object$design$param$type,param.name)
    param.value <- object$param

    n.time <- object$time$n
    index.na <- object$index.na
    X.mean <- object$design$mean
    object.Omega <- object$Omega
    object.OmegaM1 <- object$OmegaM1

    object.cluster <- object$cluster
    object.time <- object$time
    object.index.na <- object$index.na
    U.cluster <- object$cluster$levels
    U.time <- object$time$levels
    name.time <- object$time$var
    
    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    ## check format
    format[] <- match.arg(sort(unique(format)), c("wide","long"), several.ok = TRUE)  ## use 'format[] <-' instead of 'format <-' to keep the name that will be transferd to .reformat(
    if(length(format)>1){
        format <- format[1]
        attr(format,"original") <- c("wide","long")
    }
    if(format=="wide" && identical("all",tolower(type.residual))){
        message("Move to wide format to output all types of residuals. \n")
        format <- "long"
    }
    if(keep.data && format == "wide"){
        stop("Argument \'keep.data\' must be \"FALSE\" when using the wide format. \n")
    }
    ## check type.residuals
    if(identical("all",tolower(type.residual))){        
        type.residual <- c("response","studentized","pearson","normalized")
    }
    attr.ref <- attr(type.residual,"reference")
    valid.normresiduals  <- c("studentized","pearson","normalized","normalized2","normastudentized","scaled")
    valid.residuals <- c("response",valid.normresiduals,"partial","partial-center")
    if(any(grepl("-ci$", type.residual))){
        type.residual <- gsub("-ci$","",type.residual)
        ci <- TRUE ## hidden ci argument
    }else{
        ci <- FALSE
    }
    type.residual <- match.arg(type.residual, valid.residuals, several.ok = (format=="long"))
    if(any(type.residual %in% valid.normresiduals)){
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
        if(!is.null(attr.ref) && "partial-center" %in% type.residual){
            message("Attribute \'attr.ref\' of argument \'type\' is ignored when \'type\' equals \"partial-center\". \n")
            attr.ref <- NULL
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
        type.var <- c("numeric","categorical")[var %in% names(object$xfactor$mean) + 1]
        type.fit <- ifelse(keep.intercept,"static","static0")
    }else{
        if(!is.null(var)){
            warning("Argument \'var\' ignored when computing residuals other than partial residuals. \n")
        }
        keep.intercept <- TRUE
        type.var <- NULL        
        var <- NULL
        type.fit <- "static"
    }

    ## check agreement plot, format, type.residual
    if(length(type.residual)>1 && format == "wide"){
        stop("Argument \'format\' must be \"long\" when exporting several types of residuals. \n")
    }

    ## check data and create data.reference used for export and for partial residuals 
    if(!is.null(data)){
        data.reference <- as.data.frame(data)
        if(keep.data && any(colnames(data.reference) %in% name.residual)){
            stop("Argument \'data\' should not contain a column named \"",paste(name.residual[name.residual %in% colnames(data.reference)], collapse = "\" \""),"\". \n",
                 "This name is used to export the residuals. \n")
        }        
    }else{
        data.reference <- object$data.original
    }

    ## ** update design
    if("partial" %in% type.residual || "partial-center" %in% type.residual){
        ## extract data and design matrix
        if(is.null(data)){
            design <- stats::model.matrix(object, effects = effects, simplify = FALSE)
        }else{
            design <- stats::model.matrix(object, data = data, effects = effects, simplify = FALSE)
        }

        ## *** design matrix relative to a reference value
        reference <- stats::setNames(as.list(rep(NA, length = length(variableMu.name))), variableMu.name)
        ## reference: reference level of all variables not in var
        if(is.null(attr.ref)){
            if(length(setdiff(variableMuFac.name,var))>0){
                reference[setdiff(variableMuFac.name,var)] <- lapply(setdiff(variableMuFac.name,var),function(iName){factor(xfactorMu[[iName]][1], levels = xfactorMu[[iName]])})
            }
            if(length(setdiff(variableMuNum.name,var))>0){
                reference[setdiff(variableMuNum.name,var)] <- as.list(stats::setNames(rep(0, length(setdiff(variableMuNum.name,var))), setdiff(variableMuNum.name,var)))
            }
            if("partial-center" %in% type.residual){
                ## do nothing for categorical variables

                ## center numeric variables
                reference[intersect(variableMuNum.name, var)] <- lapply(intersect(variableMuNum.name, var), function(iName){mean(object$data[[iName]])})
            }
            reference <- data.frame(reference, stringsAsFactors = FALSE)
        }else if(!is.data.frame(attr.ref)){
            stop("Attribute \'reference\' form argument \'type\' must inherit from data.frame. \n")
        }else{
            reference <- attr.ref
            if(NROW(reference)!=1){
                stop("Attribute \'reference\' from argument \'type\' must be a single row data.frame. \n")
            }
            if(any(names(reference) %in% variableMu.name == FALSE)){
                stop("Incorrect column names in attribute \'reference\' from argument \'type\'. \n",
                     "Valid column names: \"",paste(variableMu.name, collapse = "\", \""),"\". \n")
            }
        }

        ## apply reference
        reference.effective <- lapply(reference, function(iRef){if(all(is.na(iRef))){NULL}else{iRef}})
        reference.effective <- as.data.frame(reference.effective[lengths(reference.effective)>0])
        for(iVar in names(reference.effective)){

            if(is.factor(data.reference[[iVar]]) && !is.factor(reference[[iVar]])){
                stop("The reference value of variable ",iVar," should be a factor. \n")
            }
            if(is.factor(data.reference[[iVar]]) && !identical(levels(reference[[iVar]]),levels(data.reference[[iVar]]))){
                stop("Levels of the reference value of variable \'",iVar,"\' should match those of the original data. \n",
                     "Levels: \"",paste(levels(data[[iVar]]), collapse = "\" \""),"\"\n")
            }
            if(iVar %in% var){
                data.reference[[iVar]] <- data.reference[[iVar]] - reference[[iVar]]
            }else{
                data.reference[[iVar]] <- reference[[iVar]]
            }
        }
        
        ## build design matrix
        design.reference <- stats::model.matrix(object, data = data.reference, effects = effects, simplify = TRUE)
        ## handle intercept term
        if(keep.intercept==FALSE && "(Intercept)" %in% colnames(design.reference)){
            design.reference[,"(Intercept)"] <- 0
        }
    }else{
        design <- stats::model.matrix(object, data = data, effects = effects, simplify = FALSE)
    }
    Y <- design$Y
    X <- design$mean
    structure <- design$vcov
    n.cluster <- length(design$index.cluster)
    precompute.XX <- design$precompute.XX
    index.cluster <- attr(design$index.cluster,"vectorwise")
    index.time <- attr(design$index.clusterTime,"vectorwise")
    pattern <- as.character(structure$pattern)
    n.pattern <-  NROW(structure$Upattern)

    ## ** update Omega
    if(!is.null(p)){
        if(any(duplicated(names(p)))){
            stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
        }
        if(any(names(param.type) %in% names(p) == FALSE)){
            stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(param.type)[names(param.type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
        }
        beta <- p[names(which(param.type=="mu"))]
        if(any(type.residual %in% valid.normresiduals)){
            Omega <- .calc_Omega(object = structure, param = p)
            precision <- lapply(Omega, solve)
        }
    }else{
        beta <- param.value[param.type=="mu"]
        if(any(type.residual %in% valid.normresiduals)){
            Omega <- object.Omega
            precision <- object.OmegaM1
        }
    }

    ## ** pre-compute
    sqrtPrecision <- list()
    if("pearson" %in% type.residual){
        sqrtPrecision$pearson <- lapply(Omega,function(iM){1/sqrt(diag(iM))})
    }
    if("studentized" %in% type.residual || "normastudentized" %in% type.residual){
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
                    tX.precision.X <- tX.precision.X + t(iX) %*% precision[[pattern[iId]]] %*% iX
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
    if(!is.null(data) || !is.null(p)){
        fitted <- (X %*% beta)[,1]
        res <-  as.vector(Y - fitted)
        M.res <- matrix(NA, nrow = NROW(X), ncol = length(type.residual), dimnames = list(NULL, name.residual))
    }else{
        fitted <- object$fitted
        res <- object$residuals
        M.res <- matrix(NA, nrow = length(res), ncol = length(type.residual), dimnames = list(NULL, name.residual))
    }
    if(ci || type.residual %in% c("partial","partial-center")){

        df.fitted <- predict(object, newdata = design.reference, type = type.fit, se = ifelse(ci,"estimation",FALSE),
                             keep.newdata = FALSE, format = "long", simplify = FALSE)

        if(ci){
            fitted <- cbind(fitted = df.fitted$estimate,
                            fitted.lower = df.fitted$lower,
                            fitted.upper = df.fitted$upper)
        }else{
            fitted <- df.fitted$estimate
        }
    }

    ## ** normalization
    if ("response" %in% type.residual) {
        M.res[,"r.response"] <- res
    }
    if ("partial" %in% type.residual || "partial-center" %in% type.residual){
        M.res[,"r.partial"] <- design.reference %*% beta + res
        attr(M.res,"reference") <- reference
    }
    if (any(type.residual %in% valid.normresiduals)) {
        for(iId in 1:n.cluster){ ## iId <- 7
            
            iIndex <- which(index.cluster==iId)
            iOrder <- order(index.time[iIndex])
            iResidual <- res[iIndex[iOrder]]
            iN.time <- length(iIndex)

            for(iType in type.residual){
                if("pearson" %in% type.residual){
                    resnorm <- iResidual * sqrtPrecision$pearson[[pattern[iId]]]
                    M.res[iIndex[iOrder],"r.pearson"] <- resnorm
                }
                if("studentized" %in% type.residual){
                    iX <- X[iIndex[iOrder],,drop=FALSE]
                    iQ <- iX %*% tX.precision.X.M1 %*% t(iX)
                    resnorm <- iResidual / sqrt(diag(Omega[[pattern[iId]]] - iQ))
                    M.res[iIndex[iOrder],"r.studentized"] <- resnorm
                }
                if("normalized" %in% type.residual){
                    ## resnorm <- as.double(sqrtPrecision$normalized[[pattern[iId]]] %*% iResidual)
                    resnorm <- as.double(iResidual %*% sqrtPrecision$normalized[[pattern[iId]]])
                    M.res[iIndex[iOrder],"r.normalized"] <- resnorm
                }
                if("normalized2" %in% type.residual){
                    resnorm <- as.double(iResidual %*% precision[[pattern[iId]]])
                    M.res[iIndex[iOrder],"r.normalized2"] <- resnorm
                }
                if("normastudentized" %in% type.residual){
                    iX <- X[iIndex[iOrder],,drop=FALSE]
                    iQ <- iX %*% tX.precision.X.M1 %*% t(iX)
                    resnorm <- as.double(iResidual %*% solve(chol(Omega[[pattern[iId]]] - iQ)))
                    M.res[iIndex[iOrder],"r.normalized2"] <- resnorm
                }
                if("scaled" %in% type.residual){
                    M.res[iIndex[iOrder][1],"r.scaled"] <- iResidual[1]/attr(Omega[[pattern[iId]]],"sd")[1]
                    if(iN.time>1){
                        for(iTime in 2:iN.time){ ## iTime <- 2
                            iVar <- Omega[[pattern[iId]]][iTime,iTime]
                            iPrecision_kk <- solve(Omega[[pattern[iId]]][1:(iTime-1),1:(iTime-1),drop=FALSE])
                            iOmega_lk <- Omega[[pattern[iId]]][iTime,1:(iTime-1),drop=FALSE]
                            iOmega_kl <- Omega[[pattern[iId]]][1:(iTime-1),iTime,drop=FALSE]
                                
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
    if(is.null(mycall$data)){
               
        M.res <- addNA(M.res, index.na = index.na,
                       level = "obs", cluster = object$cluster)
        
        if(keep.data){

            if(format == "long"){
                fitted <- addNA(fitted, index.na = index.na, 
                                level = "obs", cluster = object$cluster)

                if(is.matrix(fitted)){ ## when using partial residuals with ci for the fitted values
                    data.reference <- cbind(data.reference, fitted)
                }else{ ## normal case
                    data.reference <- cbind(data.reference, fitted = fitted)
                }
            }
        }

    }

    ## ** export
    out <- .reformat(M.res, name = names(format), format = format, simplify = simplify,
                     keep.data = keep.data, data = data.reference, index.na = object.index.na,
                     object.cluster = object.cluster, index.cluster = index.cluster,
                     object.time = object.time, index.time = index.time,                     
                     call = mycall)


    if(!simplify){
        
        attr(out,"reference") <- attr(M.res,"reference")
        attr(out,"centering") <- attr(M.res,"centering")
        attr(out,"args") <- list(type = type, format = format, keep.data = keep.data, var = var, type.var = type.var, n.time = n.time, name.time = name.time,
                                 outcome = object$outcome$var, intercept = "(Intercept)" %in% names(object$param))
        if(keep.intercept){
            attr(out,"args")$var <- c("(Intercept)",attr(out,"args")$var)
        }
        if(all(is.na(attr(object.time$var,"original")))){
            attr(out,"index.time") <- addNA(as.character(index.time), index.na = index.na,
                                            level = "obs", cluster = object$cluster)
        }
        attr(out,"args")$nameL.colres <- colnames(M.res)
        if(format == "wide"){
            attr(out,"args")$nameW.colres <- stats::setNames(names(out)[-1], object.time$levels)

        }else if(format == "long" && !is.null(attr(format,"original"))){
            attr(out,"wide") <- .reformat(M.res, name = names(format), format = "wide", simplify = simplify,
                                          keep.data = FALSE, data = data, index.na = object.index.na,
                                          object.cluster = object.cluster, index.cluster = index.cluster,
                                          object.time = object.time, index.time = index.time,                     
                                          call = mycall)
            attr(out,"args")$nameW.colres <- stats::setNames(names(attr(out,"wide"))[-1], object.time$levels)
        }
    }
    class(out) <- append("residuals_lmm",class(out))
    return(out)
}

## * residuals.clmm (code)
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
residuals.mlmm <- function(object, simplify = TRUE, ...){

    ## ** extract
    ls.out <- lapply(names(object$model), function(iBy){ ## iBy <- "A"
        residuals(object$model[[iBy]], simplify = simplify, ...)
    })

    ## ** reshape
    test.2D <- any(sapply(ls.out, inherits, "data.frame")) || any(sapply(ls.out, inherits, "matrix"))
    if(test.2D){
        out <- do.call(rbind,ls.out)
        if(simplify && is.data.frame(out)){
            object.manifest <- manifest(object)
            rm.manifest <- attributes(object.manifest)[setdiff(names(attributes(object.manifest)),c("by","cluster","time"))]
            out[unique(unlist(rm.manifest))] <- NULL
        }
    }else{
        out <- do.call(c,ls.out)
    }

    ## ** export
    return(out)
}

## * Note about normalizing the residuals
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


##----------------------------------------------------------------------
### residuals.R ends here
