### predict.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: May 12 2024 (18:23) 
##           By: Brice Ozenne
##     Update #: 1443
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * predict.lmm (documentation)
##' @title Predicted Mean Value With Uncertainty For Linear Mixed Model
##' @description Predicted mean value conditional on covariates or on covariates and other outcome values.
##'
##' @param object a \code{lmm} object.
##' @param newdata [data.frame] a dataset containing covariate values to condition on.
##' When setting the argument 'dynamic' predictions should also contain cluster, timepoint, and outcome values.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the predictions. Only relevant if differs from the fitted values.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param se [logical] should the standard error and confidence intervals for the predictions be output?
##' It can also be a logical vector of length 2 to indicate the type of uncertainty to be accounted for: estimation and residual variance.
##' In particular \code{c(TRUE,TRUE)} provides prediction intervals.
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors.
##' Not feasible for dynamic predictions when using REML.
##' @param vcov [logical] should the variance-covariance matrix of the predictions be output as an attribute.
##' @param df [logical] should a Student's t-distribution be used to model the distribution of the predicted mean. Otherwise a normal distribution is used.
##' @param type [character] evaluate the expected outcome conditional on covariates only (\code{"static"}),
##' the contribution of each variable to this 'static' prediction (\code{"terms"}),
##' or the expected outcome conditional covariates and outcome values at other timepoints (\code{"dynamic"}).
##' Based on the observed outcome and the 'dynamic' prediction for the missing outcome,
##' can also evaluate the change from first repetitition (\code{"change"}), area under the curve (\code{"auc"}), and the area under the curve minus baseline (\code{"auc-b"}).
##' @param format [character] should the prediction be output
##' in a matrix format with clusters in row and timepoints in columns (\code{"wide"}),
##' or in a data.frame/vector with as many rows as observations (\code{"long"})
##' @param keep.data [logical] should the dataset relative to which the predicted means are evaluated be output along side the predicted values?
##' Only possible in the long format.
##' @param export.vcov [logical] should the variance-covariance matrix of the prediction error be outcome as an attribute (\code{"vcov"})?
##' @param simplify [logical] simplify the data format (vector instead of data.frame) and column names (no mention of the time variable) when possible.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details Static prediction are made using the linear predictor \eqn{X\beta} while dynamic prediction uses the conditional normal distribution of the missing outcome given the observed outcomes. So if outcome 1 is observed but not 2, prediction for outcome 2 is obtain by \eqn{X_2\beta + \sigma_{21}\sigma^{-1}_{22}(Y_1-X_1\beta)}. In that case, the uncertainty is computed as the sum of the conditional variance \eqn{\sigma_{22}-\sigma_{21}\sigma^{-1}_{22}\sigma_{12}} plus the uncertainty about the estimated conditional mean (obtained via delta method using numerical derivatives).
##'
##' The model terms are computing similarly to \code{stats::predict.lm}, by centering the design matrix around the mean value of the covariates used to fit the model.
##' Then the centered design matrix is multiplied by the mean coefficients and columns assigned to the same variable (e.g. three level factor variable) are summed together.
##' 
##' @return When \code{format="long"}, a data.frame containing the following columns:\itemize{
##' \item \code{estimate}: predicted mean.
##' \item \code{se}: uncertainty about the predicted mean.
##' \item \code{df}: degree of freedom
##' \item \code{lower}: lower bound of the confidence interval of the predicted mean
##' \item \code{upper}: upper bound of the confidence interval of the predicted mean
##' }
##' When \code{format="wide"}, a matrix containing the predict means (one line per cluster, one column per timepoint).
##' 
##' @keywords method
##' 
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit Linear Mixed Model
##' eUN.lmm <- lmm(Y ~ visit + X1 + X2 + X5,
##'                repetition = ~visit|id, structure = "UN", data = dL)
##'
##' ## prediction
##' newd <- data.frame(X1 = 1, X2 = 2, X5 = 3, visit = factor(1:3, levels = 1:3))
##' predict(eUN.lmm, newdata = newd)
##' predict(eUN.lmm, newdata = newd, keep.data = TRUE)
##' predict(eUN.lmm, newdata = newd, keep.data = TRUE, se = c(TRUE,TRUE))
##'
##' ## dynamic prediction
##' newd.d1 <- cbind(newd, Y = c(NA,NA,NA))
##' predict(eUN.lmm, newdata = newd.d1, keep.data = TRUE, type = "dynamic")
##' newd.d2 <- cbind(newd, Y = c(6.61,NA,NA))
##' predict(eUN.lmm, newdata = newd.d2, keep.data = TRUE, type = "dynamic")
##' newd.d3 <- cbind(newd, Y = c(1,NA,NA))
##' predict(eUN.lmm, newdata = newd.d3, keep.data = TRUE, type = "dynamic")
##' newd.d4 <- cbind(newd, Y = c(1,1,NA))
##' predict(eUN.lmm, newdata = newd.d4, keep.data = TRUE, type = "dynamic")

## * predict.lmm (code)
##' @export
predict.lmm <- function(object, newdata, type = "static", p = NULL,
                        se = NULL, robust = FALSE, df = NULL, 
                        level = 0.95, keep.data = NULL, format = "long", export.vcov = FALSE, simplify = TRUE, ...){

    ## ** extract from object
    mycall <- match.call()
    options <- LMMstar.options()
    name.Y <- object$outcome$var
    U.time <- object$time$levels
    name.time <- attr(object$time$var,"original")
    name.cluster <- attr(object$cluster$var,"original")
    object.mean <- object$design$mean

    table.param <- object$design$param
    name.mu <- table.param$name[table.param$type=="mu"]
    index.vcov <- which(table.param$type %in% c("sigma","k","rho"))
    name.theta <- table.param$name
    n.theta <- length(name.theta)


    ## ** normalize user imput
    dots <- list(...)

    ## hidden arguments
    if("transform.sigma" %in% names(dots)){
        transform.sigma <- dots$transform.sigma
        dots$transform.sigma <- NULL
    }else{
        transform.sigma <- NULL
    }
    if("transform.k" %in% names(dots)){
        transform.k <- dots$transform.k
        dots$transform.k <- NULL
    }else{
        transform.k <- NULL
    }
    if("transform.rho" %in% names(dots)){
        transform.rho <- dots$transform.rho
        dots$transform.rho <- NULL
    }else{
        transform.rho <- NULL
    }

    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    ## type of prediction
    type.prediction <- match.arg(type, c("static0","static","terms","dynamic","change","auc","auc-b")) ## static0: no intercept
    if((type.prediction == "dynamic") && ("rho" %in% table.param$type == FALSE)){
        type.prediction <- "static"
        message("Move to static prediction as there is no correlation parameter in the model. \n")
    }
    if(type.prediction=="static0"){
        type.prediction <- "static"
        keep.intercept <- FALSE
    }else{
        keep.intercept <- TRUE
        n.time <- object$time$n
        if(n.time==1 && type %in% c("change","auc","auc-b")){
            stop("Cannot evaluate ",type," when there is a single timepoint. \n",
                 "Considering setting argument \'type\' to \"static\" or specifying the argument \'repetition\' when calling lmm. \n")
        }
        if(type.prediction == "dynamic"){
            M.contrast <- diag(1, n.time, n.time)
            dimnames(M.contrast) <- list(U.time, U.time)
        }else if(type.prediction == "change"){
            M.contrast <- matrix(c(-1,rep(0,n.time-1)), byrow = TRUE, nrow = n.time, ncol = n.time) + diag(1, n.time, n.time)
            dimnames(M.contrast) <- list(U.time, U.time)
        }else if(type.prediction %in% c("auc","auc-b")){
            time.num <- as.numeric(U.time)
            dtime.num <- diff(c(utils::head(time.num,1),time.num))/2 + diff(c(time.num,utils::tail(time.num,1)))/2
            if(type.prediction == "auc-b"){
                dtime.num[1] <- dtime.num[1] - sum(dtime.num)
            }
            M.contrast <- matrix(dtime.num, byrow = TRUE, nrow = n.time, ncol = n.time)
            dimnames(M.contrast) <- list(U.time, U.time)
            if(any(is.na(time.num))){
                stop("Cannot evaluate the area under the curve of the outcome. \n",
                     "When calling lmm, argument \'repetition\'(=~time|cluster) must contain a numeric time variable. \n",
                     "Or a factor variable whose levels can be converted as numeric")
            }
            if(is.unsorted(time.num)){
                warning("The levels of the time variable do not correspond to numeric values in increasing order. \n",
                        "Can be an issue when evaluating the area under the curve.")
            }
        }                
    }
       
    ## dataset
    if(is.null(newdata)){
        index.na <- object$index.na
        newdata.index.cluster <- attr(object$design$index.cluster,"vectorwise")
        newdata.index.time <- attr(object$design$index.clusterTime,"vectorwise")
        if(is.null(keep.data)){
            keep.data <- FALSE
        }
    }else{
        index.na <- NULL

        if(is.matrix(newdata)){
            if(type %in% c("dynamic","change","auc","auc-b")){
                stop("Argument \'newdata\' cannot be a matrix when asking for dynamic predictions. \n",
                     "It should be a data.frame. \n")
            }
            if(format == "wide"){
                stop("Argument \'newdata\' cannot be a matrix when requesting a wide format. \n",
                     "It should be a data.frame. \n")
            }
            if(keep.data == TRUE){
                stop("Argument \'keep.data\' should be FALSE when argument \'newdata\' is a matrix. \n")
            }
            if(length(se)==2 && se[2]>0){
                stop("Argument \'se\' should be have length 1 or its second element should be FALSE when argument \'newdata\' is a matrix. \n",
                     "(cannot associate observation to the repetition structure with matrix input - consider using data.frame) \n")
            }
            if(is.null(keep.data)){
                keep.data <- FALSE
            }
            ## used by residuals for type = partial-center
        }else{
            newdata <- as.data.frame(newdata)
            if(type.prediction %in% c("dynamic","change","auc","auc-b") == FALSE && (length(se)==2 && se[2] && export.vcov==FALSE) && name.cluster %in% names(newdata) == FALSE ){
                ## add cluster variable if missing and no duplicated time
                if(any(!is.na(name.time)) && all(name.time %in% names(newdata)) && all(duplicated(newdata[name.time])==FALSE)){
                    newdata[[name.cluster]] <- 1
                }else{
                    stop("Incorrect argument 'newdata': missing cluster variable \"",name.cluster,"\". \n")
                }
            }
            if(is.null(keep.data)){
                keep.data <- FALSE
            }
        }
        if(format == "wide"){            
            newdata.design <- stats::model.matrix(object, newdata = newdata, effects = "index", na.rm = FALSE)
            newdata.index.cluster <- attr(newdata.design$index.cluster, "vectorwise")
            newdata.index.time <- attr(newdata.design$index.clusterTime, "vectorwise")        
        }
    }

    ## standard error
    if(is.null(se)){
        se <- c((!simplify || keep.data) && !(format == "wide") && is.null(p), FALSE)
    }else{
        if(!is.logical(se) & !is.numeric(se)){
            stop("Argument \'se\' should be TRUE or FALSE, or a logical vector of length 2. \n")
        }
        if(length(se)==1){
            se <- c(se,FALSE)
        }else if(length(se)!=2){
            stop("Argument \'se\' should have length 1 or 2. \n")
        }
        if(any(se>0) && !is.null(p)){
            stop("Cannot evaluate the uncertainty when the argument \'p\' is specified")
        }    
    }

    ## parameter value and transformation
    init <- .init_transform(p = p, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho,
                            table.param = object$design$param)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform
    if(is.null(p)){
        theta <- object$param
    }else{
        theta <- init$p
    }    
    mu <- theta[name.mu]

    ## degrees of freedom
    if(is.null(df)){
        df <- !is.null(object$df) & se[1]
    }else if(!is.logical(df) & !is.numeric(df)){
        stop("Argument \'df\' should be TRUE or FALSE. \n")
    }else{
        df <- as.logical(df)
        if(se[1]==FALSE && df){
            message("Argument \'df\' ignored when the first element of argument \'se\' is FALSE. \n",
                    "(can only evaluate degrees of freedom relative to the estimation error)\n")
            df <- FALSE
        }
    }
        
    ## check format
    format[] <- match.arg(format, c("wide","long"))  ## use 'format[] <-' instead of 'format <-' to keep the name that will be transferd to .reformat(
    if(keep.data && format == "wide" && is.null(attr(keep.data, "var"))){
        mean.type <- lava::manifest(object, effects = "mean.type")
        attr(keep.data, "var") <- names(mean.type)[mean.type=="baseline"]
    }    
    if(format == "wide" && (df==TRUE||sum(se)>0)){
        if("df" %in% names(mycall) || "se" %in% names(mycall)){
            message("When using the wide format arguments \'se\' and \'df\' are ignored. \n",
                    "(i.e. set to NULL and FALSE, respectively). \n")
        }
        df <- FALSE
        se <- c(FALSE, FALSE)
    }

    ## impute cluster when missing (if static) and unambiguous, i.e. no repeated times (id dynamic)
    if(inherits(newdata,"data.frame") && !is.na(name.cluster)){
        if(type.prediction %in% c("dynamic","change","auc","auc-b") && is.null(newdata[[name.cluster]]) && all(!is.na(name.time)) && all(name.time %in% names(newdata))){
            if(any(duplicated(newdata[,name.time, drop=FALSE]))){
                stop("Duplicated time values found in column ",name.time,".\n",
                     "Consider specifying the cluster variable in argument \'newdata\'. \n")
            }else{
                newdata[[name.cluster]] <- "1"
            }
        }else if(name.cluster %in% names(newdata) && (is.factor(newdata[[name.cluster]]) || is.numeric(newdata[[name.cluster]]))){
            newdata[[name.cluster]] <- as.character(newdata[[name.cluster]])
        }
    }

    ## check data
    if(type.prediction %in% c("dynamic","change","auc","auc-b")){

        if(!is.null(newdata) && name.Y %in% names(newdata) == FALSE){
            stop("The outcome column \"",name.Y,"\" in argument \'newdata\' is missing and necessary when doing dynamic predictions. \n")
        }
        if(all(is.na(name.cluster))){
            stop("The cluster variable should be explicit when doing dynamic predictions. \n",
                 "Consider re-fitting lmm specifying the argument repetition as ~time|cluster. \n")
        }
        if(name.cluster %in% names(newdata) == FALSE){
            stop("The cluster column \"",name.cluster,"\" in argument \'newdata\' is missing and necessary when doing dynamic predictions. \n")
        }
        if(all(is.na(name.time))){
            stop("The time variable should be explicit when doing dynamic predictions. \n",
                 "Consider re-fitting lmm specifying the argument repetition as ~time|cluster. \n")
        }
        if(any(name.time %in% names(newdata) == FALSE)){
            stop("The time column \"",paste(name.time, collapse=", "),"\" in argument \'newdata\' is missing and necessary when doing dynamic predictions. \n")
        }
        test.level <- sapply(name.time, function(iVar){all(newdata[[iVar]] %in% attr(U.time,"original")[,iVar])})
        if(any(test.level == FALSE)){
            stop("The time column \"",names(which(test.level==FALSE)[1]),"\" in argument \'newdata\' should match the existing times. \n",
                 "Existing times: \"",paste0(attr(U.time,"original")[,names(which(test.level==FALSE)[1])], collapse = "\" \""),"\".\n")
        }
        test.duplicated <- by(newdata[,name.time,drop=FALSE],newdata[[name.cluster]], function(iT){any(duplicated(iT))})
        if(any(unlist(test.duplicated))){
            stop("The time column \"",paste(name.time, collapse = "\", \""),"\" in argument \'newdata\' should not have duplicated values within clusters. \n")
        }

    }else if(type.prediction == "static"){
        
        if(se[2] && all(!is.na(name.time)) && any(name.time %in% names(newdata) == FALSE)){
            stop("The time column \"",paste(name.time[name.time %in% names(newdata) == FALSE], collapse = "\", \""),"\" in missing from argument \'newdata\'. \n",
                 "It is necessary for uncertainty calculation involving the residual variance. \n")
        }

    }

    if(keep.data && any(c("estimate","se","df","lower","upper") %in% names(newdata))){
        stop("Argument \'newdata\' should not contain a column named \"estimate\", \"se\", \"lower\", \"upper\", or \"df\" as those names are used internally. \n")
    }

    ## ** parameters
    if(se[2] || type.prediction %in% c("dynamic","change","auc","auc-b")){
        if(is.null(p) && test.notransform){
            reparametrize <- object$reparametrize
        }else{
            reparametrize <- .reparametrize(p = theta[index.vcov], type = table.param$type[index.vcov], level = table.param$level[index.vcov], 
                                            sigma = table.param$sigma[index.vcov], k.x = table.param$sigma[index.vcov], k.y = table.param$sigma[index.vcov],
                                            Jacobian = TRUE, dJacobian = FALSE, inverse = FALSE, 
                                            transform.sigma = transform.sigma,
                                            transform.k = transform.k,
                                            transform.rho = transform.rho,
                                            transform.names = FALSE)
        }
    }

    ## vcov of the estimates
    if(se[1]>0 && type.prediction == "static"){
        vcov.mu <- vcov(object, effects = "mean", p = p, robust = robust)
    }
    if(df || se[2]>0 || type.prediction %in% c("dynamic","change","auc","auc-b")){
        vcov.theta <- vcov(object, effects = "all", p = p, robust = robust, df = 2*df, ## 2* to extract dVcov and no only df
                           transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
        if(df){
            dVcov.theta <- attr(vcov.theta,"dVcov")
            attr(vcov.theta, "df") <- NULL
            attr(vcov.theta, "dVcov") <- NULL
        }        
    }

    ## ** design matrix
    if(is.matrix(newdata)){
        X <- newdata
    }else{
        X <- stats::model.matrix(object, newdata = newdata, effects = "mean", na.rm = FALSE)
    }
    if(!keep.intercept && "(Intercept)" %in% colnames(X)){
        X[,"(Intercept)"] <- 0
    }
    n.obs <- NROW(X)

    ## ** terms        
    if(type.prediction == "terms"){
        Xmean <- colMeans(object.mean)
        Xc <- sweep(X, FUN = "-", MARGIN = 2, STATS = Xmean)
        Xcmu <- sweep(Xc, FUN = "*", MARGIN = 2,  STATS = mu)

        index.n0 <- which(attr(object.mean,"assign")!=0)
        if(length(index.n0)==0){
            Xterm <- matrix(nrow = NROW(newdata), ncol = 0)
        }else{
            Xterm <- do.call(cbind,tapply(index.n0,attr(object.mean,"assign")[index.n0],function(iCol){
                list(rowSums(Xcmu[,iCol,drop=FALSE]))
            }))
            colnames(Xterm) <- attr(object.mean,"term.labels")
        }

        if(any(attr(object.mean,"assign")==0)){
            attr(Xterm, "constant") <- sum(Xmean*mu)
        }
        return(Xterm)
    }

    ## ** identify variance patterns
    if(type.prediction %in% c("dynamic","change","auc","auc-b") || se[2]){

        newdesign <- stats::model.matrix(object, newdata = newdata, effect = "variance", simplify = FALSE, na.rm = FALSE)
        index.cluster <- newdesign$index.cluster
        index.clusterTime <- newdesign$index.clusterTime
        n.cluster <- length(index.cluster)

        Omega <- .calc_Omega(newdesign$vcov, param = theta, keep.interim = TRUE)
        if(se[1]){
            dOmega <- .calc_dOmega(newdesign$vcov, param = theta, Omega = Omega,
                                   Jacobian = reparametrize$Jacobian,
                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)            
        }
        
        if(!is.na(name.cluster) && length(unique(newdata[[name.cluster]]))!=n.cluster){
            stop("Something went wrong when extracting the residual variance-covariance matrices from newdata. \n")
        }

    }

    ## ** prediction with uncertainty
    if(se[1]>0 & (type.prediction %in% c("dynamic","change","auc","auc-b") | !simplify | df) ){
        grad.var <- matrix(0, nrow = n.obs, ncol = n.theta,
                           dimnames = list(NULL, name.theta))
    }else{
        grad.var <- NULL
    }
        

    if(type.prediction == "static"){
        ## compute predictions
        prediction <- (X %*% mu)[,1]

        ## compute uncertainty about the predictions
        if(sum(se)>0){
            prediction.var <- rep(0, times = n.obs)
            if(export.vcov){
                attr(prediction.var,"vcov") <- matrix(0, n.obs, n.obs)
            }
            if(se[1]){ ## same as diag(X %*% vcov.mu %*% t(X))
                prediction.var <- prediction.var + rowSums((X %*% vcov.mu) * X)
                if(!simplify || df){
                    grad.var[,name.mu] <- X
                }
                if(export.vcov){
                    attr(prediction.var,"vcov") <- attr(prediction.var,"vcov") + X %*% vcov.mu %*% t(X)
                }
            }
            if(se[2]){
                prediction.var[unlist(index.cluster)] <- prediction.var[unlist(index.cluster)] + unlist(lapply(Omega,diag))
                if(export.vcov){
                    for(iC in 1:n.cluster){ 
                        attr(prediction.var,"vcov")[index.cluster[[iC]],index.cluster[[iC]]] <- attr(prediction.var,"vcov")[index.cluster[[iC]],index.cluster[[iC]]] + Omega[[iC]]
                    }
                }
            }
                        
        }else{
            prediction.var <- rep(NA, times = n.obs)
        }

        
    }else if(type.prediction %in% c("dynamic","change","auc","auc-b")){

        ## prepare
        prediction <- rep(NA, times = n.obs)
        prediction.var <- rep(NA, times = n.obs)
        if(se[1]>0){
            grad.var <- matrix(0, nrow = n.obs, ncol = n.theta,
                               dimnames = list(NULL, name.theta))
        }
        if(export.vcov){
            attr(prediction.var,"vcov") <- matrix(0, n.obs, n.obs) ## ok when factor.residual but should be NA otherwise where no inputation is used. This is done below 

        }
        iX <- matrix(0, nrow = n.time, ncol = length(mu), dimnames = list(NULL, name.mu))
        iY <- stats::setNames(rep(0, n.time), U.time)
                
        for(iC in 1:n.cluster){ ## iC <- 1

            iNewdata <- newdata[index.cluster[[iC]],,drop=FALSE] ## subset, making sure that time are in increasing order via index.cluster
            iOmega <- Omega[[newdesign$vcov$pattern[iC]]]

            iIndex.con <- which(!is.na(iNewdata[[name.Y]])) ## position of the condition set in the specific dataset
            iTime.con <- index.clusterTime[[iC]][iIndex.con] ## timepoints of the condition set in the specific dataset
            iPos.con <- index.cluster[[iC]][iIndex.con] ## position of the condition set in the original dataset
            iX.con <- X[iPos.con,,drop=FALSE]
            
            iIndex.pred <- which(is.na(iNewdata[[name.Y]])) ## position of the predictions required by the user in the specific dataset
            iTime.pred <- index.clusterTime[[iC]][iIndex.pred] ## timepoints of the predictions required by the user in the specific dataset
            iPos.pred <- index.cluster[[iC]][iIndex.pred]## position of the predictions required by the user in the original dataset
            iX.pred <- X[iPos.pred,,drop=FALSE]

            if(type.prediction == "dynamic"){
                iIndex.contrast <- iIndex.pred
                iPos.contrast <- iPos.pred
            }else if(type.prediction %in% c("change","auc","auc-b")){
                iIndex.contrast <- 1:length(index.cluster[[iC]])
                iPos.contrast <- index.cluster[[iC]]
            }
            iM.contrast <- M.contrast[index.clusterTime[[iC]][iIndex.contrast],,drop=FALSE]
            
            if(type.prediction == "change"){
                if(1 %in% iIndex.pred){
                    iIndex.varcontrast <- iIndex.contrast[-1] ## do not store uncertainty relative to change between baseline and baseline which is precisely 0
                    iPos.varcontrast <- iPos.contrast[-1]
                }else{
                    iIndex.varcontrast <- iIndex.pred
                    iPos.varcontrast <- iPos.pred
                }
            }else if(sum(se)>0){
                iIndex.varcontrast <- iIndex.contrast
                iPos.varcontrast <- iPos.contrast
            }
            if(sum(se)>0){
                ## if not all necessary timepoint will be set to NA at the end
                iM.varcontrast <- M.contrast[index.clusterTime[[iC]][iIndex.varcontrast],index.clusterTime[[iC]],drop=FALSE]
            }
            
            if(length(iPos.con)>0 && length(iPos.pred)==0 && type.prediction %in% c("change","auc","auc-b")){ ## no prediction but contrast observations
                ##-- compute predictions
                prediction[iPos.contrast] <- iM.contrast %*% iNewdata[[name.Y]]

            }else if(length(iPos.pred)>0 && length(iPos.con)==0){ ## static prediction

                ##-- compute predictions (all observations)
                ## Here we use that iIndex.pred = index.cluster[[iC]] since length(iIndex.con)=0
                ## so iX = iX.pred and prediction[iPos.pred] is prediction[index.cluster[[iC]]]
                iX[index.clusterTime[[iC]],] <- iX.pred
                iX[-index.clusterTime[[iC]],] <- 0  ## make sure no confusion between individuals
                prediction[iPos.pred] <- iM.contrast %*% iX %*% mu
                
                ##-- compute uncertainty about the predictions (all observations except for the baseline when evaluating change - no uncertainty)                
                if(se[1]){
                    grad.var[iPos.varcontrast,name.mu] <- iM.varcontrast %*% iX.pred
                }

                if(se[2]){
                    prediction.var[iPos.varcontrast] <- rowSums((iM.varcontrast %*% iOmega)*iM.varcontrast)
                    if(export.vcov){
                        attr(prediction.var,"vcov")[iPos.varcontrast,iPos.varcontrast] <- iM.varcontrast %*% iOmega %*% t(iM.varcontrast)
                    }
                }else if(se[1] || type.prediction=="change"){
                    prediction.var[iPos.varcontrast] <- 0
                }
                
            }else if(length(iPos.pred)>0){ ## dynamic prediction

                iY[index.clusterTime[[iC]]] <- iNewdata[[name.Y]]
                iY[-index.clusterTime[[iC]]] <- 0 ## make sure no confusion between individuals

                iY.con <- iY[iIndex.con]
                iOmegaM1.con <- solve(iOmega[iIndex.con,iIndex.con,drop=FALSE])
                iResiduals.norm <- iOmegaM1.con %*% (iY.con - iX.con %*% mu)
                iOmega.predcon <- iOmega[iIndex.pred,iIndex.con,drop=FALSE]

                ##-- compute prediction
                iY[iTime.pred] <-  iX.pred %*% mu + iOmega.predcon %*% iResiduals.norm
                prediction[iPos.contrast] <- iM.contrast %*% iY

                ##-- compute uncertainty about the predictions
                if(sum(se)>0){
                    iOmega.OmegaM1 <- iOmega.predcon %*% iOmegaM1.con
                }
                
                if(se[1]>0){
                    idOmega <- dOmega[[newdesign$vcov$pattern[iC]]]
                
                    iGrad <- matrix(0, nrow = length(index.cluster[[iC]]), ncol = n.theta,
                                    dimnames = list(NULL, name.theta))
                    iGrad[iIndex.pred,name.mu] <- iX.pred - iOmega.OmegaM1 %*% iX.con
                    for(iParam in names(idOmega)){ ## iParam <- name.theta[index.vcov][1]
                        ## NOTE: may not contain all vcov parameters when only a subset of the timepoints is passed as dataset
                        iGrad[iIndex.pred,iParam] <- (idOmega[[iParam]][iIndex.pred,iIndex.con,drop=FALSE] - iOmega.OmegaM1 %*% idOmega[[iParam]][iIndex.con,iIndex.con,drop=FALSE]) %*% iResiduals.norm
                    }
                    grad.var[iPos.varcontrast,] <- iM.varcontrast %*% iGrad
                }
                if(se[2]){
                    iOmega.conditional <- iOmega[iIndex.pred, iIndex.pred, drop = FALSE] - iOmega.OmegaM1 %*% t(iOmega.predcon)
                    prediction.var[iPos.varcontrast] <- rowSums((iM.varcontrast[,iIndex.pred,drop=FALSE] %*% iOmega.conditional)*iM.varcontrast[,iIndex.pred,drop=FALSE])
                    if(export.vcov){
                        attr(prediction.var,"vcov")[iPos.varcontrast,iPos.varcontrast] <- iM.varcontrast[,iIndex.pred,drop=FALSE] %*% iOmega.conditional %*% t(iM.varcontrast[,iIndex.pred,drop=FALSE])
                    }
                }else if(se[1] || type.prediction=="change"){
                    prediction.var[iPos.varcontrast] <- 0
                }
            }

            if(type.prediction %in% c("change","auc","auc-b") && (length(iTime.pred) + length(iTime.con) != n.time)){
                iTime.missing <- setdiff(1:n.time, c(iTime.pred,iTime.con))
                test.NAcontrast <- rowSums(abs(iM.contrast[,iTime.missing,drop=FALSE]))>1e-12
                if(any(test.NAcontrast)){ ## add NA if not enough information is provided,
                    ## e.g. only rows for second and third follow-up but no baseline row in the dataset when computing the change
                    prediction[iPos.contrast[test.NAcontrast]] <- NA
                    prediction.var[iPos.contrast[test.NAcontrast]] <- NA
                    grad.var[iPos.contrast[test.NAcontrast],] <- NA
                }
            }
        }

        if(se[1]){
            prediction.var <- prediction.var + rowSums((grad.var %*% vcov.theta) * grad.var)
            if(export.vcov){
                attr(prediction.var,"vcov") <- attr(prediction.var,"vcov") + grad.var %*% vcov.theta %*% t(grad.var)
            }
        }else if(sum(se)==0 & type.prediction == "change"){
            grad.var <- !is.na(prediction.var) ## imputed changes ['hidden' output]
            prediction.var[] <- NA
        }
            
    }

    ## ** df
    if(df){
        prediction.df <- pmax(.dfX(X.beta = grad.var, vcov.param = vcov.theta, dVcov.param = dVcov.theta), options$min.df)
        prediction.df[is.na(prediction) | is.na(prediction.var)] <- NA        
    }else{
        prediction.df <- rep(Inf, length(prediction))
    }

    ## ** restaure NA
    prediction <- restaureNA(unname(prediction), index.na = index.na,
                             level = "obs", cluster = object$cluster)

    if(format == "long" && sum(se)>0){
        prediction.se <- sqrt(restaureNA(unname(prediction.var), index.na = index.na,
                                         level = "obs", cluster = object$cluster))
        prediction.df <- restaureNA(unname(prediction.df), index.na = index.na,
                                    level = "obs", cluster = object$cluster)
    
        alpha <- 1-level
        M.pred <- cbind(estimate = prediction, se = prediction.se, df = prediction.df,
                        lower = prediction + stats::qt(alpha/2, df = prediction.df) * prediction.se,
                        upper = prediction + stats::qt(1-alpha/2, df = prediction.df) * prediction.se)
    }else{
        M.pred <- cbind(estimate = prediction)
    }

    ## ** export
    if(is.null(newdata) && keep.data){
        ## even when no NA, use the initial dataset instead of the augmented one
        newdata <- object$data.original
    }
    out <- .reformat(M.pred, name = names(format), format = format, simplify = simplify,
                     keep.data = keep.data, data = newdata, index.na = index.na,
                     object.cluster = object$cluster, index.cluster = newdata.index.cluster,
                     object.time = object$time, index.time = newdata.index.time,                     
                     call = mycall)

    if(simplify==FALSE){
        attr(out,"grad") <- grad.var
        attr(out,"args") <- list(se = se, df = df, type = type.prediction, level = level, keep.data = keep.data, format = format, simplify = simplify)
    }
    if(export.vcov){
        attr(out,"vcov") <- attr(prediction.var,"vcov")
    }

    return(out)
}

## * predict.mlmm (code)
##' @export
predict.mlmm <- function(object, ...){
    stop("No \'predict\' method for mlmm objects, consider using lmm instead of mlmm. \n")
}

## * .dfX
.dfX <- function(X.beta, vcov.param, dVcov.param, return.vcov = FALSE){

    if(!is.matrix(X.beta) && is.vector(X.beta)){
        X.beta <- rbind(X.beta)
    }
    n.obs <- NROW(X.beta)
    name.beta <- colnames(X.beta)
    n.beta <- length(name.beta)

    vcov.beta <- vcov.param[name.beta,name.beta,drop=FALSE]
    n.param <- NCOL(vcov.param)
    name.param <- colnames(vcov.param)
    prediction.vcov <- X.beta %*% vcov.beta %*% t(X.beta)
    
    ## ** variance covariance matrix of the variance of the mean parameters
    ## delta method:
    ## Cov[Sigma_\beta,Sigma_\beta'] = [dSigma_\beta/d\beta] [Sigma_\theta] [dSigma_\beta'/d\beta']
    n.beta2 <- n.beta^2
    Mname.beta2 <- expand.grid(name.beta,name.beta)
    name.beta2 <- nlme::collapse(Mname.beta2, as.factor = TRUE)
    
    Mpair_dVcov.param <- do.call(rbind,lapply(1:n.beta2, function(iIndex){dVcov.param[Mname.beta2[iIndex,1],Mname.beta2[iIndex,2],]})) ## iIndex <- 240
    vcov.vcovbeta <- Mpair_dVcov.param %*% vcov.param %*% t(Mpair_dVcov.param)
    rownames(vcov.vcovbeta)  <-  name.beta2
    colnames(vcov.vcovbeta)  <-  name.beta2
    
    ## GS <- matrix(0, nrow = n.beta2, ncol = n.beta2, dimnames = list(name.beta2,name.beta2))
    ## for(iP in 1:n.beta2){ ## iP <- 15
    ##     for(iiP in 1:iP){ ## iiP <- 29
    ##         iParam <-  name.beta2[iP]
    ##         iiParam <-  name.beta2[iiP]
    ##         GS[iParam,iiParam] <- dVcov.param[Mname.beta2[iP,1],Mname.beta2[iP,2],] %*% vcov.param %*% dVcov.param[Mname.beta2[iiP,1],Mname.beta2[iiP,2],]
    ##         if(iParam != iiParam){
    ##             GS[iiParam,iParam] <- vcov.vcovbeta[iParam,iiParam]
    ##         }
    ##     }
    ## }
    ## GS[iP,iiP] - vcov.vcovbeta[iP,iiP]

    ## ** gradient of the variance of the predictions relative to the variance of the mean parameters
    ## d(X\Sigma t(X))_ij is computed by X\delta_ij t(X)
    dXTX.dT <- matrix(NA, nrow = n.obs, ncol = n.beta2, dimnames = list(NULL,name.beta2))

    ## matrix cookbook: A_ki B_lj = (A \delta_ij \trans{B})_kl
    ## so               A_ki A_kj = (A \delta_ij \trans{A})_kk
    ## so               t(A)_ik A_kj = (A \delta_ij \trans{A})_kk
    
    for(iP in 1:n.beta2){ ## iP <- 35
        ## iDelta <- matrix(0, nrow = n.beta, ncol = n.beta, dimnames = list(name.beta,name.beta))
        ## iDelta[Mname.beta2[iP,1],Mname.beta2[iP,2]] <- 1
        ## dXTX.dT[,iP] <- diag(X.beta %*% iDelta %*% t(X.beta))
        dXTX.dT[,iP] <- (X.beta[,Mname.beta2[iP,1]] * X.beta[,Mname.beta2[iP,2]])
    }

    ## ** Satterthwaite approximation of the degrees of freedom
    ## NOTE: df(beta) is 2*diag(vcov.beta)^2 / diag(vcov.vcovbeta[Mname.beta2[,1]==Mname.beta2[,2],Mname.beta2[,1]==Mname.beta2[,2]])
    denum <- rowSums((dXTX.dT %*% vcov.vcovbeta) * dXTX.dT)
    out <- 2*diag(prediction.vcov)^2 / denum
    out[denum==0] <- Inf
    ## out <- sapply(1:n.obs, function(iObs){
    ##     2 * prediction.vcov[iObs,iObs]^2 / (dXTX.dT[iObs,] %*% vcov.vcovbeta %*% dXTX.dT[iObs,])
    ## })

    ## ** export
    if(return.vcov){
        attr(out,"vcov") <- vcov.vcovbeta
        attr(out,"contrast") <- X.beta
    }
    return(out)
}


##----------------------------------------------------------------------
### predict.R ends here
