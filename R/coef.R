### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: jul 11 2024 (16:36) 
##           By: Brice Ozenne
##     Update #: 828
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
##'
##' @param object a \code{lmm} object.
##' @param effects [character] Should all coefficients be output (\code{"all"}),
##' or only coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only coefficients relative to the variance structure (\code{"variance"}),
##' or only coefficients relative to the correlation structure (\code{"correlation"}).
##' Can also be \code{"ranef"} to output random effect (only for \code{CS} structure).
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
##' @seealso
##' \code{\link{confint.lmm}} or \code{\link{model.tables.lmm}} for a data.frame containing estimates with their uncertainty. \cr
##' 
##' @keywords methods
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
##' @export
coef.lmm <- function(object, effects = NULL, p = NULL,
                     transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE, ...){

    mycall <- match.call()
    
    ## ** extract from object
    param.name <- object$design$param$name
    param.type <- stats::setNames(object$design$param$type,param.name)
    param.level <- stats::setNames(object$design$param$level,param.name)
    param.sigma <- stats::setNames(object$design$param$sigma,param.name)
    param.k.x <- stats::setNames(object$design$param$k.x,param.name)
    param.k.y <- stats::setNames(object$design$param$k.y,param.name)

    object.reparametrize.name <- names(object$reparametrize$p)
    object.reparametrize.value <- object$reparametrize$p
    object.reparametrize.newname <- object$reparametrize$newname

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
    if("ranef" %in% effects){
        if(length(effects)>1){
            stop("Argument \'effects\' should be of length 1 when it contains \"ranef\". \n")
        }
        return(ranef(object, p = p))
    }
    effects <- match.arg(effects, c("mean","fixed","variance","correlation","ranef"), several.ok = TRUE)
    effects[effects== "fixed"] <- "mean"

    ## initialize parameter values
    init <- .init_transform(p = p, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho,
                            table.param = object$design$param)
    test.notransform <- init$test.notransform
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho

    if(transform.rho=="cov"){
        if("transform.k" %in% names(mycall) == FALSE){
            transform.k <- "var"
        }
        if("transform.sigma" %in% names(mycall) == FALSE){
            transform.sigma <- "square"
        }
    }
    if(transform.k %in% c("logsd","var","logvar") && "transform.sigma" %in% names(mycall) == FALSE){
        transform.sigma <- switch(transform.k,
                                  "logsd" = "log",
                                  "var" = "square",
                                  "logvar" = "logsquare")
                                      
    }
    transform <- init$transform
    if(!is.null(p)){
        theta <- init$p
    }else{
        theta <- object$param
    }

    effects2 <- effects
    if(transform.rho == "cov"){
        if(all("correlation" %in% effects == FALSE)){
            stop("Cannot use the argument \'transform.rho\' set to \"cov\" when \"correlation\" is not in argument \'effect\'. \n")
        }
        if(all("variance" %in% effects == FALSE)){
            effects2 <- c("variance",effects2)
        }
    }

    ## apply transformation request by the user
    if(is.null(p) && test.notransform){
        theta.trans <- theta
        theta.trans[match(object.reparametrize.name, names(theta))] <- object.reparametrize.value
        if(transform.names){
            names(theta.trans)[match(object.reparametrize.name, names(theta))] <- object.reparametrize.newname
        }        
    }else if((transform.sigma == "none" || "variance" %in% effects2 == FALSE) && (transform.k == "none" || "variance" %in% effects2 == FALSE) && (transform.rho == "none" || "correlation" %in% effects2 == FALSE)){
        theta.trans <- theta
    }else{
        reparam <- .reparametrize(p = theta[object.reparametrize.name],  
                                  type = param.type[object.reparametrize.name],
                                  sigma = param.sigma[object.reparametrize.name],
                                  k.x = param.k.x[object.reparametrize.name],
                                  k.y = param.k.y[object.reparametrize.name],
                                  level = param.level[object.reparametrize.name],                                              
                                  Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                  transform.sigma = transform.sigma,
                                  transform.k = transform.k,
                                  transform.rho = transform.rho,
                                  transform.names = transform.names)
        theta.trans <- theta
        theta.trans[match(object.reparametrize.name, names(theta))] <- reparam$p
        if(transform.names){
            names(theta.trans)[match(object.reparametrize.name, names(theta))] <- reparam$newname
        }        
    }

    ## ** extract
    keep.type <- unlist(lapply(effects, switch,
                               "mean" = "mu",
                               "variance" = c("sigma","k"),
                               "correlation" = "rho"))
    out <- theta.trans[param.type %in% keep.type]

    ## ** export
    return(out)
}

## * coef.lmmCC (code)
##' @export
coef.lmmCC <- function(object, effects = NULL, ...){

    if(object$time$n==4 && (is.null(effects) || effects == "change")){

        dots <- list(...)
        if(length(dots)>0){
            stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
        }
    
        Mcon <- cbind(c(-1,1,0,0),c(0,0,-1,1))
        Sigma.change <- t(Mcon) %*% stats::sigma(object) %*% Mcon
        out <- c(cor = stats::cov2cor(Sigma.change)[1,2],
                 beta = Sigma.change[1,2]/Sigma.change[1,1])
        
    }else{

        class(object) <- setdiff(class(object),"lmmCC")
        out <- coef(object, effects = effects, ...)

    }

    ## ** export
    return(out)

}

## * coef.LRT_lmm
##' @export
coef.LRT_lmm <- function(object, ...){
    message("No effect size available for likelihood ratio tests.")
    return(NULL)
}

## * coef.mlmm (documentation)
##' @title Extract Coefficients From Multiple Linear Mixed Models
##' @description Extract coefficient or constrast coefficients from multiple linear mixed models.
##'
##' @param object a \code{mlmm} object.
##' @param effects [character] By default will output the estimates relative to the hypotheses being tested (\code{"contrast"}).
##' But can also output all model coefficients (\code{"all"}),
##' or only coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only coefficients relative to the variance structure (\code{"variance"}),
##' or only coefficients relative to the correlation structure (\code{"correlation"}).
##' Not relevant when \code{type="contrast"} or \code{type="ls.contrast"}.
##' @param p [list of numeric vector] list of model coefficients to be used. Only relevant if differs from the fitted values.
##' @param type [character] Should the coefficients be extracted (\code{"coef"}) or the contrast matrix (\code{"contrast"} or \code{"ls.contrast"}).
##' \code{"contrast"} will extract the contrast matrix acrossed models (typically identity when testing an exposure effect in each model separately)
##' whereas \code{"ls.contrast"} will extract the contrast matrix applied to each model as a list.
##' @param ordering [character] should the output be ordered by type of parameter (\code{parameter}) or by model (\code{by}).
##' @param ... passed to \code{coef.Wald_lmm}.

## * coef.mlmm (code)
##' @export
coef.mlmm <- function(object, effects = "contrast", p = NULL, type = "coef", ordering = "parameter", ...){

    ## ** normalize user input

    ## effects
    if(!is.null(effects)){
        effects <- match.arg(effects, c("contrast","mean","fixed","variance","correlation","all"), several.ok = TRUE)
    }

    ## p
    if(!is.null(p)){
        if(!is.list(p)){
            stop("Argument \'p\' should either be NULL or a list. \n")
        }
        if(is.null(names(p))){
            stop("Argument \'p\' should either be NULL or a named list. \n")
        }
        if(any(names(p) %in% names(object$model) == FALSE)){
            stop("Incorrect names for argument \'p\': \"",paste(setdiff(names(p),names(object$model)), collapse = "\", \""),"\". \n", 
                 "Should be among \"",paste(names(object$model), collapse = "\", \""),"\". \n")
        }
        if(!is.null(object$univariate) & all(object$univariate$type=="mu")){
            effects2 <- "mean"
        }else{
            effects2 <- "all"
        }
    }else{
        effects2 <- effects
    }

    ordering <- match.arg(ordering, c("by","parameter"))
    type <- match.arg(type, c("coef","contrast","ls.contrast"))
    

    ## ** extract
    if(type == "contrast"){
        out <- object$glht$all[[1]]$linfct
    }else if(type == "ls.contrast"){
        out <- object$glht$all[[1]]$linfct.original
        names(out) <- names(object$model)
    }else if(!is.null(effects) && length(effects)==1 && effects=="contrast" && is.null(p)){
        out <- coef.Wald_lmm(object, backtransform = object$args$backtransform, ...)
        if(length(out) == NROW(object$univariate)){
            ## only reorder when no pooling
            if(ordering=="by"){
                out <- out[order(object$univariate[["by"]])]
            }else if(is.list(object$univariate$parameter)){
                out <- out[order(object$univariate$type,sapply(object$univariate$parameter, paste, collapse=";"))]
            }else{
                out <- out[order(object$univariate$type,object$univariate$parameter)]
            }
        }
    }else{
        if(is.null(p)){
            ls.out <- lapply(object$model, coef, effects = effects, ...)
        }else{
            ls.out <- mapply(object = object$model, p = p, FUN = coef, effects = effects2, ..., SIMPLIFY = FALSE)
        }
        if(!is.null(effects) && length(effects)==1 && effects=="contrast"){
            out <- mapply(iC = coef(object, type = "ls.contrast"), iMu = ls.out, FUN = function(iC,iMu){iC[,names(iMu),drop=FALSE] %*% cbind(iMu)})
        }else if(ordering == "by"){
            out <- ls.out
        }else if(ordering == "parameter"){
            Uname <- unique(unlist(lapply(ls.out,names)))
            ls.out2 <- stats::setNames(lapply(Uname, function(iName){ ## iName <- "X1"
                unlist(lapply(ls.out,function(iVec){unname(iVec[iName])}))
            }), Uname)
            out <- ls.out2
        }
    }

    ## ** export
    return(out)
}
## * coef.Wald_lmm
##' @title Extract Coefficients From Wald Tests for Linear Mixed Model
##' @description Extract coefficients from Wald tests applied to a linear mixed model.
##'
##' @param object a \code{Wald_lmm} object.
##' @param type [character] Should the coefficients be extracted (\code{"coef"}) or the contrast matrix (\code{"contrast"} or \code{"ls.contrast"}).
##' \code{"contrast"} will try to simplify the output into a matrix whereas \code{"ls.contrast"} will keep the original format (list of list of matrix).
##' @param backtransform [logical] should the estimate, standard error, and confidence interval be back-transformed?
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @export
coef.Wald_lmm <- function(object, type = "coef", method = "none", backtransform = object$args$backtransform, ...){

    options <- LMMstar.options()
    pool.method <- options$pool.method
    adj.method <- options$adj.method
    table.univariate <- object$univariate
    
    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** type
    type <- match.arg(type, c("coef","contrast","ls.contrast"))
    
    ## *** method
    if(any(method %in% c(adj.method,pool.method)==FALSE)){
        stop("Unknown value for argument \'type\': \"",paste(setdiff(method, c(adj.method,pool.method)),collapse = "\", \""),"\". \n",
             "Possible values: \"",paste(c(adj.method,pool.method), collapse = "\", \""),"\". \n")
    }


    if(length(method)>1){
        ##  handle multiple pooling technics
        if(type != "coef"){
            message("Argument \'type\' ignored when argument \"method\" is a vector. \n")
        }
        if(sum(method %in% adj.method)>1){
            stop("Incorrect argument \'method\' \n",
                 "coef.Wald_lmm cannot handle several methods to adjust for multiple comparisons.")
        }
        if("p.rejection" %in% method & is.null(attr(method,"method")) & sum(method %in% adj.method)==1){
            attr(method,"method") <- intersect(method, adj.method)
        }

        ls.coef <- lapply(method, function(iMethod){
            if(iMethod == "p.rejection"){
                attr(iMethod,"method") <- attr(method,"method")
                attr(iMethod,"qt") <- attr(method,"qt")
            }
            iOut <- coef.Wald_lmm(object = object, type = "coef", method = iMethod, backtransform = backtransform, ...)
            return(iOut)
        })
        out <- do.call(c, ls.coef)
        return(out)
    }else{
        name.method <- names(method)
    }
    
    ## ** extract from object
    if(type %in% c("contrast","ls.contrast")){
        out <- lapply(object$glht, function(iGlht){lapply(iGlht,"[[","linfct")})
        if(type == "contrast"){
            if(length(out)==1){
                out <- out[[1]]
                if(length(out)==1){
                    out <- out[[1]]
                }
            }else{
                out <- lapply(out, function(iOut){if(length(iOut)==1){iOut[[1]]}else{iOut}})
            }
        }
    }else if(type == "coef"){

        if(method %in% pool.method){

            grid <- unique(table.univariate[,c("type","test"),drop=FALSE])
            n.grid <- NROW(grid)
            table.out <- do.call(rbind,lapply(1:n.grid , function(iGrid){
                if(n.grid>1){
                    iIndex.table <- intersect(which(table.univariate$type==grid$type[iGrid]),
                                              which(table.univariate$test==grid$test[iGrid]))
                }else{
                    iIndex.table <- 1:NROW(table.univariate)
                }
                if(method == "p.rejection"){
                    iOut <- proportion.mlmm(object = object, index = iIndex.table,
                                            name.method = name.method, method = attr(method,"method"), qt = attr(method,"qt"),
                                            null = NA, ci = FALSE, df = FALSE, alpha = NA)
                }else{
                    iOut <- poolWald.mlmm(object = object, index = iIndex.table,
                                          method = method, name.method = name.method,
                                          ci = FALSE, df = FALSE, alpha = NA)
                }
            }))
            out <- stats::setNames(table.out[,"estimate"],rownames(table.out))
            
        }else{

            if(is.null(table.univariate)){
                out <- NULL
            }else if(!backtransform){
                out <- stats::setNames(table.univariate$estimate, rownames(table.univariate))
            }else{ ## backtransformation
                tableBack.univariate <- .backtransform(table.univariate, type.param = table.univariate$type,  
                                                       backtransform = TRUE, backtransform.names = object$args$backtransform.names[[1]],
                                                       transform.mu = "none",
                                                       transform.sigma = object$args$transform.sigma,
                                                       transform.k = object$args$transform.k,
                                                       transform.rho = object$args$transform.rho)

                vec.backtransform <- attr(table.univariate,"backtransform")
                if(!is.null(vec.backtransform)){
                    ## case where a contrast is performed on transformed coefficients (e.g. sigma:male vs sigma:female)
                    ## the back transformed version exp(log(sigma:male) - log(sigma:female)) differs from the original version sigma:male - sigma:female
                    ## thus without further indication the original version is output
                    tableBack.univariate[names(vec.backtransform),"estimate"] <- unname(vec.backtransform)
                }
                out <- stats::setNames(tableBack.univariate$estimate, rownames(tableBack.univariate))
            }
        }
    }

    ## ** export
    return(out)
}



##----------------------------------------------------------------------
### coef.R ends here
