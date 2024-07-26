### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: jul 26 2024 (18:39) 
##           By: Brice Ozenne
##     Update #: 1070
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
##' or the random effects (\code{"ranef"}) when using a random effect model.
##' @param p [numeric vector] value of the model coefficients to be used. Only relevant if differs from the fitted values.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see detail##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param simplify [logical] Should the type of each parameter (mu/sigma/k/rho) and the corresponding variance parameters (sigma/k.x/k.y) be output as attributes?
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
##' @return A vector with the value of the model coefficients. \cr
##' When using \code{simplify=FALSE} the character strings in attribute \code{"type"} refer to:
##' \itemize{
##' \item \code{"mu"}: mean parameters.
##' \item \code{"sigma"}: standard deviation parameters,
##' \item \code{"k"}: ratio between standard deviation,
##' \item \code{"rho"}: correlation parameter
##' }
##' The character strings in attribute \code{"sigma"} refer, for each parameter, to a possible corresponding standard deviation parameter.
##' Those in attribute \code{"k.x"} and \code{"k.y"} refer to the ratio parameter.
##' \code{NA} indicates no corresponding standard deviation or ratio parameter.
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
##' coef(eUN.lmm, effects = "variance", transform.k = "sd")
##' coef(eUN.lmm, effects = "all", simplify = FALSE)

## * coef.lmm (code)
##' @export
coef.lmm <- function(object, effects = NULL, p = NULL,
                     transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE,
                     simplify = TRUE, ...){

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
    ## *** dots
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    
    ## *** effects
    if(is.null(effects)){
        if((is.null(transform.sigma) || identical(transform.sigma,"none")) && (is.null(transform.k) || identical(transform.k,"none")) && (is.null(transform.rho) || identical(transform.rho,"none"))){
            effects <- options$effects
        }else{
            effects <- c("mean","variance","correlation")
        }
    }else if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }else if("ranef" %in% effects){
        if(length(effects)>1){
            stop("Argument \'effects\' should be of length 1 when it contains \"ranef\". \n")
        }
        return(nlme::ranef(object, p = p))
    }else{
        effects <- match.arg(effects, c("mean","fixed","variance","correlation","ranef"), several.ok = TRUE)
        effects[effects== "fixed"] <- "mean"
    }

    ## *** transformation & initialize parameter value
    init <- .init_transform(p = p, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho,
                            table.param = object$design$param)
    test.notransform <- init$test.notransform
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
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

    ## ** apply transformation request by the user
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
    if(!simplify){
        attr(out,"type") <- unname(param.type[param.type %in% keep.type])
        attr(out,"sigma") <- unname(param.sigma[param.type %in% keep.type])
        attr(out,"k.x") <- unname(param.k.x[param.type %in% keep.type])
        attr(out,"k.y") <- unname(param.k.y[param.type %in% keep.type])
    }    
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
##' whereas \code{"ls.contrast"} will extract the contrast matrix applied to each model as a list.##' 
##' @param method [character] type of adjustment for multiple comparisons, one of \code{"none"}, \code{"bonferroni"}, ..., \code{"fdr"}, \code{"single-step"}, \code{"single-step2"}.
##' Alternatively, a method for combining the estimates, one of \code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.gls1"}, \code{"pool.rubin"}.
##' Ignored if \code{type="contrast"} or \code{type="ls.contrast"}.
##' @param ordering [character] should the output be ordered by type of parameter (\code{"parameter"}) or by model (\code{"by"}).
##' Not relevant when \code{effects="contrast"}.
##' @param backtransform [logical] should the estimate, standard error, and confidence interval be back-transformed?
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}.
##' Ignored if \code{type="contrast"}, \code{type="ls.contrast"}, or \code{effects="contrast"}.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"}.
##' Ignored if \code{type="contrast"}, \code{type="ls.contrast"}, or \code{effects="contrast"}.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"}.
##' Ignored if \code{type="contrast"}, \code{type="ls.contrast"}, or \code{effects="contrast"}.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' Ignored if \code{type="contrast"}, \code{type="ls.contrast"}, or \code{effects="contrast"}.
##' @param ... Not used. For compatibility with the generic method.

## * coef.mlmm (code)
##' @export
coef.mlmm <- function(object, effects = "contrast", p = NULL, type = "coef", method = "none", ordering = "by", backtransform = object$args$backtransform,
                      transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE, ...){

    options <- LMMstar.options()
    pool.method <- options$pool.method
    adj.method <- options$adj.method

    ## ** normalize user input

    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character vector")
    }
    valid.effects <- c("contrast","mean","fixed","variance","correlation","all")
    if(any(effects %in% valid.effects == FALSE)){
        if(all(effects == "ls.contrast")){
            stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
                 "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n",
                 "Consider to pass the value \"ls.effects\" to the argument \'type\'. \n")
        }else{
            stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
                 "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
        }
    }    
    if("effects" %in% effects && length(effects)>1){
        stop("Argument \'effects\' must have length 1 when containing the element \'effects\'. \n")
    }
    if("all" %in% effects && length(effects)>1){
        stop("Argument \'effects\' must have length 1 when containing the element \'all\'. \n")
    }
    
    ## *** p
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
            stop("Argument \'type\' must be equal to \"coef\" when the length of argument \"method\" is strictly greater than one. \n")
        }
        if(any(effects != "contrast")){
            stop("Argument \'effects\' must be equal to \"coef\" when the length of argument \"method\" is strictly greater than one. \n")
        }
        if(sum(method %in% adj.method)>1){
            stop("Incorrect argument \'method\': cannot handle several methods to adjust for multiple comparisons. \n")
        }
        if("p.rejection" %in% method & is.null(attr(method,"method")) & sum(method %in% adj.method)==1){
            attr(method,"method") <- intersect(method, adj.method)
        }
    }else{
        name.method <- names(method)
    }

    ## *** ordering
    ordering <- match.arg(ordering, c("by","parameter"))

    ## ** extract
    if(type == "contrast"){
        out <- object$glht$all[[1]]$linfct
    }else if(type == "ls.contrast"){
        out <- object$glht$all[[1]]$linfct.original
        names(out) <- names(object$model)
    }else if(all(effects=="contrast")){ ## contrast of mean/variance/correlation parameters
        
        ls.out <- lapply(method, function(iMethod){ ## iMethod <- method[1]
    
            if(iMethod %in% pool.method){

                if(method == "p.rejection"){
                    iOut <- proportion.mlmm(object = object, p = p, 
                                            name.method = name.method, method = attr(method,"method"), qt = attr(method,"qt"),
                                            null = NA, ci = FALSE, df = FALSE, alpha = NA)
                }else{
                    iOut <- poolWald.mlmm(object = object, p = p, 
                                          method = method, name.method = name.method,
                                          ci = FALSE, df = FALSE, alpha = NA)
                }
            
                iOut <- stats::setNames(table.out[,"estimate"],rownames(table.out))
            
            }else{ 
                if(is.null(p)){
                    iOut <- coef.Wald_lmm(object, type = "coef", backtransform = backtransform)
                }else{
                    contrast <- coef(object, type = "contrast") ## extract linear combination from each model
                    ls.contrast <- coef(object, type = "ls.contrast") ## contrast combinations across models
                    iCoef <- coef(object, p = p, effects = "all", ordering = "by",
                                  transform.sigma = object$args$transform.sigma, transform.k = object$args$transform.k,
                                  transform.rho = object$args$transform.rho, transform.names = FALSE)
                    iCCoef <- do.call(rbind,mapply(iC = ls.contrast, iTheta = iCoef, function(iC,iTheta){iC %*% iTheta}, SIMPLIFY = FALSE))
                    iOut <- (contrast %*% iCCoef)[,1]
                    
                    if(backtransform && object$args$backtransform){
                        iOut[] <- .backtransform(data.frame(estimate = iOut), type.param = object$univariate$type,  
                                                               backtransform = TRUE, backtransform.names = object$args$backtransform.names[[1]],
                                                               transform.mu = "none",
                                                               transform.sigma = object$args$transform.sigma,
                                                               transform.k = object$args$transform.k,
                                                               transform.rho = object$args$transform.rho)$estimate
                    }
                }

                if(ordering=="by"){
                    ## nothing to do
                }else if(ordering=="parameter"){
                    if(all(lengths(object$univariate$parameter)==1)){ ## single variable
                        iLevel <- unique(unlist(object$univariate$parameter)) ## keep ordering univariate table that should remember how factors are ordered.
                        iOut <- iOut[order(factor(unlist(object$univariate$parameter), levels = iLevel))]
                    }else{ ## multiple variables, e.g. when a contrast is based on two model parameters X1 + 2*X2
                        iLevel <- unique(sapply(object$univariate$parameter, paste, sep = "XX_XX"))
                        iOut <- iOut[order(factor(sapply(object$univariate$parameter, paste, sep = "XX_XX"),iLevel))]
                    }
                }
             
            }
            return(iOut)
        })

        if(length(method)==1){
            out <- ls.out[[1]]
        }else{
            out <- do.call(rbind, ls.out)
        }

    }else{ ## mean/variance/correlation/all parameters
        if(is.null(p)){
            out <- lapply(object$model, coef, effects = effects,
                          transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        }else{
            out <- lapply(names(object$model), function(iBy){
                coef(object$model[[iBy]], p = p[[iBy]], effects = effects,
                     transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
            })
        }

        ## only reorder when no pooling
        if(ordering=="by"){
            ## nothing to do
        }else if(ordering=="parameter"){
            all.parameters <- unique(unlist(lapply(out,names)))
            out <- stats::setNames(lapply(all.parameters, function(iParam){
                sapply(out,"[[",iParam)
            }), all.parameters)
        }
        
    }

    ## ** export
    return(out)
}

## * coef.rbindWald_lmm
##' @title Extract Coefficients From Wald Tests for Linear Mixed Model
##' @description Extract coefficients from Wald tests applied to a linear mixed model.
##'
##' @param object a \code{Wald_lmm} object.
##' @param method [character vector] should the linear contrasts involved in the Wald test (\code{"none"} or any adjustment for multiple comparisons),
##' or a combination of the contrasts (\code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.gls1"}, \code{"pool.rubin"}, \code{"p.rejection"}) be output?
##' Can also out the value of the linear mixed model parameters (\code{"all"}: after transformation, \code{"all.original"}: without transformation),
##' or the linear mixed model type of parameters (\code{"all.type"}) 
##' @param ordering [character] should the output be ordered by name of the linear contrast (\code{"contrast"}) or by model (\code{"model"}).
##' Only relevant when \code{effects="contrast"}, \code{effects="all"}, or \code{effects="all.original"}.
##' @param backtransform [logical] should the estimate, standard error, and confidence interval be back-transformed?
##' @param simplify [logical] should the output be a vector or a list with one element specific to each possible ordering (i.e. contrast or model).
##' Not relevant when argument \code{effects} refers to combination of contrasts.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details Argument \bold{effects}: when evaluating the proportion of rejected hypotheses (\code{effects="p.rejection"})
##' a \code{"single-step"} method will be used by default to evaluate the critical quantile.
##' This can be changed by adding adjustment method, e.g. \code{effects=c("bonferronin","p.rejection"}, in the argument.
##' 
##' @export
coef.rbindWald_lmm <- function(object, method = "none", ordering = NULL, backtransform = NULL, simplify = TRUE, ...){

    options <- LMMstar.options()
    pool.method <- options$pool.method
    adj.method <- options$adj.method

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** object
    if(object$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling rbind.Wald_lmm. \n")
        return(invisible(NULL))
    }

    ## *** method
    if(!is.character(method) || !is.vector(method)){
        stop("Argument \'method\' must be a character.")
    }
    valid.method <- c("none","all","all.original","all.type",pool.method,adj.method)
    if(any(method %in% valid.method == FALSE)){
        stop("Incorrect value for argument \'method\': \"",paste(setdiff(method,valid.method), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.method, collapse ="\", \""),"\". \n")
    }
    if(length(method)>1 && any(c("all","all.original","all.type") %in% method)){
        stop("Argument \'method\' must have length 1 when \"all\", \"all.original\", or \"all.type\" is in \"method\". \n")
    }
    if(sum(method %in% c("none",adj.method))>1){
        stop("Argument \'method\' must refer no more than one adjustment for multiple comparisons. \n",
             "Proposed adjustments: \"",paste(intersect(method, c("none",adj.method)), collapse ="\", \""),"\". \n")
    }

    ## *** ordering
    if(!is.null(ordering)){
        ordering <- match.arg(ordering, c("contrast","model"))
        if(any(method %in% pool.method)){
            message("Argument \'ordering\' is ignored when argument \'method\' is \"",intersect(pool.method,method)[1],"\". \n")
            ordering <- NULL
        }
    }
    
    ## *** backtransform
    if(is.null(backtransform)){
        backtransform <- any(object$univariate$tobacktransform)
    }else if(is.character(backtransform)){
        backtransform <-  eval(parse(text=backtransform))
    }else if(is.numeric(backtransform)){
        backtransform <- as.logical(backtransform)
    }

    ## *** simplify
    if(any(method %in% pool.method)){
        simplify <- TRUE
    }

    ## ** extract from object
    if(any(method %in% c("all","all.original","all.type"))){

        if(method == "all"){
            value.out <- stats::setNames(object$glht[[1]]$coef,object$glht[[1]]$coef.name)
        }else if(method == "all.original"){
            value.out <- object$glht[[1]]$coef.notrans
        }else if(method == "all.type"){
            value.out <- object$glht[[1]]$coef.type
        }

        if(!is.null(ordering)){
            ordering.out <- object$glht[[1]][[c(model = "model", contrast = "term")[ordering]]]
        }else if(simplify == FALSE){
            attr(value.out,"parameter") <- object$glht[[1]]$term
            attr(value.out,"model") <- object$glht[[1]]$model
        }

    }else{

        if(any(method %in% c("none",adj.method))){
            table.univariate <- object$univariate

            if(is.function(backtransform) || identical(backtransform,TRUE)){

                if(is.function(backtransform)){

                    df.out <- .backtransform(table.univariate[,"estimate",drop=FALSE], type.param = table.univariate$type,
                                             backtransform = TRUE, backtransform.names = NULL,
                                             transform.mu = backtransform,
                                             transform.sigma = backtransform,
                                             transform.k = backtransform,
                                             transform.rho = backtransform)

            }else{

                ## force no back-transform, e.g. when comparing two correlation coefficients
                df.out <- .backtransform(table.univariate[,"estimate",drop=FALSE], type.param = ifelse(table.univariate$tobacktransform,table.univariate$type,"mu"),  
                                         backtransform = TRUE, backtransform.names = NULL,
                                         transform.mu = "none",
                                         transform.sigma = object$args$transform.sigma,
                                         transform.k = object$args$transform.k,
                                         transform.rho = object$args$transform.rho)
            
            }
                value.out <- stats::setNames(df.out$estimate,rownames(df.out))

            }else{
                df.out <- table.univariate[,"estimate",drop = FALSE]
                value.out <- stats::setNames(df.out$estimate,rownames(df.out))
            }

            if(!is.null(ordering)){
                ordering.out <- table.univariate[[c(model = "model", contrast = "term")[ordering]]]
            }else if(simplify == FALSE){
                attr(value.out,"parameter") <- table.univariate$term
                attr(value.out,"model") <- table.univariate$model
            }
        }else{
            value.out <- NULL
        }

        if("p.rejection" %in% method){
            value.out <- c(value.out,
                           proportion.mlmm(object = object, p = p, 
                                           name.method = name.method, method = attr(method,"method"), qt = attr(method,"qt"),
                                           null = NA, ci = FALSE, df = FALSE, alpha = NA)
                           )
        }
        if(method %in% setdiff(pool.method,"p.rejection")){
            value.out <- c(value.out,
                           poolWald.mlmm(object = object, p = p, 
                                         method = method, name.method = name.method,
                                         ci = FALSE, df = FALSE, alpha = NA)
                           )
        }
    }
    
    ## ** export
    if(simplify == FALSE && !is.null(ordering)){
        out <- tapply(value.out, INDEX = ordering.out, FUN = base::identity, simplify = FALSE)[unique(ordering.out)]
    }else if(!is.null(ordering)){
        out <- value.out[order(ordering.out)]
    }else{
        out <- value.out
    }
    return(out)
}

## * coef.Wald_lmm
##' @title Extract Coefficients From Wald Tests for Linear Mixed Model
##' @description Extract coefficients from Wald tests applied to a linear mixed model.
##'
##' @param object a \code{Wald_lmm} object.
##' @param method [character] should the linear contrasts involved in the Wald test (\code{"none"}),
##' the value of the linear mixed model parameters (\code{"all"}: after transformation, \code{"all.original"}: without transformation),
##' or the linear mixed model type of parameters (\code{"all.type"}) be output?
##' @param backtransform [logical] should the estimate, standard error, and confidence interval be back-transformed?
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @export
coef.Wald_lmm <- function(object, method = "none", backtransform = NULL, ...){

    options <- LMMstar.options()
    adj.method <- options$adj.method
    if(object$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling anova. \n")
        return(invisible(NULL))
    }

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** method
    if(!is.character(method) || !is.vector(method)){
        stop("Argument \'method\' must be a character.")
    }
    if(length(method)!=1){
        stop("Argument \'method\' must have length 1.")
    }
    if(object$args$type=="auto"){
        valid.method <- c("none")
    }else{
        valid.method <- c("none","all","all.type","all.original",adj.method)
    }
    if(any(method %in% valid.method == FALSE)){
        stop("Incorrect value for argument \'method\': \"",paste(setdiff(method,valid.method), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.method, collapse ="\", \""),"\". \n")
    }    
    if(method %in% adj.method){
        method <- "none"
    }

    ## *** backtransform
    if(is.null(backtransform)){
        backtransform <- any(object$univariate$tobacktransform)
    }else if(is.character(backtransform)){
        backtransform <-  eval(parse(text=backtransform))
    }else if(is.numeric(backtransform)){
        backtransform <- as.logical(backtransform)
    }
    
    ## ** extract from object
    if(method == "all"){
        out <- stats::setNames(object$glht[[1]]$coef,object$glht[[1]]$coef.name)
    }else if(method == "all.type"){
        out <- object$glht[[1]]$coef.type
    }else if(method == "all.original"){
        out <- object$glht[[1]]$coef.notrans
    }else if(method == "none"){
        table.univariate <- object$univariate
    
        if(is.function(backtransform) || identical(backtransform,TRUE)){

            if(is.function(backtransform)){

                df.out <- .backtransform(table.univariate[,"estimate",drop=FALSE], type.param = table.univariate$type,
                                         backtransform = TRUE, backtransform.names = NULL,
                                         transform.mu = backtransform,
                                         transform.sigma = backtransform,
                                         transform.k = backtransform,
                                         transform.rho = backtransform)

            }else{

                ## force no back-transform, e.g. when comparing two correlation coefficients
                df.out <- .backtransform(table.univariate[,"estimate",drop=FALSE], type.param = ifelse(table.univariate$tobacktransform,table.univariate$type,"mu"),  
                                         backtransform = TRUE, backtransform.names = NULL,
                                         transform.mu = "none",
                                         transform.sigma = object$args$transform.sigma,
                                         transform.k = object$args$transform.k,
                                         transform.rho = object$args$transform.rho)
            
            }
            out <- stats::setNames(df.out$estimate,rownames(df.out))
        }else{
            df.out <- table.univariate[,"estimate",drop = FALSE]
            out <- stats::setNames(df.out$estimate,rownames(df.out))
        }
    }

    ## ** export
    return(out)
}



##----------------------------------------------------------------------
### coef.R ends here
