### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: okt  3 2024 (17:11) 
##           By: Brice Ozenne
##     Update #: 1260
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
##' @description Extract estimated parameters from a linear mixed model.
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
##' @param simplify [logical] Omit from the output the attribute containing the type of each parameter (mu/sigma/k/rho) and the corresponding variance parameters (sigma/k.x/k.y).
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
                     transform.sigma = NULL, transform.k = NULL, transform.rho = "none", transform.names = TRUE,
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
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
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
    }else{
        if(!is.character(effects) || !is.vector(effects)){
            stop("Argument \'effects\' must be a character vector. \n")
        }
        valid.effects <- c("ranef","mean","fixed","variance","correlation","all")
        if(any(effects %in% valid.effects == FALSE)){
            stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
                 "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
        }
        if(all("ranef" %in% effects)){
            if(length(effects)>1){
                stop("Argument \'effects\' must have length 1 when containing the element \"ranef\". \n")
            }
        }else if(all("all" %in% effects)){
            if(length(effects)>1){
                stop("Argument \'effects\' must have length 1 when containing the element \"all\". \n")
            }else{
                effects <- c("mean","variance","correlation")
            }
        }else{
            effects[effects == "fixed"] <- "mean"
        }
    }

    ## *** transformation & initialize parameter value
    if(is.null(transform.k) && (!is.null(transform.rho) && transform.rho == "none")){
        transform.k <- "none"
    }
    if(is.null(transform.sigma) && (!is.null(transform.k) && transform.k == "none") && (!is.null(transform.rho) && transform.rho == "none")){
        transform.sigma <- "none"
    } 
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

    ## *** simplify
    if(!is.numeric(simplify) && !is.logical(simplify)){
        stop("Argument \'simplify\' must be numeric or logical. \n")
    }
    if(length(simplify)!=1){
        stop("Argument \'simplify\' must have length 1. \n")
    }
    if(simplify %in% c(0,1) == FALSE){
        stop("Argument \'simplify\' must be TRUE/1 or FALSE/0. \n")
    }

    ## ** special case
    if(all(effects %in% "ranef")){
        return(nlme::ranef(object, p = p, simplify = simplify))
    }

    ## ** apply transformation request by the user
    if(is.null(p) && test.notransform){
        theta.trans <- theta
        theta.trans[match(object.reparametrize.name, names(theta))] <- object.reparametrize.value
        if(transform.names && !is.null(object.reparametrize.newname)){
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
        if("options" %in% names(dots) && !is.null(dots$options)){
            options <- dots$options
        }else{
            options <- LMMstar.options()
        }
        dots$options <- NULL
        if(length(dots)>0){
            stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
        }
    
        Mcon <- cbind(c(-1,1,0,0),c(0,0,-1,1))
        Sigma.change <- t(Mcon) %*% stats::sigma(object) %*% Mcon
        out <- c(cor = stats::cov2cor(Sigma.change)[1,2],
                 beta = Sigma.change[1,2]/Sigma.change[1,1])
        
    }else{

        class(object) <- setdiff(class(object),"lmmCC")
        out <- stats::coef(object, effects = effects, options = options, ...)

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
##' @description Combine estimated parameters or linear contrasts applied on parameters from group-specific linear mixed models.
##'
##' @param object a \code{mlmm} object.
##' @param effects [character] By default will output the estimates relative to the hypotheses being tested (\code{"Wald"}).
##' But can also output all model coefficients (\code{"all"}),
##' or only coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only coefficients relative to the variance structure (\code{"variance"}),
##' or only coefficients relative to the correlation structure (\code{"correlation"}).
##' @param method [character vector] should the estimated value for the linear contrasts be output (one of \code{"none"}, \code{"bonferroni"}, ..., \code{"single-step2"})
##' and/or pooled linear contrast estimate(s) (\code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.gls1"}, \code{"pool.rubin"}, \code{"p.rejection"})?
##' Only relevant when \code{effects = "Wald"}.
##' @param p [list of numeric vector] values for the model parameters to be used to evaluate the estimates relative to the hypotheses being tested.
##' Only relevant if differs from the fitted values.
##' @param ordering [character] should the output be ordered by name of the linear contrast (\code{"contrast"}) or by model (\code{"model"}).
##' @param backtransform [logical] should the estimate be back-transformed?
##' Only relevant when \code{effects="Wald"}.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}.
##' Ignored when \code{effects="Wald"}.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"}.
##' Ignored when \code{effects="Wald"}.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"}.
##' Ignored when \code{effects="Wald"}.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' Ignored when \code{effects="Wald"}.
##' @param simplify [logical] should the output be a vector or a list with one element specific to each possible ordering (i.e. contrast or model).
##' Only relevant when argument \code{method} refers to multiple comparisons and not to a pooling method.
##' @param ... Not used. For compatibility with the generic method.

## * coef.mlmm (code)
##' @export
coef.mlmm <- function(object, effects = "Wald", method = "none", p = NULL, ordering = "model",
                      backtransform = object$args$backtransform, transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE,
                      simplify = TRUE, ...){


    ## ** normalize user input

    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    pool.method <- options$pool.method
    adj.method <- options$adj.method

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character vector")
    }
    valid.effects <- c("Wald","mean","fixed","variance","correlation","all")
    if(any(effects %in% valid.effects == FALSE)){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }        
    if("Wald" %in% effects && length(effects)>1){
        stop("Argument \'effects\' must have length 1 when containing the element \'Wald\'. \n")
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

    ## *** ordering
    ordering <- match.arg(ordering, c("contrast","model"))

    ## ** extract
    if("Wald" %in% effects){
        message <- NULL
        if(!is.null(transform.sigma) && transform.sigma!="none"){
            message <- c(message,"transform.sigma")
        }
        if(!is.null(transform.k) && transform.k!="none"){
            message <- c(message,"transform.k")
        }
        if(!is.null(transform.rho) && transform.rho!="none"){
            message <- c(message,"transform.rho")
        }
        if(!is.null(transform.names) && transform.names!=TRUE){
            message <- c(message,"transform.names")
        }
        if(length(message)==1){
            message("Argument \'",message,"\' is ignored when argument \"effects\" is \"Wald\". \n")
        }else if(length(message)>1){
            message("Arguments \'",paste(message, collapse = "\', \'"),"\' are ignored when argument \"effects\" is \"Wald\". \n")
        }
        if(!is.null(p)){
            ls.lmm <- object$model
            name.lmm <- names(ls.lmm)
            ls.Cmat <- model.tables(object, effects = "contrast", transform.names = FALSE, simplify = FALSE)
            ls.anova <- stats::setNames(lapply(name.lmm, function(iName){ ## iName <- name.lmm[1]

                iWald <- stats::anova(ls.lmm[[iName]], effects = ls.Cmat[[iName]], rhs = rhs[[iName]], robust = robust, df = df,
                                      univariate = TRUE, multivariate = FALSE,
                                      transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                      options = options)

                return(iWald)
            }), name.lmm)
            out <- do.call("rbind.Wald_lmm",
                   args = c(list(model = ls.anova[[1]], effects = contrast.rbind, rhs = rhs.by, name = name.model, sep = sep), unname(ls.anova[-1]))
                   )
    
        }
        out <- coef.rbindWald_lmm(object, effects = "Wald", method = method, ordering = ordering, backtransform = backtransform, simplify = simplify)
    }else{
        message <- NULL
        if(!is.null(method) && method!="none"){
            message <- c(message,"method")
        }
        if(!is.null(backtransform) && backtransform!=object$args$backtransform){
            message <- c(message,"backtransform")
        }
        if(length(message)==1){
            message("Argument \'",message,"\' is ignored when argument \"effects\" is not \"Wald\". \n")
        }else if(length(message)>1){
            message("Arguments \'",paste(message, collapse = "\', \'"),"\' are ignored when argument \"effects\" is not \"Wald\". \n")
        }
        out <- lapply(names(object$model), function(iBy){
            stats::coef(object$model[[iBy]], effects = effects, p = p[[iBy]], transform.sigma = transform.sigma, 
                        transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names, 
                        simplify = simplify, options = options)
        })
        if(ordering == "contrast"){
            browser()
        }
    }

    ## ** export
    return(out)
}

## * coef.rbindWald_lmm
##' @title Extract Coefficients From Combined Wald Tests Applied to Linear Mixed Models
##' @description Combine estimated values across linear contrasts applied on parameters from different linear mixed models. 
##'
##' @param object a \code{rbindWald_lmm} object.
##' @param effects [character] should the linear contrasts involved in the Wald test be output (\code{"Wald"}),
##' or the value of the linear mixed model parameters (\code{"all"})?
##' @param method [character vector] should the estimated value for the linear contrasts be output (one of \code{"none"}, \code{"bonferroni"}, ..., \code{"single-step2"})
##' and/or pooled linear contrast estimate(s) (\code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.gls1"}, \code{"pool.rubin"}, \code{"p.rejection"})?
##' Only relevant when \code{effects = "Wald"}.
##' @param ordering [character] should the output be ordered by name of the linear contrast (\code{"contrast"}) or by model (\code{"model"}).
##' @param backtransform [logical] should the estimates be back-transformed?
##' @param transform.names [logical] should the name of the coefficients be updated to reflect the transformation that has been used?
##' Only relevant when \code{effects="all"}.
##' @param simplify [logical] should the output be a vector or a list with one element specific to each possible ordering (i.e. contrast or model).
##' Only relevant when argument \code{method} refers to multiple comparisons and not to a pooling method.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details Argument \bold{effects}: when evaluating the proportion of rejected hypotheses (\code{effects="p.rejection"})
##' a \code{"single-step"} method will be used by default to evaluate the critical quantile.
##' This can be changed by adding adjustment method, e.g. \code{effects=c("bonferronin","p.rejection"}, in the argument.
##' 
##' @export
coef.rbindWald_lmm <- function(object, effects = "Wald", method = "none", ordering = NULL, transform.names = TRUE, backtransform = NULL, simplify = TRUE, ...){

    table.param <- stats::model.tables(object, effects = "param")
        
    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    pool.method <- options$pool.method
    adj.method <- options$adj.method

    ## *** object
    if(object$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling rbind.Wald_lmm. \n")
        return(invisible(NULL))
    }

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character. \n")
    }
    if(length(effects)!=1){
        stop("Argument \'effects\' must have length 1. \n")
    }
    valid.effects <- c("Wald","all")
    if(effects %in% valid.effects == FALSE){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }

    ## *** method
    if(!is.character(method) || !is.vector(method)){
        stop("Argument \'method\' must be a character. \n")
    }
    valid.method <- c("none",pool.method,adj.method)
    if(any(method %in% valid.method == FALSE)){
        stop("Incorrect value for argument \'method\': \"",paste(setdiff(method,valid.method), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.method, collapse ="\", \""),"\". \n")
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

    ## *** qt (hidden argument)
    if("p.rejection" %in% method){
        if(!is.null(attr(method,"qt"))){
            qt <- qt
            if(is.character(qt) && (qt %in% c("none","bonferroni","single-step","single-step2")==FALSE)){
                stop("Invalid \"qt\" attribe for argument \'method\'. \n",
                     "Should be one of \"none\", \"bonferroni\", \"single-step\", \"single-step2\". \n")
            }
        }else if(any(method %in% c("none",adj.method))){
            qt <- intersect(method, c("none",adj.method))
            if(is.character(qt) && (qt %in% c("none","bonferroni","single-step","single-step2")==FALSE)){
                stop("Incompatible values for argument \'method\': \"p.rejection\" and \"",qt,"\" \n",
                     "Consider using \"p.rejection\" and one of \"none\", \"bonferroni\", \"single-step\", \"single-step2\". \n")
            }
        }else{
            qt <- NULL
        }
    }else{
        qt <- NULL
    }
    
    ## ** extract from object
    if(effects == "all"){

        ## *** select value 
        if(backtransform){
            value.out <- table.param$value
        }else{
            value.out <- table.param$trans.value
        }

        ## *** rename
        if(backtransform || transform.names == FALSE){
            if(simplify==FALSE && !is.null(ordering)){
                newname <- table.param[[switch(ordering, "model" = "name", "contrast" = "model")]]                
            }else{
                newname <- table.param$Uname
            }
        }else{
            if(simplify==FALSE && !is.null(ordering)){
                newname <- table.param[[switch(ordering, "model" = "trans.name", "contrast" = "model")]]
            }else{
                newname <- table.param$trans.Uname
            }
        }
        names(value.out) <- newname

        ## *** ordering
        if(!is.null(ordering)){
            ordering.var <- switch(ordering, "model" = "model", "contrast" = "name")
            ordering.out <- factor(table.param[[ordering.var]], unique(table.param[[ordering.var]]))            
        }else if(simplify == FALSE){ ## no need for both as when both the attributes are not kept
            attr(value.out,"parameter") <- table.param$name
            attr(value.out,"model") <- table.param$model
        }

    }else if(effects == "Wald"){
        
        table.univariate <- object$univariate

        ## *** adjustment for multiple comparisons
        if(any(method %in% c("none",adj.method))){

            value.out <- coef.Wald_lmm(object, effects = "Wald", backtransform = backtransform, simplify = simplify, options = options)
            
            if(!is.null(ordering)){
                ordering.var <- switch(ordering, "model" = "model", "contrast" = "term")
                ordering.out <- factor(table.univariate[[ordering.var]], unique(table.univariate[[ordering.var]])) 
            }
            if(simplify == FALSE){
                attr(value.out,"parameter") <- table.univariate$term
                attr(value.out,"model") <- table.univariate$model
            }
        }else{

            value.out <- NULL

        }

        ## *** pooling
        if(any(method %in% pool.method)){
            table.pool <- pool.rbindWald_lmm(object, method = method, qt = qt, null = FALSE, level = NA, df = FALSE)
            value.out <- c(value.out,
                           stats::setNames(table.pool$estimate,rownames(table.pool))
                           )
        }
    }
    
    ## ** export
    if(simplify == FALSE && !is.null(ordering)){        
        out <- tapply(value.out, INDEX = ordering.out, FUN = base::identity, simplify = FALSE)
    }else if(!is.null(ordering)){
        out <- value.out[order(ordering.out)]
    }else{
        out <- value.out
    }
    return(out)
}

## * coef.Wald_lmm
##' @title Extract Coefficients From Wald Tests Applied to a Linear Mixed Model
##' @description Extract estimated value of linear contrasts applied on parameters from a linear mixed model. 
##'
##' @param object a \code{Wald_lmm} object.
##' @param effects [character] should the linear contrasts involved in the Wald test be output (\code{"Wald"}),
##' or the value of the linear mixed model parameters (\code{"all"})?
##' @param backtransform [logical] should the estimates be back-transformed?
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' Only relevant when \code{effects="all"}.
##' @param simplify [logical] omit from the output the attribute containing the type of each parameter or contrast (mu/sigma/k/rho).
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @export
coef.Wald_lmm <- function(object, effects = "Wald", backtransform = NULL, transform.names = TRUE, simplify = TRUE, ...){

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    adj.method <- options$adj.method

    ## *** object
    if(object$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling anova. \n")
        return(invisible(NULL))
    }

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character value. \n")
    }
    if(length(effects)!=1){
        stop("Argument \'effects\' must have length 1. \n")
    }
    valid.effects <- c("Wald","all")
    if(any(effects %in% valid.effects == FALSE)){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }        

    ## *** backtransform
    if(effects == "all"){
        if(is.null(backtransform)){
            backtransform <- FALSE
        }else if(!is.logical(backtransform)){
            stop("Argument \'backtransform\' must be TRUE or FALSE when argument \'effects\ is \"all\". \n")
        }
    }else if(effects == "contrast"){
        if(is.null(backtransform)){
            backtransform <- any(object$univariate$tobacktransform)
        }else if(is.character(backtransform)){
            backtransform <-  eval(parse(text=backtransform))
        }else if(is.numeric(backtransform)){
            backtransform <- as.logical(backtransform)
        }
    }

    ## *** simplify
    if(!is.numeric(simplify) && !is.logical(simplify)){
        stop("Argument \'simplify\' must be numeric or logical. \n")
    }
    if(length(simplify)!=1){
        stop("Argument \'simplify\' must have length 1. \n")
    }
    if(simplify %in% c(0,1) == FALSE){
        stop("Argument \'simplify\' must be TRUE/1 or FALSE/0. \n")
    }
    
    ## ** extract from object
    if(effects == "all"){
        table.param <- stats::model.tables(object, effects = "param")
        if(backtransform){
            out <- stats::setNames(table.param$value, table.param$name)
        }else{
            if(transform.names){
                out <- stats::setNames(table.param$trans.value, table.param$trans.name)
            }else{
                out <- stats::setNames(table.param$trans.value, table.param$name)
            }
        }
        if(!simplify){
            attr(out,"type") <- table.param$type
        }
    }else if(effects == "Wald"){
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
        if(!simplify){
            attr(out,"type") <- table.univariate$type
        }
    }

    ## ** export
    return(out)
}



##----------------------------------------------------------------------
### coef.R ends here
