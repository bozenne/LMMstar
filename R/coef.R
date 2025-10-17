### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: okt 17 2025 (16:01) 
##           By: Brice Ozenne
##     Update #: 1403
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
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see detail.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
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

## * coef.Wald_lmm
##' @title Extract Coefficients From Wald Tests
##' @description Extract estimated value of linear contrasts involved in Wald tests. 
##'
##' @param object a \code{Wald_lmm} object.
##' @param effects [character] should the linear contrasts involved in the Wald test be output (\code{"Wald"}),
##' or the value of the linear mixed model parameters (\code{"all"})?
##' @param method [character] how linear contrast estimates should be pooled (\code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.gls1"}, \code{"pool.rubin"}, \code{"p.rejection"})?
##' Only relevant when \code{effects = "Wald"}. See \code{\link{confint.Wald_lmm}} for details.
##' @param backtransform [logical] should the estimates be back-transformed?
##' @param transform.sigma,transform.k,transform.rho [character] for internal use (delta-method via estimate).
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' Only relevant when \code{effects="all"}.
##' @param simplify [logical] omit from the output an attribute containing the parameter type or contrast term.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @export
coef.Wald_lmm <- function(object, effects = "Wald", method = "none", backtransform = NULL,
                          transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE,
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

    ## *** method
    if(effects == "all"){
        if(!identical(method,"none")){
            message("Argument \'method\' is ignored when argument \'effects\ is \"all\". \n")
        }
        table.param <- stats::model.tables(object, effects = "param")
    }

    ## *** backtransform
    if(effects == "all"){
        if(!is.null(backtransform) && (length(backtransform)!=1 || any(!is.logical(backtransform)) && any(backtransform %in% 0:1 == FALSE))){
            stop("Argument \'backtransform\' must be TRUE or FALSE when argument \'effects\ is \"all\". \n")
        }
        if(identical(transform.sigma, object$args$transform.sigma) && identical(transform.k, object$args$transform.k) && identical(transform.rho, object$args$transform.rho)){
            backtransform <- FALSE
            transform.sigma <- NULL
            transform.k <- NULL
            transform.rho <- NULL
        }else if(identical(transform.sigma, "none") && identical(transform.k, "none") && identical(transform.rho, "none")){
            backtransform <- TRUE
            transform.sigma <- NULL
            transform.k <- NULL
            transform.rho <- NULL
        }else if(is.null(backtransform) && is.null(transform.sigma) && is.null(transform.k) && is.null(transform.rho)){
            backtransform <- TRUE
        }
    }else if(effects == "Wald"){
        if(is.null(backtransform)){
            backtransform <- any(object$univariate$tobacktransform)
        }else if(is.character(backtransform)){
            backtransform <-  eval(parse(text=backtransform))
        }else if(is.numeric(backtransform)){
            backtransform <- as.logical(backtransform)
        }
    }

    ## *** transform
    if("Wald" %in% effects && !is.null(transform.names) && transform.names!=TRUE){
        message("Argument \'transform.names\' is ignored when argument \"effects\" is \"Wald\". \n")
    }
    if(!is.null(transform.sigma) || !is.null(transform.k) || !is.null(transform.rho)){
        if(effects == "Wald"){
            txt.print <- paste0("transform.",c("sigma","k","rho"))[c(!is.null(transform.sigma),!is.null(transform.k),!is.null(transform.rho))]
            message("Argument(s) \'",paste(txt.print, collapse = "\', \'"),"\' are ignored when argument \'effects\' equals to \"contrast\". \n")
        }
        if(object$args$p.null){
            p <- stats::setNames(table.param$value, table.param$name)
        }else{
            p <- NULL
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
        if(!is.null(transform.sigma) || !is.null(transform.k) || !is.null(transform.rho)){
            out <- coef(lmm(object), effects = "all", p = p,
                        transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names, simplify = simplify)
        }else if(backtransform){
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

        Mout <- confint(object, method = method, backtransform = backtransform, column = c("term","estimate"))
        out <- stats::setNames(Mout[,"estimate"], rownames(Mout))
        if(!simplify){
            attr(out,"term") <- Mout[,"term"]
        }
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

## * coef.rbindWald_lmm
##' @title Extract Coefficients From Combined Wald Tests
##' @description Combine estimated values across linear contrasts applied on parameters from different linear mixed models. 
##'
##' @param object a \code{rbindWald_lmm} object.
##' @param effects [character] By default will output the estimates relative to the hypotheses being tested (\code{"Wald"}).
##' But can also output all model coefficients (\code{"all"}),
##' @param method [character vector] should the estimated value for the linear contrasts be output (one of \code{"none"}, \code{"bonferroni"}, ..., \code{"single-step2"})
##' and/or pooled linear contrast estimate(s) (\code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.gls1"}, \code{"pool.rubin"}, \code{"p.rejection"})?
##' Only relevant when \code{effects = "Wald"}.
##' @param ordering [character] should the output be ordered by name of the linear contrast (\code{"contrast"}) or by model (\code{"model"}).
##' @param backtransform [logical] should the estimate be back-transformed?
##' @param transform.sigma,transform.k,transform.rho [character] for internal use (delta-method via estimate).
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' Ignored when \code{effects="Wald"}.
##' @param simplify [logical] omit from the output an attribute containing the parameter type/model or contrast term.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @keywords methods
##' @return A numeric vector
##' 
##' @export
coef.rbindWald_lmm <- function(object, effects = "Wald", method = "none", ordering = NULL, backtransform = NULL,
                               transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE,
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
    if(effects == "all" && !identical(method,"none")){
        message("Argument \'method\' is ignored when argument \'effects\ is \"all\". \n")
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
    if(effects == "all"){
        if(!is.null(backtransform) && (length(backtransform)!=1 || any(!is.logical(backtransform)) && any(backtransform %in% 0:1 == FALSE))){
            stop("Argument \'backtransform\' must be TRUE or FALSE when argument \'effects\ is \"all\". \n")
        }
        if(identical(transform.sigma, object$args$transform.sigma) && identical(transform.k, object$args$transform.k) && identical(transform.rho, object$args$transform.rho)){
            backtransform <- FALSE
            transform.sigma <- NULL
            transform.k <- NULL
            transform.rho <- NULL
        }else if(identical(transform.sigma, "none") && identical(transform.k, "none") && identical(transform.rho, "none")){
            backtransform <- TRUE
            transform.sigma <- NULL
            transform.k <- NULL
            transform.rho <- NULL
        }else if(is.null(backtransform)){
            if(is.null(transform.sigma) && is.null(transform.k) && is.null(transform.rho)){
                backtransform <- TRUE
            }else{
                backtransform <- all(c(transform.sigma,transform.k,transform.rho) == "none")
            }
        }
    }else if(effects == "Wald"){
        if(is.null(backtransform)){
            backtransform <- any(object$univariate$tobacktransform)
        }else if(is.character(backtransform)){
            backtransform <-  eval(parse(text=backtransform))
        }else if(is.numeric(backtransform)){
            backtransform <- as.logical(backtransform)
        }
    }

    ## *** transform.names
    if("Wald" %in% effects && !is.null(transform.names) && transform.names!=TRUE){
        message("Argument \'transform.names\' is ignored when argument \"effects\" is \"Wald\". \n")
    }
    if(!is.null(transform.sigma) || !is.null(transform.k) || !is.null(transform.rho)){
        if(effects == "Wald"){
            txt.print <- paste0("transform.",c("sigma","k","rho"))[c(!is.null(transform.sigma),!is.null(transform.k),!is.null(transform.rho))]
            message("Argument(s) \'",paste(txt.print, collapse = "\', \'"),"\' are ignored when argument \'effects\' equals to \"contrast\". \n")
        }
    }
    
    ## *** simplify
    if(!is.numeric(simplify) && !is.logical(simplify)){
        stop("Argument \'simplify\' must be numeric or logical. \n")
    }
    if(length(simplify)!=1){
        stop("Argument \'simplify\' must have length 1. \n")
    }
    if(!is.logical(simplify) && simplify %in% 0:1 == FALSE){
        stop("Argument \'simplify\' must be TRUE or FALSE. \n")
    }
    
    ## ** extract from object
    if(effects == "all"){
        table.param <- stats::model.tables(object, effects = "param", ordering  = ordering)

        if(backtransform){
            out <- stats::setNames(table.param$value, table.param$Uname)
        }else if(transform.names==FALSE){
            out <- stats::setNames(table.param$trans.value, table.param$Uname)
        }else{
            out <- stats::setNames(table.param$trans.value, table.param$trans.Uname)
        }

        if(!simplify){
            attr(out,"type") <- table.param$type
            attr(out,"model") <- table.param$model
        }

    }else if(effects == "Wald"){

        Mout <- confint(object, method = method, ordering = ordering,
                        backtransform = backtransform, column = c("term","estimate"))
        out <- stats::setNames(Mout[,"estimate"], rownames(Mout))
        if(!simplify){
            attr(out,"term") <- Mout[,"term"]
        }
        
    }
    
    ## ** export
    return(out)
}


## * coef.mlmm (documentation)
##' @title Extract Coefficients From Multiple Linear Mixed Models
##' @description Combine estimated parameters or linear contrasts applied on parameters from group-specific linear mixed models.
##'
##' @param object a \code{mlmm} object.
##' @inheritParams coef.rbindWald_lmm
##' 
##' @keywords methods
##' @return A numeric vector

## * coef.mlmm (code)
##' @export
coef.mlmm <- coef.rbindWald_lmm



##----------------------------------------------------------------------
### coef.R ends here
