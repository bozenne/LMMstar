### iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  4 2021 (10:04) 
## Version: 
## Last-Updated: Jun  7 2021 (14:34) 
##           By: Brice Ozenne
##     Update #: 14
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

iid.lmm <- function(object, effects = "all", type.information = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

    options <- LMMstar.options()
    x.transform.sigma <- object$reparametrize$transform.sigma
    x.transform.k <- object$reparametrize$transform.k
    x.transform.rho <- object$reparametrize$transform.rho

    ## ** normalize user imput
    dots <- list(...)
    dots$complete <- NULL ## for multcomp which passes an argument complete when calling vcov
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","variance","correlation"), several.ok = TRUE)

    if(is.null(type.information)){
        type.information <- options$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, options = options,
                            x.transform.sigma = x.transform.sigma, x.transform.k = x.transform.k, x.transform.rho = x.transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform

    ## ** get information and score
    object.vcov <- stats::vcov(object, effects = effects, type.information = type.information, df = FALSE,
                               transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
    object.score <- lava::score(object, effects = effects, indiv = TRUE, 
                                transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)

    ## ** compute iid
    out <- object.score %*% object.vcov
    return(out)
}
##----------------------------------------------------------------------
### iid.R ends here
