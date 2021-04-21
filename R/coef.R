### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: Apr 21 2021 (16:10) 
##           By: Brice Ozenne
##     Update #: 37
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * coef.lmm (code)
##' @export
coef.lmm <- function(object, effects = "all", type.object = "lmm", strata = NULL, transform = NULL){
    options <- LMMstar.options()
    
    ## ** normalize user imput
    type.object <- match.arg(type.object, c("lmm","gls"))
    if(identical(effects,"all")){
        effects <- c("mean","variance")
    }
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }
    if(is.null(transform)){
        transform <- options$transform
    }else if(transform %in% c(0,1,2) == FALSE){
        stop("Argument \'transform\' must be 0 (standard error parameters, correlation parameters), \n",
             "                               1 (log transformation of the standard error parameters, atanh transformation of the correlation parameters), \n",
             "                               2 (variance parameters, correlation parameters). \n")
    }

    ## ** extract
    if(type.object=="lmm"){

        out <- NULL
        if("mean" %in% effects){
            out <- c(out, object$param$mu)
        }
        if("variance" %in% effects){
            if(!is.null(object$param$sigma)){
                if(transform == 0){
                    out <- c(out, object$param$sigma)
                }else  if(transform == 1){
                    out <- c(out, log(object$param$sigma))
                }else if(transform == 2 && is.null(object$param$k)){
                    out <- c(out, object$param$sigma^2)
                }
            }
            if(!is.null(object$param$k)){
                if(transform == 0){
                    out <- c(out, object$param$k)
                }else if(transform==1){
                    out <- c(out, log(object$param$k))
                }else if(transform == 2){
                    strata.type <- object$param$type[object$param$type %in% c("sigma","k")]
                    strata.index <- object$param$strata[object$param$type %in% c("sigma","k")]
                    n.strata <- unique(strata.index)
                    for(iStrata in 1:n.strata){ ## iStrata <- 1 
                        iType <- strata.type[strata.index==iStrata]
                        iName.sigma <- names(iType[iType=="sigma"])
                        iName.k <- names(iType[iType=="k"])
                        out <- c(out, setNames(object$param$sigma[iName.sigma]^2*c(1,object$param$k[iName.k]^2),paste0(iName.sigma,":",object$time$var,object$time$levels)))
                    }
                }
            }
            if(!is.null(object$param$cor)){
                if(transform==1){
                    out <- c(out, atanh(object$param$cor))
                }else{
                    out <- c(out, object$param$cor)
                }
            }
        }
        if(!is.null(strata)){
            out <- out[object$param$strata[names(out)] %in% strata]
        }
        return(out)

    }else if(type.object=="gls"){
        if(!is.null(transform)){
            stop("Cannot handle argument \'transform\' when argument \'type.object\' is \"gls\". \n")
        }
        if(length(effects)!=1 || "variance" %in% effects){
            stop("Cannot handle argument \'effects\' when argument \'type.object\' is \"gls\". \n")
        }

        if(is.null(strata) && is.null(object$variable$strata)){
            return(coef(object$gls[[1]]))
        }else{
            return(lapply(object$gls[which(object$strata$levels %in% strata)], coef))
        }

    }
}


##----------------------------------------------------------------------
### coef.R ends here
