### vcov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:28) 
## Version: 
## Last-Updated: Apr 16 2021 (12:36) 
##           By: Brice Ozenne
##     Update #: 34
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * vcov.lmm
##' @export
vcov.lmm <- function(object, effects = "all", type = "lmm", strata = NULL, data = NULL, p = NULL, transform = NULL){
    options <- LMMstar.options()
    x.transform <- attr(e.lmm$design$X.var, "transform")

    ## ** normalize user imput
    type <- match.arg(type, c("lmm","gls"))
    if(identical(effects,"all")){
        effects <- c("mean","variance")
    }
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }
    if(is.null(transform)){transform <- options$transform}

    ## ** extract or recompute variance covariance matrix
    if(type=="lmm"){

        if(!is.null(data) || !is.null(p) || (transform == x.transform)){
            vcov <- solve(information(object, data = data, p = p, transform = transform))
        }else{
            vcov <- object$vcov
        }

        keep.name <- names(coef(object, effects = effects, type = type, strata = strata))
        return(vcov[keep.name,keep.name,drop=FALSE])

    }else if(type=="gls"){
        if(!is.null(data)){
            stop("Cannot handle argument \'data\' when argument \'type\' is \"gls\". \n")
        }
        if(!is.null(p)){
            stop("Cannot handle argument \'p\' when argument \'type\' is \"gls\". \n")
        }
        if(!is.null(transform)){
            stop("Cannot handle argument \'transform\' when argument \'type\' is \"gls\". \n")
        }
        
        .getVcov <- function(oo, effects){
            if(identical(effects,"mean")){
                return(vcov(oo))
            }else if(identical(effects,"variance")){
                return(oo$apVar)
            }else{
                out.names <- c(colnames(vcov(oo)),colnames(oo$apVar))
                out <- as.matrix(Matrix::bdiag(vcov(oo), oo$apVar))
                dimnames(out) <- list(out.names,out.names)
                return(out)
            }
        }

        if(is.null(strata) && (object$strata$n == 1)){
            .getVcov(object$gls[[1]], effects = effects)
        }else{
            return(lapply(object$gls[which(object$strata$levels %in% strata)], .getVcov, effects = effects))
        }

    }
}

##----------------------------------------------------------------------
### vcov.R ends here
