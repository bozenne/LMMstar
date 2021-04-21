### vcov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:28) 
## Version: 
## Last-Updated: Apr 20 2021 (15:53) 
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

## * vcov.lmm
##' @export
vcov.lmm <- function(object, effects = "all", type.object = "lmm", strata = NULL, data = NULL, p = NULL, transform = NULL, type.information = NULL){
    options <- LMMstar.options()
    x.transform <- attr(e.lmm$design$X.var, "transform")

    ## ** normalize user imput
    type.object <- match.arg(type.object, c("lmm","gls"))
    if(identical(effects,"all")){
        effects <- c("mean","variance")
    }
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }
    if(is.null(transform)){transform <- options$transform}
    if(is.null(type.information)){type.information <- options$type.information}

    ## ** extract or recompute variance covariance matrix
    if(type.object=="lmm"){

        if(!is.null(data) || !is.null(p) || (transform == x.transform)){
            vcov <- solve(information(object, data = data, p = p, transform = transform, type.information = type.information))
        }else{
            vcov <- object$vcov
        }

        keep.name <- names(coef(object, effects = effects, type.object = type.object, strata = strata))
        return(vcov[keep.name,keep.name,drop=FALSE])

    }else if(type.object=="gls"){
        if(!is.null(data)){
            stop("Cannot handle argument \'data\' when argument \'type.object\' is \"gls\". \n")
        }
        if(!is.null(p)){
            stop("Cannot handle argument \'p\' when argument \'type.object\' is \"gls\". \n")
        }
        if(!is.null(transform)){
            stop("Cannot handle argument \'transform\' when argument \'type.object\' is \"gls\". \n")
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
