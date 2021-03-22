### vcov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:28) 
## Version: 
## Last-Updated: mar 22 2021 (22:21) 
##           By: Brice Ozenne
##     Update #: 23
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
vcov.lmm <- function(object, effects = "mean", type = "lmm", strata = NULL, data = NULL, p = NULL){

    ## ** normalize user imput
    type <- match.arg(type, c("lmm","gls"))
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }

    ## ** extract
    if(type=="lmm"){

        if(!is.null(data) || !is.null(p)){
            vcov <- solve(information(object, data = data, p = p))
        }else{
            vcov <- object$vcov
        }

        keep.name <- names(coef(object, effects = effects, type = type, strata = strata))
        return(vcov[keep.name,keep.name,drop=FALSE])

    }else if(type=="gls"){

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
