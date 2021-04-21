### getVarCov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:57) 
## Version: 
## Last-Updated: Apr 20 2021 (15:55) 
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

## * getVarCov.lmm
getVarCov.lmm <- function(object, individual = NULL, type.object = c("lmm","gls"), simplifies = TRUE, strata = NULL){

    ## ** normalize user imput
    type.object <- match.arg(type.object, c("lmm","gls"))
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }

    if(!is.null(individual)){
        individual <- match.arg(individual, object$cluster$levels, several.ok = TRUE)
    }

    if(type.object == "lmm"){
    
        if(is.null(individual)){
            index.fulltime <- which(sapply(attr(object$design$X.var,"Upattern.time"),identical,as.character(object$time$levels)))
            if(!is.null(strata)){
                index.strata <- which(attr(object$design$X.var,"UX.strata") %in% strata)
            }else{
                index.strata <- 1:attr(object$design$X.var,"nUpattern")
                strata <- attr(object$design$X.var,"UX.strata")
            }
            out <- setNames(object$Omega[intersect(index.fulltime,index.strata)],strata)
        }else{
            out <- object$Omega[object$design$index.vargroup[individual]]
        }
        if(is.list(out) && length(out)==1 && simplifies){
            return(out[[1]])
        }else{
            return(out)
        }
        
    }else if(type.object == "gls"){
        if(object$strata$n==1){
            if(is.null(individual)){
                return(getVarCov(object$gls[[1]]))
            }else{
                return(getVarCov(object$gls[[1]], individual = individual))
            }
        }else{
            if(is.null(individual)){
                return(lapply(object$gls, getVarCov))
            }else{
                out <- setNames(vector(mode = "list", length = length(individual)),individual)
                for(iStrata in 1:object$strata$n){ ## iStrata <- 1
                    iIndiv <- intersect(individual,names(object$design$index.strata[object$design$index.strata==iStrata]))
                    if(length(iIndiv)>0){
                        out[match(iIndiv,names(out))] <- getVarCov(object$gls[[iStrata]], individual = iIndiv)
                    }
                }
                return(out)
            }
        }

    }
        
}

##----------------------------------------------------------------------
### getVarCov.R ends here
