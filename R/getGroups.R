### getGroups.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 14 2026 (15:55) 
## Version: 
## Last-Updated: apr 15 2026 (15:24) 
##           By: Brice Ozenne
##     Update #: 49
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * getGroups.structure (documentation)
##' @title Extract Information About a Pattern
##' @description Extract information (cluster, strata, data, design matrix, linear predictor) about covariance, variance, or correlation pattern.
##' 
##' @param object [structure] ID, IND, CS, UN, DUN, ...
##' @param form [character] should information be extracted w.r.t a variance-covariance pattern (\code{"covariance"}),
##' a variance pattern (\code{"variance"}), or a correlation pattern (\code{"correlation"}).
##' @param level [character] name of the pattern.
##' @param data [character] information to extract: \itemize{
##' \item \code{"cluster"}: index of the clusters with the pattern
##' \item \code{"data"}: dataset relative to the pattern (number of rows: number of repetitions for a specific cluster with the pattern)
##' \item \code{"lp"}: linear predictor relative to the pattern (length: number of repetitions for a specific cluster with the pattern)
##' \item \code{"strata"}: index of the strata for the pattern
##' \item \code{"X"}: design matrix relative to the pattern (number of rows: number of repetitions for a specific cluster with the pattern)
##' \item \code{"Xpattern"}: array with the name of the variance-covariance parameters to be multiplied to obtain the residual variance-covariance matrix.
##' @param sep Not used, only for compatibilty with the generic method.
##' 

## * getGroups.structure (code)
##' @export
getGroups.structure <- function(object, form, level, data, sep){

    ## ** check input

    ## *** form
    form <- match.arg(form, c("covariance","variance","correlation"))
    form.norm <- switch(form,
                        "covariance" = "name",
                        "variance" = "var",
                        "correlation" = "cor")

    ## *** data
    data <- match.arg(data, c("cluster","data","lp","pattern","strata","X","Xpattern"))
    if(form=="covariance"){
        if(data == "data"){
            stop("Argument \'form\' can only be \"variance\" or \"correlation\" when argument \'data\' equals \"data\". \n")
        }
        if(data != "pattern"){
            level.var <- object$Upattern[object$Upattern$name == level,"var"]
            level.cor <- object$Upattern[object$Upattern$name == level,"cor"]
        }
    }
    
    ## *** level
    if(missing(level)){
        if(data != "pattern"){
            stop("Missing value for argument \'level\'. \n")
        }
    }else{

        if(data == "pattern"){
            warning("Value for argument \'level\' ignored when argument \'data\' equals \"pattern\". \n")
        }

        if(any(level %in% object$Upattern[[form.norm]] == FALSE)){
            Ulevel <- unique(object$Upattern[[form.norm]])
            if(all(level %in% 1:length(Ulevel))){
                level <- Ulevel[match(level,1:length(Ulevel))]
            }else{
                stop("Incorrect value for argument \'level\'. \n",
                     "Valid values: \"",paste(Ulevel, collapse="\", \""),"\". \n")
            }
        }

        if(length(level)>1 && data %in% c("Xpattern","lp","X")){
            stop("Argument \'level\' should have length 1 when argument \'data\' equal \"Xpattern\", \"lp\", or \"X\". \n")
        }
    }
    
    ## ** extract information
    if(data == "pattern"){
        out <- object$Upattern[[form.norm]]
    }else if(data == "strata"){
        out <- object$Upattern[object$Upattern[[form.norm]] %in% level,"index.strata"]
    }else if(data == "cluster"){
        out <- unlist(object$Upattern[object$Upattern[[form.norm]] %in% level,"index.cluster"])
    }else if(data=="Xpattern"){
        if(form=="covariance"){
            out <- array(c(object$var$Xpattern[[level.var]],
                           object$cor$Xpattern[[level.cor]]),
                         dim = c(dim(object$var$Xpattern[[level.var]])[1:2], dim(object$var$Xpattern[[level.var]])[3] + dim(object$cor$Xpattern[[level.var]])[3]),
                         dimnames = list(NULL,NULL,c(dimnames(object$var$Xpattern[[level.var]])[[3]], dimnames(object$cor$Xpattern[[level.var]])[[3]])))
        }else{ ## variance or correlation
            out <- object[[form.norm]]$Xpattern[[level]]
        }
    }else if(data=="lp"){
        if(form=="covariance"){
            out <- cbind(var = object$var$pattern2lp[[level.var]],
                         cor = object$cor$pattern2lp[[level.cor]])
        }else{
            out <- object[[form.norm]]$pattern2lp[[level]]
        }
    }else if(data=="X"){
        if(form=="covariance"){
            out <- cbind(object$var$lp2X[object$var$pattern2lp[[level.var]],,drop=FALSE],
                         object$cor$lp2X[object$cor$pattern2lp[[level.cor]],,drop=FALSE])
        }else{
            out <- object[[form.norm]]$lp2X[object[[form.norm]]$pattern2lp[[level]],,drop=FALSE]
        }
    }else if(data=="data"){
        out <- object[[form.norm]]$lp2data[object[[form.norm]]$pattern2lp[[level]],,drop=FALSE]
    }    

    ## ** export
    return(out)
}

##----------------------------------------------------------------------
### getGroups.R ends here
