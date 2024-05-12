### reformat.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 21 2023 (09:28) 
## Version: 
## Last-Updated: May 12 2024 (14:50) 
##           By: Brice Ozenne
##     Update #: 77
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .reformat (helper)
##' @description Convert results from \code{fitted.lmm}, \code{predict.lmm}, \code{residuals.lmm}
##' to the long or wide format.
##' @noRd
##' 
.reformat <- function(object, name, format, simplify,
                      keep.data, data, index.na,                      
                      object.cluster, index.cluster,
                      object.time, index.time,
                      call){

    sep <- LMMstar.options()$sep["reformat"]
    if(format == "long"){
        if(keep.data){
            out <- cbind(data, object)
        }else if(simplify>0 && NCOL(object)==1){
            out <- as.vector(object)
        }else{
            out <- as.data.frame(object)
        }

    }else if(format=="wide"){
        U.time <- object.time$levels
        if(is.matrix(object)){
            if(NCOL(object)>1){
                stop("Incorrect argument \'format\': cannot convert the result in the wide format with multiple types of ",name," \n")
            }
        }

        ## normalize cluster and time variables (in case of NAs)
        if(is.null(data) || (is.null(call$data) && is.null(call$newdata))){
            U.cluster <- object.cluster$levels
            
            if(length(index.na)==0){
                indexAll.cluster <- U.cluster[index.cluster]
                indexAll.time <- factor(U.time[index.time], U.time)
            }else{                
                indexAll.cluster <- rep(NA, NROW(object))
                indexAll.cluster[-unname(index.na)] <- U.cluster[match(index.cluster, object.cluster$index)]
                indexAll.cluster[unname(index.na)] <- as.character(attr(index.na,"cluster"))

                indexAll.time <- factor(rep(NA, NROW(object)), U.time)
                indexAll.time[-unname(index.na)] <- U.time[match(index.time, object.time$index)]
                indexAll.time[unname(index.na)] <- as.character(attr(index.na,"time"))
            }           
            
        }else{
            if(!is.null(data[[object.cluster$var]])){
                indexAll.cluster <- data[[object.cluster$var]]
            }else{
                indexAll.cluster <- index.cluster
            }
            indexAll.time <- U.time[index.time]
        }

        ## move to the wide format
        object2list <- list(as.vector(object))
        if(!is.null(name)){
            names(object2list) <- name
            sep <- ""
        }else{
            names(object2list) <- colnames(object)
        }

        df.object <- cbind(data.frame(cluster = indexAll.cluster, XXtimeXX = indexAll.time), object2list, stringsAsFactors = FALSE)
        names(df.object)[1] <- object.cluster$var
        if((simplify <= 0) && any(U.time %in% unique(df.object$XXtimeXX) == FALSE)){
            ## add missing times
            df.object <- rbind(df.object,
                               cbind(data.frame(cluster = indexAll.cluster[1], XXtimeXX = setdiff(U.time, unique(df.object$XXtimeXX))),
                                     stats::setNames(list(NA), names(object2list)))
                               )
        }
        out <- stats::reshape(data = df.object[order(factor(df.object$XXtimeXX, levels = object.time$levels)),], direction = "wide",
                              timevar = "XXtimeXX", idvar = object.cluster$var, v.names = names(object2list), times = U.time, sep = sep)
        if(!is.null(name)){ ## in case the user specify name <- " " to only keep the time levels (otherwise leads to " time1" as column names instead of "time1")
            names(out)[-1] <- trimws(names(out)[-1], which = "left")
        }

        ## use nicer column names
        if(simplify>0 && length(object.time$levels)==1){
            names(out)[-1] <- names(object2list)
        }

        ## restaure covariates
        keep.var <- attr(keep.data,"var")
        if(keep.data && length(keep.var)>0){
            out <- merge(x = unique(data[,c(object.cluster$var, keep.var),drop=FALSE]),
                         y = out,
                         by = object.cluster$var,
                         sort = FALSE)
        }
    }

    rownames(out) <- NULL
    return(out)
}


##----------------------------------------------------------------------
### reformat.R ends here
