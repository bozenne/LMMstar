### repetition.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  7 2024 (10:55) 
## Version: 
## Last-Updated: aug  7 2024 (13:18) 
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

## * repetition (documentation)
##' @title Number of Repetitions Within Cluster
##' @description Create a vector containing, for each line of the dataset, the number of occurences of the corresponding cluster up to the current line.
##' Can stratify the number of occurences on one or several variables.
##'
##' @param formula [formula] Specify the structure of the data: the time/repetition variable and the grouping variable,
##' e.g. ~1|id, ~ time|id, or ~time+region|id.
##' @param data [data.frame] dataset containing the observations.
##' @param format [character] the type of the output: a numeric vector (\code{"numeric"}), a character vector (\code{"character"}), or a factor vector (\code{"factor"}).
##' @param keep.time [logical] should the value of the time variable at the repetition be kept in the output (e.g. baseline.1, baseline.2, followUp.1 instead of 1,2,3).
##' Only relevant when argument \code{formula} contain a time/repetition variable and \code{format="character"} or \code{format="factor"}.
##' @param sep [character vector of length 2] character strings used to combine time variables (first element) and the name of the time variable with its values (second element).
##'
##' @return A numeric/character/factor vector, depending on argument \code{format}.
##'
##' @examples
##' data(sleepL, package = "LMMstar")
##' sleepL$dday <- paste0("d",sleepL$day)
##' sleepL$rep <- repetition(~1|id, data = sleepL)
##' sleepL$repDay <- repetition(~dday|id, data = sleepL)
##' sleepL$repDay.num <- repetition(~dday|id, data = sleepL, format = "numeric")
##' head(sleepL,15)


## * repetition (code)
##' @export
repetition <- function(formula, data,
                       format = "factor", keep.time = TRUE, sep = c(":",".")){

    ## ** check and normalize user input
    
    ## *** data
    if(!inherits(data,"data.frame")){
        stop("Incorrect type for argument \'data\': it must be or inherit from \"data.frame\". \n")
    }
    data <- as.data.frame(data)
    
    ## *** formula
    if(!inherits(formula,"formula")){
        stop("Incorrect type for argument \'formula\': it must be or inherit from \"formula\". \n")
    }

    terms.formula <- stats::terms(formula)
    var.formula <- all.vars(terms.formula)
    if(any(var.formula %in% names(data) == FALSE)){
        stop("Argument \'formula\' is inconsistent with argument \'data\'. \n",
             "Could not find column \"",paste(setdiff(var.formula,names(data)), collapse = "\" \""),"\" in argument \'data\'. \n", sep="")
    }
    if(length(setdiff(var.formula, all.vars(stats::delete.response(terms.formula))))>0){
        stop("There should not be any variable on the left hand side of argument \'formula\'. \n",
             "Should be a formula like ~1|cluster, ~time|cluster, or ~time+region|cluster. \n")
    }
    if(length(var.formula)<1){
        stop("Argument \'formula\' must contain at least a cluster variable. \n",
             "Should be a formula like ~1|cluster, ~time|cluster, or ~time+region|cluster. \n", sep="")
    }

    txt.formula <- deparse(formula)
    if(length(grepRaw(pattern = "|", txt.formula, all = TRUE, fixed = TRUE))!=1){
        stop("Argument \'formula\' must contain the symbold \"|\" exactly one. \n",
             "It separates the repetition and cluster variables, e.g.: ~time|cluster. \n")
    }

    ## *** format
    format <- match.arg(format, c("character","factor","numeric"))

    ## *** format
    if(length(sep)!=2){
        stop("Argument \'sep\' should have length 2. \n")
    }
    if(!is.character(sep) || !is.vector(sep)){
        stop("Argument \'sep\' should be a character vector. \n")
    }

    ## ** identify time and id
    var.time <- all.vars(formula[[2]][[2]])
    var.id <- all.vars(formula[[2]][[3]])
    if(length(var.id)==0){
        if(length(time)==0){
            stop("Argument \'formula\' must contain exactly one cluster variable. \n",
                 "Should be a formula like ~1|cluster. \n")
        }else{
            stop("Argument \'formula\' must contain exactly one cluster variable. \n",
                 "Should be a formula like ~",paste(var.time, collapse="+"),"|cluster. \n")
        }
    }else if(length(var.id)>1){
        if(length(time)==0){
            stop("Argument \'formula\' must contain exactly one cluster variable. \n",
                 "Should be a formula like ~1|",var.id[1],". \n")
        }else{
            stop("Argument \'formula\' must contain exactly one cluster variable. \n",
                 "Should be a formula like ~",paste(var.time, collapse="+"),"|",var.id[1],". \n")
        }
    }

    if(length(var.time)>0){
        if(all(sep %in% unlist(lapply(data[var.time],unique)))){
            warning("Values taken by the repetition variable(s) \"",paste(var.time, collapse ="\", \""),"\" include symbols \"",sep[1],"\" and \"",sep[2],"\" used in argument \'sep\'. \n",
                    "Can confuse the repetition function and lead to incorrect output. \n",
                    "Consider setting argument \'sep\' to another value, e.g. sep = c(\"",paste0("X",sep[1],"X"),"\",\"",paste0("X",sep[2],"X"),"\"). \n")
        }else if(sep[1] %in% unlist(lapply(data[var.time],unique))){
            warning("Values taken by the repetition variable(s) \"",paste(var.time, collapse ="\", \""),"\" include the symbol \"",sep[1],"\" used in argument \'sep\'. \n",
                    "Can confuse the repetition function and lead to incorrect output. \n",
                    "Consider setting argument \'sep\' to another value, e.g. sep = c(\"",paste0("X",sep[1],"X"),"\",\"",sep[2],"\"). \n")
        }else if(sep[2] %in% unlist(lapply(data[var.time],unique))){
            warning("Values taken by the repetition variable(s) \"",paste(var.time, collapse ="\", \""),"\" include symbol \"",sep[2],"\" used in argument \'sep\'. \n",
                    "Can confuse the repetition function and lead to incorrect output. \n",
                    "Consider setting argument \'sep\' to another value, e.g. sep = c(\"",sep[1],"\",\"",paste0("X",sep[2],"X"),"\"). \n")
        }
    }

    ## ** create repetition variable
    n.obs <- NROW(data)
    out <- numeric(length = n.obs)

    ls.index.id <- tapply(1:n.obs,data[[var.id]],identity, simplify = FALSE)
    if(length(var.time)==0){
        M.rep <- do.call(rbind,lapply(ls.index.id, function(iVec){cbind(iVec, 1:length(iVec))}))
        out[M.rep[,1]] <- M.rep[,2]
        if(format == "character"){
            out <- as.character(out)
        }else if(format == "factor"){
            out <- as.factor(out)
        }
        
    }else{
        vec.time <- nlme::collapse(data[var.time], sep = sep[1])
        table.idXtime <- table(data[[var.id]],vec.time)
        max.time <- apply(table.idXtime,2,max)
        
        out[unlist(ls.index.id)] <- do.call(base::c,lapply(ls.index.id, function(iVec){
            iRep <- do.call("+",lapply(unique(vec.time[iVec]), function(iTime){
                iCum <- cumsum(iTime == vec.time[iVec])
                iCum[diff(c(0,iCum))==0] <- 0
                iCum
            }))
            return(paste(vec.time[iVec],iRep,sep=sep[2]))
        }))

        if((format == "numeric" || format == "factor") || (keep.time == FALSE)){
            order <- unlist(mapply(iRep = max.time, iName = names(max.time), FUN = function(iRep,iName){paste(iName, 1:iRep, sep = sep[2])}))
            out <- factor(out, levels = order)
            if((format == "numeric") || (keep.time == FALSE)){
                out <- as.numeric(out)
            }
        }
        if(keep.time == FALSE){
            if(format == "character"){
                out <- as.character(out)
            }else if(format == "factor"){
                out <- as.factor(out)
            }
        }
    }

    ## ** export
    return(out)
}

##----------------------------------------------------------------------
### repetition.R ends here
