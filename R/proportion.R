### proportion.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 15 2022 (14:09) 
## Version: 
## Last-Updated: nov  8 2023 (16:02) 
##           By: Brice Ozenne
##     Update #: 36
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * proportion
#' @title Proportion of Significant Findings
#' @description Evaluate the proportion of test above the statistical significance level
#' 
#' @param object \code{Wald_lmm} object
#' @param n.sample [numeric,>=0] number of bootstrap sample used to assess the uncertainty.
#' If 0, then only the point estimate is computed.
#' @param trace [logical] shoudl the execution of the boostrap be trace.
#' @param ... additional arguments passed to \code{confint.Wald_lmm}
#'
#' @return a data.frame with the estimated proportion (estimate column), standard error and confidence interval (when boostrap is used).
#' 
#' @keywords utilities
#' 
#' @export
`proportion` <-
  function(object, n.sample, trace, ...) UseMethod("proportion")

## * proportion.mlmm
#' @export
proportion.mlmm <- function(object, n.sample = 100, trace = TRUE, ...){

    ## ** prepare
    if(n.sample>0){
        call <- attr(object,"call")
        data <- try(as.data.frame(eval(call$data)), silent = TRUE)
        if(inherits(data,"try-error")){
            stop("Could not extract the data from the call. \n",
                 data)
        }
        cluster.var <- attr(object$object$cluster.var,"original")
        cluster <- object$object$cluster
        n.cluster <- length(cluster)
        index.cluster <- split(1:NROW(data),data[[cluster.var]])
    }else if(n.sample < 0){
        stop("Argument \'n.sample\' must be a non-negative integer. \n")
    }else{
        trace <- FALSE
    }
    alpha <- attr(object,"level")
    

    ## ** warper
    warper <- function(iSample){ ## iSample <- 1
        if(iSample==0){
            iO <- object
        }else{
            iCluster <- cluster[sample.int(n.cluster, replace = TRUE)]
            iLs.indexCluster <- index.cluster[iCluster]
            iData <- data[unlist(iLs.indexCluster),,drop=FALSE]
            iData[[cluster.var]] <- unlist(lapply(1:n.cluster, function(iC){rep(paste0("C",iC), times =  length(iLs.indexCluster[[iC]]))}))
            iCall <- call
            iCall$data <- iData
            iO <- eval(iCall)
        }
        iCI <- confint(iO, columns = c("estimate","se","df","lower","upper","statistic","null","p.value"), ...)
        iTc <- attr(iCI, "quantile")
        iWald <- iCI$statistic
        iDf <- iCI$df
        iIntegral <- sapply(1:NROW(iCI), function(iStat){stats::pt(iTc - iWald[iStat], df = iDf[iStat]) - stats::pt(-iTc - iWald[iStat], df = iDf[iStat])})
        iOut <- 1-mean(iIntegral)
        if(iSample==0){
            attr(iOut,"level") <- attr(iCI, "level")
        }
        return(iOut)
    }
    
    ## ** iterate
    out.estimate <- rep(NA, n.sample+1)
    out.boot <- rep(NA, n.sample)

    if(trace){
        pb <- utils::txtProgressBar(max = n.sample+1, style =  3)
    }

    for(iSample in 0:n.sample){
        if(iSample==0){
            out.estimate <- warper(iSample)
        }else{
            out.boot[iSample] <- warper(iSample)
        }
        if(trace){
            utils::setTxtProgressBar(pb,iSample+1)
        }
    }

    if(trace){
        close(pb)
    }

    ## ** post process
    ## e.boot <- lapply(0:n.sample, warper)
    alpha <- 1-attr(out.estimate,"level")
    out <- data.frame(estimate = as.double(out.estimate),
                      se = stats::sd(out.boot, na.rm = TRUE),
                      df = NA,
                      lower = as.double(stats::quantile(out.boot,alpha/2, na.rm = TRUE)),
                      upper = as.double(stats::quantile(out.boot,1-alpha/2, na.rm = TRUE)),
                      p.value = NA)
    rownames(out) <- "proportion"
    return(out)
}

##----------------------------------------------------------------------
### proportion.R ends here
