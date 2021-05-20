### LMMstar.options.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 16 2021 (12:01) 
## Version: 
## Last-Updated: May 19 2021 (15:00) 
##           By: Brice Ozenne
##     Update #: 22
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Function LMMstar.options
## * Documentation - BuyseTest.options
#' @title Global options for LMMstar package
#' @name LMMstar.options
#' @include 0-onLoad.R
#'
#' @description Update or select global options for the LMMstar package.
#'
#' @param ... options to be selected or updated
#' @param reinitialise should all the global parameters be set to their default value
#'         
#' @export
LMMstar.options <- function(..., reinitialise = FALSE){
  
    if (reinitialise == TRUE) {
        assign(".LMMstar-options", 
               list(df = TRUE,
                    method.fit = "REML",
                    method.numDeriv = "simple", 
                    transform.sigma = "none",
                    transform.k = "none",
                    transform.rho = "none",
                    type.information = "observed"),
               envir = LMMstar.env)
    
    return(invisible(get(".LMMstar-options", envir = LMMstar.env)))
    
  }else{
    
    args <- list(...)
    object <- get(".LMMstar-options", envir = LMMstar.env)

      if(length(args)==0){ ## read all
          return(object)
      }else if (!is.null(names(args))) { ## write

          if(any(names(args) %in% names(object) == FALSE)){
              stop("Incorrect element selected: \"",paste0(names(args)[names(args) %in% names(object) == FALSE], collapse = "\" \""),"\"\n",
                   "Available elements: \"",paste0(setdiff(names(object),names(args)), collapse = "\" \""),"\"\n")
          }

          if("df" %in% names(args) && !is.logical(args$df)){
              stop("Argument \'df\' must be of type logical. \n")
          }
          if("method.fit" %in% names(args)){
              args$method.fit <- match.arg(args$method.fit, c("ML","REML"))
          }
          if("method.numDeriv" %in% names(args)){
              args$method.numDeriv <- match.arg(args$method.numDeriv, c("simple","Richardson","complex"))
          }
          if("transform.sigma" %in% names(args)){
              args$transform.sigma <- match.arg(args$transform.sigma, c("none","log","square","logsquare"))
          }
          if("transform.k" %in% names(args)){
              args$transform.k <- match.arg(args$transform.k, c("none","log","square","logsquare","sd","logsd","var","logvar"))
          }
          if("transform.rho" %in% names(args)){
              args$transform.rho <- match.arg(args$transform.rho, c("none","atanh","cov"))
          }
          if("type.information" %in% names(args)){
              args$type.information <- match.arg(args$type.information, c("expected","observed"))
          }
          object[names(args)] <- args
      
          assign(".LMMstar-options", 
                 object, 
                 envir = LMMstar.env)
      
          return(invisible(object))
      
      } else {# read
          args <- unlist(args)
          if(any(args %in% names(object) == FALSE)){
              stop("Incorrect element selected: \"",paste0(args[args %in% names(object) == FALSE], collapse = "\" \""),"\"\n",
                   "Available elements: \"",paste0(setdiff(names(object),args), collapse = "\" \""),"\"\n")
          }
          return(object[args])
      }
    
  }
  
  
}


##----------------------------------------------------------------------
### LMMstar.options.R ends here
