### LMMstar.options.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 16 2021 (12:01) 
## Version: 
## Last-Updated: May 14 2021 (14:41) 
##           By: Brice Ozenne
##     Update #: 14
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
               list(transform.sigma = "none",
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

          if("type" %in% names(args)){
              args$type <- match.arg(args$type, c("expected","observed"))
          }
          if("transform.sigma" %in% names(args)){
              args$transform.sigma <- match.args(args$transform.sigma, c("none","log","square"))
          }
          if("transform.k" %in% names(args)){
              args$transform.k <- match.args(args$transform.k, c("none","sd"))
          }
          if("transform.rho" %in% names(args)){
              args$transform.rho <- match.args(args$transform.rho, c("none","atanh"))
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
