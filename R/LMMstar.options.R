### LMMstar.options.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 16 2021 (12:01) 
## Version: 
## Last-Updated: Apr 20 2021 (22:46) 
##           By: Brice Ozenne
##     Update #: 10
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
               list(transform.sigmaDeriv = FALSE,
                    type.information = "expected"), 
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
          if(("transform" %in% names(args)) && (args$transform %in% 0:2 == FALSE)){
              stop("Argument \'transform\' must be 0 (standard error parameters, correlation parameters), \n",
                   "                               1 (log transformation of the standard error parameters, atanh transformation of the correlation parameters), \n",
                   "                               2 (variance parameters, correlation parameters). \n")
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
