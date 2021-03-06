### LMMstar.options.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 16 2021 (12:01) 
## Version: 
## Last-Updated: jul  7 2021 (18:02) 
##           By: Brice Ozenne
##     Update #: 53
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
#' @include 0-onload.R
#'
#' @description Update or select global options for the LMMstar package.
#'
#' @param ... options to be selected or updated
#' @param reinitialise should all the global parameters be set to their default value
#'
#' @details The options are: \itemize{
#' \item backtransform.confint [logical]: should variance/covariance/correlation estimates be back-transformed when they are transformed on the log or atanh scale. Used by \code{confint}.
#' \item columns.confint [character vector]: columns to ouput when displaying the model coefficients using \code{confint}.
#' \item columns.summary [character vector]: columns to ouput when displaying the model coefficients using \code{summary}.
#' \item df [logical]: should approximate degrees of freedom be computed for Wald and F-tests. Used by \code{lmm}, \code{anova}, \code{predict}, and \code{confint}.
#' \item drop.X [logical]: should columns causing non-identifiability of the model coefficients be dropped from the design matrix. Used by \code{lmm}.
#' \item method.fit [character]: objective function when fitting the Linear Mixed Model (REML or ML). Used by \code{lmm}.
#' \item method.numDeriv [character]: type used to approximate the third derivative of the log-likelihood (when computing the degrees of freedom). Can be \code{"simple"} or \code{"Richardson"}. See \code{numDeriv::jacobian} for more details. Used by \code{lmm}.
#' \item precompute.moments [logical]: Should the cross terms between the residuals and design matrix be pre-computed. Useful when the number of subject is substantially larger than the number of mean paramters.
#' \item trace [logical]: Should the progress of the execution of the \code{lmm} function be displayed?
#' \item tranform.sigma, tranform.k, tranform.rho: transformation used to compute the confidence intervals/p-values for the variance and correlation parameters. See the detail section of the coef function for more information.
#' Used by \code{lmm}, \code{anova} and \code{confint}.
#' \item type.information [character]: Should the expected or observed information (\code{"expected"} or \code{"observed"}) be used to perform statistical inference? Used by \code{lmm}, \code{anova} and \code{confint}.
#' }
#' @export
LMMstar.options <- function(..., reinitialise = FALSE){
  
    if (reinitialise == TRUE) {
        assign(".LMMstar-options", 
               list(backtransform.confint = TRUE,
                    columns.confint = c("estimate","se","df","lower","upper","p.value"),
                    columns.summary = c("estimate","se","df","lower","upper","p.value",""),
                    df = TRUE,
                    drop.X = TRUE,
                    precompute.moments = TRUE,
                    method.fit = "REML",
                    method.numDeriv = "simple", 
                    trace = FALSE,
                    transform.sigma = "log",
                    transform.k = "log",
                    transform.rho = "atanh",
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

          ok.column <- c("estimate","se","statistic","df","lower","upper","null","p.value")
          if("columns.confint" %in% names(args) && any(args$columns.confint %in% ok.column == FALSE)){
              stop("Argument \'columns.confint\' must be a character vector with values among \"",paste(ok.column, collapse = "\" \""),"\". \n")
          }
          if("columns.summary" %in% names(args) && any(args$columns.summary %in% c(ok.column,"") == FALSE)){
              stop("Argument \'columns.summary\' must be a character vector with values among \"",paste(c(ok.column,""), collapse = "\" \""),"\". \n")
          }
          if("df" %in% names(args) && !is.logical(args$df)){
              stop("Argument \'df\' must be of type logical. \n")
          }
          if("drop.X" %in% names(args) && !is.logical(args$drop.X)){
              stop("Argument \'drop.X\' must be of type logical. \n")
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
