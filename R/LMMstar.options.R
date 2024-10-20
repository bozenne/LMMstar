### LMMstar.options.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 16 2021 (12:01) 
## Version: 
## Last-Updated: okt 20 2024 (16:42) 
##           By: Brice Ozenne
##     Update #: 186
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * LMMstar.options (documentation) 
#' @title Global options for LMMstar package
#' @include 0-onload.R
#'
#' @description Update or select global options for the LMMstar package.
#'
#' @param ... options to be selected or updated
#' @param reinitialise should all the global parameters be set to their default value
#'
#' @details The options are: \itemize{
#' \item adj.method [character vector]: possible methods to adjust for multiple comparisons. NOT MEANT TO BE CHANGED BY THE USER.
#' \item backtransform.confint [logical]: should variance/covariance/correlation estimates be back-transformed when they are transformed on the log or atanh scale. Used by \code{confint}.
#' \item columns.anova [character vector]: columns to ouput when using \code{anova} with argument \code{ci=TRUE}.
#' \item columns.confint [character vector]: columns to ouput when using \code{confint}.
#' \item columns.summarize [character vector]: columns to ouput when displaying descriptive statistics using \code{summarize}.
#' \item columns.summary [character vector]: columns to ouput when displaying the model coefficients using \code{summary}.
#' \item df [logical]: should approximate degrees-of-freedom be computed for Wald and F-tests. Used by \code{lmm}, \code{anova}, \code{predict}, and \code{confint}.
#' \item drop.X [logical]: should columns causing non-identifiability of the model coefficients be dropped from the design matrix. Used by \code{lmm}.
#' \item effects [character]: parameters relative to which estimates, score, information should be output.
#' \item min.df [integer]: minimum possible degree-of-freedom. Used by \code{confint}.
#' \item method.fit [character]: objective function when fitting the Linear Mixed Model (REML or ML). Used by \code{lmm}.
#' \item method.numDeriv [character]: type used to approximate the third derivative of the log-likelihood (when computing the degrees-of-freedom). Can be \code{"simple"} or \code{"Richardson"}. See \code{numDeriv::jacobian} for more details. Used by \code{lmm}.
#' \item n.sampleCopula [integer]: number of samples used to compute confidence intervals and p-values adjusted for multiple comparisons via \code{"single-step2"}. Used by \code{confint.Wald_lmm}.
#' \item optimizer [character]: method used to estimate the model parameters. Either \code{"FS"}, an home-made fisher scoring algorithm, or a method from \code{optimx:optimx} like \code{"BFGS"}or \code{Nelder-Mead}.
#' \item param.optimizer [numeric vector]: default option for the \code{FS} optimization routine: maximum number of gradient descent iterations (\code{n.iter}), maximum acceptable score value (\code{tol.score}), maximum acceptable change in parameter value (\code{tol.param}), method to initialize the correlation parameters (\code{init.cor}).
#' \item pool.method [character vector]: possible methods to pool estimates. NOT MEANT TO BE CHANGED BY THE USER.
#' \item precompute.moments [logical]: Should the cross terms between the residuals and design matrix be pre-computed. Useful when the number of subject is substantially larger than the number of mean paramters.
#' \item sep [character vector]: character used to combined two strings of characters in various functions (lp: .vcov.model.matrix, k.cov/k.strata: .skeletonK, pattern: .findUpatterns, rho.name/rho.strata: .skeletonRho, reformat: .reformat ).
#' \item trace [logical]: Should the progress of the execution of the \code{lmm} function be displayed?
#' \item tranform.sigma, tranform.k, tranform.rho: transformation used to compute the confidence intervals/p-values for the variance and correlation parameters. See the detail section of the coef function for more information.
#' Used by \code{lmm}, \code{anova} and \code{confint}.
#' \item type.information [character]: Should the expected or observed information (\code{"expected"} or \code{"observed"}) be used to perform statistical inference? Used by \code{lmm}, \code{anova} and \code{confint}.
#' }
#'
#' @return A list containing the default options.
#'
#' @keywords utilities


 
## * LMMstar.options (code)
#' @export
LMMstar.options <- function(..., reinitialise = FALSE){
    
    ## print("start")
    ## print(sys.calls())
    ## print("stop")
    
    ## ** load current options (if necessary)
    if(reinitialise == FALSE){
        if(!is.null(names(args))){
            object <- get(".LMMstar-options", envir = LMMstar.env)
        }else{
            object <- try(get(".LMMstar-options", envir = LMMstar.env))
        }
    }

    ## ** define default options (if necessary)
    if(reinitialise || inherits(object,"try-error")){
        default <- list(adj.method = c(stats::p.adjust.methods,"single-step", "Westfall", "Shaffer", "free", "single-step2"),
                    backtransform.confint = TRUE,
                    columns.anova = c("estimate","se","df","lower","upper","p.value",""),
                    columns.confint = c("estimate","lower","upper"),
                    columns.summarize = c("observed","missing","pc.missing","mean","sd","min","q1","median","q3","max"),
                    columns.summary = c("estimate","se","df","lower","upper","p.value",""),
                    df = TRUE,
                    drop.X = TRUE,
                    effects = "mean",
                    min.df = 1, ## smallest degree-of-freedom to be used when performing Satterthwaite approximation
                    method.fit = "REML",
                    method.numDeriv = "simple",
                    n.sampleCopula = 1e5,
                    optimizer = "FS",
                    param.optimizer = c(n.iter = 100, tol.score = 1e-4, tol.param = 1e-5, n.backtracking = 10, init.cor = 1),
                    pool.method = c("average","pool.se","pool.gls","pool.gls1","pool.rubin","p.rejection"),
                    precompute.moments = TRUE,
                    sep = c(lp = ":", ## (.vcov.matrix.lmm) separator between the linear predictor when aggregated across repetitions
                            k.cov = ".", ## (.skeletonK) separator between the letter k and the covariate levels, e.g. k?2.1 
                            k.strata = ":", ## (.skeletonK) separtor between the covariate level(s), e.g. k.2?1
                            pattern = ":", ## (.findUpatterns) separtor between the index of the variance and index of the correlation pattern when forming the overall pattern name, e.g. 1?1 
                            Gpattern.var = ":", ## separtor between variable when naming pattern groups
                            Gpattern.level = ",", ## separtor between variable levels when naming pattern groups
                            rho.name = ".", ## (.skeletonRho) separator between the covariate levels, e.g. rho(1?1)
                            rho.strata = ":", ## (.skeletonRho) separator between the covariate levels and the strata, e.g. rho(1.1)?:1
                            reformat = "_"),
                    trace = FALSE,
                    transform.sigma = "log",
                    transform.k = "log",
                    transform.rho = "atanh",
                    type.information = "observed")
    }

    if (reinitialise == TRUE) {
        ## ** reinitialize
        assign(".LMMstar-options", value = default, envir = LMMstar.env)
        return(invisible(get(".LMMstar-options", envir = LMMstar.env)))
    }else{
    
        args <- list(...)
        if(inherits(object,"try-error")){
            object <- default
        }

        if(length(args)==0){

            ## ** read all
            return(object)

        }else if (is.null(names(args))) {

            ## ** read some
            args <- unlist(args)
            if(any(args %in% names(object) == FALSE)){
                stop("Incorrect element selected: \"",paste0(args[args %in% names(object) == FALSE], collapse = "\" \""),"\"\n",
                     "Available elements: \"",paste0(setdiff(names(object),args), collapse = "\" \""),"\"\n")
            }
            return(object[args])

        }else{

            ## ** write
            if(any(names(args) %in% names(object) == FALSE)){
                stop("Incorrect element selected: \"",paste0(names(args)[names(args) %in% names(object) == FALSE], collapse = "\" \""),"\"\n",
                     "Available elements: \"",paste0(setdiff(names(object),names(args)), collapse = "\" \""),"\"\n")
            }

            ok.column <- c("estimate","se","df","lower","upper","null","statistic","p.value")
            if("columns.anova" %in% names(args) && any(args$columns.confint %in% ok.column == FALSE)){
                stop("Argument \'columns.anova\' must be a character vector with values among \"",paste(ok.column, collapse = "\" \""),"\". \n")
            }
            if("columns.confint" %in% names(args) && any(args$columns.confint %in% ok.column == FALSE)){
                stop("Argument \'columns.confint\' must be a character vector with values among \"",paste(ok.column, collapse = "\" \""),"\". \n")
            }
            if("columns.summary" %in% names(args) && any(args$columns.summary %in% c(ok.column,"") == FALSE)){
                stop("Argument \'columns.summary\' must be a character vector with values among \"",paste(c(ok.column,""), collapse = "\" \""),"\". \n")
            }
            ok.column2 <- c("observed","missing","pc.missing",
                            "mean","mean.lower","mean.upper","predict.lower","predict.upper",
                            "sd","sd.lower","sd.upper", "sd0","sd0.lower","sd0.upper",
                            "min","q1","median","q3","median.upper","median.lower","IQR","max",
                            "correlation")
            if("columns.summarize" %in% names(args) && any(args$columns.summarize %in% c(ok.column2,"") == FALSE)){
                stop("Argument \'columns.summarize\' must be a character vector with values among \"",paste(c(ok.column2,""), collapse = "\" \""),"\". \n")
            }
            if("df" %in% names(args) && !is.logical(args$df)){
                stop("Argument \'df\' must be of type logical. \n")
            }
            if("drop.X" %in% names(args) && !is.logical(args$drop.X)){
                stop("Argument \'drop.X\' must be of type logical. \n")
            }
            if("optimizer" %in% names(args)){
                optimx.method <- c("BFGS", "CG", "Nelder-Mead", "nlminb", "bobyqa")
                args$optimizer <- match.arg(args$optimizer, c("FS",optimx.method)) ## FS = fisher scoring
                if(args$optimizer %in% optimx.method){
                    requireNamespace("optimx")
                }
            }
            if("param.optimizer" %in% names(args)){
                if(!is.vector(args$param.optimizer)){
                    stop("Argument \'param.optimizer\' should be a vector. \n")
                }
                valid.args <- names(object$param.optimizer)
                if(is.null(args$param.optimizer) || any(names(args$param.optimizer) %in% valid.args == FALSE)){
                    stop("Argument \'param.optimizer\' must contain elements named \"",paste(valid.args,collapse ="\" \""),"\". \n")
                }
                if(("n.iter" %in% names(args$param.optimizer)) && (args$param.optimizer["n.iter"]<=0)){
                    stop("Element \"n.iter\" in argument \'param.optimizer\' should be strictly positive. \n")
                }
                if(("n.backtracking" %in% names(args$param.optimizer)) && (args$param.optimizer["n.backtracking"]<=0)){
                    stop("Element \"n.backtracking\" in argument \'param.optimizer\' should be strictly positive. \n")
                }
                if(("tol.score" %in% names(args$param.optimizer)) && (args$param.optimizer["tol.score"]<=0)){
                    stop("Element \"tol.score\" in argument \'param.optimizer\' should be strictly positive. \n")
                }
                if(("tol.param" %in% names(args$param.optimizer)) && (args$param.optimizer["tol.param"]<=0)){
                    stop("Element \"tol.param\" in argument \'param.optimizer\' should be strictly positive. \n")
                }                
                if(("init.cor" %in% names(args$param.optimizer))){
                    if(is.numeric(init.cor)){
                        if(args$param.optimizer["init.cor"] %in% 1:2 == FALSE){
                            stop("Element \"init.cor\" in argument \'param.optimizer\' should either be 1 (i.e. \"average\") or 2 (i.e. \"overall\"). \n")
                        }
                    }else{
                        if(is.factor(init.cor)){
                            args$param.optimizer["init.cor"] <- as.character(args$param.optimizer["init.cor"])
                        }
                        if(args$param.optimizer["init.cor"] %in% c("average","overall") == FALSE){
                            stop("Element \"init.cor\" in argument \'param.optimizer\' should either be \"average\" or \"overall\". \n")
                        }
                        args$param.optimizer["init.cor"] <- as.numeric(factor(args$param.optimizer["init.cor"], levels = c("average","overall")))
                    }
                }
                param.optimizer.save <- args$param.optimizer
                args$param.optimizer <- get(".LMMstar-options", envir = LMMstar.env)$param.optimizer
                args$param.optimizer[names(param.optimizer.save)] <- param.optimizer.save
            }
            if("method.fit" %in% names(args)){
                args$method.fit <- match.arg(args$method.fit, c("ML","REML"))
            }
            if("method.numDeriv" %in% names(args)){
                args$method.numDeriv <- match.arg(args$method.numDeriv, c("simple","Richardson","complex"))
            }
            if("sep" %in% names(args)){
                sep.save <- args$sep
                check <- match.arg(sort(names(sep.save)), c("lp","kname.cov","kname.strata","pattern","pattern.var","pattern.level","rho.name","reformat"), several.ok = TRUE)
                args$sep <- object$sep
                args$sep[names(sep.save)] <- sep.save              
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
        }
    
    }
  
  
}


##----------------------------------------------------------------------
### LMMstar.options.R ends here
