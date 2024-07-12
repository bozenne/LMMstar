### proportion.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 15 2022 (14:09) 
## Version: 
## Last-Updated: jul 11 2024 (15:50) 
##           By: Brice Ozenne
##     Update #: 82
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * proportion.mlmm
proportion.mlmm <- function(object, index, name.method, method, qt = NULL, null, ci, df, alpha){

    options <- LMMstar.options()
    
    ## ** normalize user input
    if(!is.null(qt)){
        if(s.numeric(qt) || length(qt)!=1){
            stop("Argument \'qt\' should be a numeric value (with length 1). \n")
        }
        critical.threshold <- qt
    }else{
        critical.threshold <- NULL
    }

    ## ** name for the pooled estimator
    if(!is.null(name.method)){
        pool.name <- name.method
    }else{
        pool.name <- poolName.mlmm(object, index = index, method = "p.rejection")
    }

    ## ** estimate proportion
    object.ci <- confint(object, method = method, columns = c("estimate","se","df","lower","upper","statistic","null","p.value"))
    if(is.null(critical.threshold)){
        if(!is.null(attr(object.ci, "quantile"))){
            critical.threshold <- attr(object.ci, "quantile")
            if(is.list(critical.threshold)){
                critical.threshold <- critical.threshold[[object$univariate[index,"test"][1]]]
            }
        }else{
            stop("Unknown critical threshold: consider specifying argument \'qt\' or changing argument \'method\'. \n")
        }
    }

    integral <- stats::pt(critical.threshold - object.ci$statistic, df = object.ci$df) - stats::pt(-critical.threshold - object.ci$statistic, df = object.ci$df)
    estimate <- 1 - mean(integral)

    ## ** null
    if(is.null(null)){
        method <- attr(object.ci,"method")
        rho.linfct <- stats::cov2cor(vcov(object))
        n.test <- NROW(object.ci)

        n.sample <- options$n.sampleCopula
        myMvd <- copula::mvdc(copula = copula::normalCopula(param=rho.linfct[lower.tri(rho.linfct)], dim = NROW(rho.linfct), dispstr = "un"),
                              margins = rep("t", NROW(rho.linfct)),
                              paramMargins = as.list(stats::setNames(object.ci$df,rep("df",NROW(rho.linfct)))))
        sample.copula <- copula::rMvdc(n.sample, myMvd)

        null.integral <- do.call(cbind, lapply(1:n.test, function(iTest){
            stats::pt(critical.threshold[iTest] - sample.copula[,iTest], df = object.ci$df[iTest]) - stats::pt(-critical.threshold[iTest] - sample.copula[,iTest], df = object.ci$df[iTest])
        }))
        null <- mean(1 - rowMeans(null.integral))
    }

    ## ** variance
    if(!ci){

        integral.se  <- NA

    }else{

        browser()

    }

    ## ** degree of freedom
    if(!ci || df == FALSE){

        integral.df <- Inf

    }else{

        browser()

    }

    ## ** post process
    alpha <- 1-attr(object.ci, "level")

    out <- data.frame(estimate = as.double(estimate),
                      se = integral.se,
                      df = integral.df,
                      stratistic = NA,
                      lower = NA,
                      upper = NA,
                      null = null,
                      p.value = NA)
    out$statistic <- as.double((out$estimate-out$null)/out$se)
    out$lower <- out$estimate + out$se * stats::qt(alpha/2, df = out$df)
    out$upper <- out$estimate + out$se * stats::qt(1-alpha/2, df = out$df)
    out$p.value = 2*(1-stats::pt(abs(out$statistic), df = out$df))

    rownames(out) <- pool.name
    attr(out,"method") <- method
    attr(out,"quantile") <- critical.threshold
    return(out)
}

##----------------------------------------------------------------------
### proportion.R ends here
