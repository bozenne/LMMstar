### aggregate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul 28 2024 (19:14) 
## Version: 
## Last-Updated: okt 17 2025 (17:31) 
##           By: Brice Ozenne
##     Update #: 424
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * aggregate
##' @title Pool Estimates
##' @description Pool estimates from one or several linear mixed models. \code{confint} is the prefered method as this is meant for internal use.
##' 
##' @param x output of \code{anova.lmm}, \code{rbind.Wald_lmm}, or \code{mlmm}.
##' @param method [character] method(s) used to pool the estimates
##' @param rhs [logical or numeric vector] either the right-hand side of the null hypotheses to be tested
##' or whether the value of the summary statistic under the null hypothesis should be computed - mostly relevant when \code{method="p.rejection"} where this can be time consuming.
##' @param columns [character] column(s) to be exported
##' @param qt [numeric or character] critical quantile to be considered when evaluating the proportion of rejected hypotheses.
##' Can also be the name of a multiple testing method from which the quantile should be computed.
##' @param level [numeric, 0-1] nominal coverage of the confidence intervals.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the summary statistic. Otherwise a normal distribution is used.
##' @param tol [numeric, >0] threshold below which a pseudo-inverse is used when inverted a matrix, i.e., the eigen-values are truncated.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @seealso \code{\link{confint.Wald_lmm}}, \code{\link{confint.rbindWald_lmm}}, \code{\link{confint.mlmm}}
##'
##' @details Use a first order delta method to evaluate the standard error.
##' When \code{method="p.rejection"} the p-value and distribution under the null is evaluated by simulating data using a Gaussian or t-copula, which can be time consumming.
##' The number of sample is controlled by the argument \code{n.sampleCopula} in \code{\link{LMMstar.options}}.


## * aggregate.lmm
##' @export
aggregate.Wald_lmm <- function(x, method, rhs = NULL, columns = NULL, qt = NULL, level = 0.95, df = NULL, tol = 1e-10, ...){

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    pool.method <- options$pool.method

    ## *** method
    if(any(method %in% pool.method == FALSE)){
        stop("Incorrect value(s) \"",paste(method[method %in% pool.method == FALSE], collapse = "\" \""),"\" for argument \'method\'. \n",
             "Valid values: \"",paste(setdiff(pool.method, method), collapse = "\" \""),"\"\n")        
    }
    if(is.null(names(method))){
        names(method) <- method
    }

    ## *** columns
    valid.columns <-  c("estimate", "se", "df", "lower", "upper", "quantile", "null", "statistic", "p.value",
                        "type", "n.param", "transformed", "tobacktransform", "model", "name", "term", "method")
    save.columns <- columns
    if(identical(columns,"all")){
        columns <- valid.columns
    }else if(!is.null(columns)){
        columns <- tolower(columns)
        if(any(columns %in% valid.columns == FALSE)){
            stop("Incorrect value(s) \"",paste(columns[columns %in% valid.columns == FALSE], collapse = "\" \""),"\" for argument \'columns\'. \n",
                 "Valid values: \"",paste(setdiff(valid.columns, columns), collapse = "\" \""),"\"\n")
        }
        if(!is.null(names(columns)) && all(names(columns)=="add")){
            columns <- union(setdiff(options$columns.anova,""), unname(columns))
        }
        if(!is.null(names(columns)) && all(names(columns)=="remove")){
            columns <- setdiff(setdiff(options$columns.anova,""), unname(columns))
        }
    }else{
        columns <- setdiff(options$columns.anova,"")
    }

    ## *** level
    if(all(c("se", "df", "lower", "upper", "quantile", "statistic", "p.value") %in% columns == FALSE)){
        level <- NA
    }

    ## *** df
    if(all(c("df", "lower", "upper", "quantile", "p.value") %in% columns == FALSE)){
        df <- FALSE
    }else if(is.null(df)){
        df <- x$args$df
    }

    ## *** null
    ## can either be the null hypotheses to be tested or logical: conformity to this requirement is tested in .aggregate()
    if(is.null(rhs)){
        rhs <- any(c("null","statistic","p.value") %in% columns)
    }

    ## ** extract information from object

    ## *** estimates
    tableUni <- stats::model.tables(x, method = "none", columns = c("type","n.param","name","term","transformed","tobacktransform",
                                                                    "estimate","se","statistic","df","null"),
                                    options = options)
    if(all(tableUni$type=="rho")){
        transform.rho <- "none"
    }else if(all(tableUni$type!="rho")){
        transform.rho <- NULL
    }else{
        stop("Cannot pool linear combinations of parameters when correlation and non-correlation parameters are involved. \n")
    }

    ## *** variance
    if(any(method %in% c("pool.se","pool.gls","pool.gls1","pool.rubin")) || ((!is.na(level) || identical(as.logical(rhs),TRUE)) && "p.rejection" %in% method)){
        if(any(method %in% c("pool.se","pool.gls","pool.gls1","p.rejection")) && !is.na(level)){
            args.effects <- c("Wald","gradient")
        }else{
            args.effects <- "Wald"
        }
        Wald.Sigma <- stats::vcov(x, effects = args.effects, transform.rho = transform.rho, transform.names = FALSE, options = options)
        Wald.dSigma <- attr(Wald.Sigma,"gradient")
        attr(Wald.dSigma,"gradient") <- NULL        
    }else{
        Wald.Sigma <- NULL
        Wald.dSigma <- NULL
    }
    if(!is.na(level) && any(method %in% c("average", "pool.se","pool.gls","pool.gls1","p.rejection"))){
        contrast <- stats::model.tables(x, effects = "contrast", transform.names = FALSE, simplify = FALSE, options = options)[[1]]
        All.Sigma <- stats::vcov(x, effects = list("all",c("all","gradient"))[[df+1]], transform.names = FALSE, transform.rho = transform.rho, options = options)
        All.dSigma <- attr(All.Sigma,"gradient")
        attr(All.Sigma,"gradient") <- NULL
    }else{
        All.Sigma <- NULL
        All.dSigma <- NULL
    }

    ## *** critical quantile
    if("p.rejection" %in% method){
        if(!is.null(qt)){
            if(is.character(qt) && length(qt)!=1){
                stop("Argument \'qt\' should be have length 1 when character. \n")
            }else if(is.numeric(qt) && (length(qt)!=1 && length(qt)!=NROW(x$univariate))){
                stop("Argument \'qt\' should either have length 1 or the number of contrasts (here ",NROW(x$univariate),") when numeric. \n")
            }
            if(!is.numeric(qt) && (!is.character(qt) || (qt %in% c("none","bonferroni","single-step","single-step2")==FALSE))){
                stop("Argument \'qt\' should either be numeric or one of \"none\", \"bonferroni\", \"single-step\", \"single-step2\". \n")
            }
            if(is.numeric(qt)){
                if(length(qt)==1){
                    critical.threshold <- rep(qt, NROW(x$univariate))
                }else{
                    critical.threshold <- qt
                }
            }else{
                critical.threshold <- stats::confint(x, method = qt, columns = "quantile", options = options)$quantile
            }
        }else{
            critical.threshold <- stats::confint(x, columns = "quantile", options = options)$quantile
        }
    }else{
        critical.threshold <- NULL
    }

    ## ** pool
    ## null argument for method="p.rejection": evaluating the null and the p-value can be time consuming (done by simulation)
    outPool <- .aggregate(table = tableUni, contrast = contrast, null = rhs, level = level, df = df, threshold = critical.threshold, method = method,
                          Wald.Sigma = Wald.Sigma, Wald.dSigma = Wald.dSigma,
                          All.Sigma = All.Sigma, All.dSigma = All.dSigma,
                          tol = tol, options = options)
    pool.contrast <- attr(outPool,"contrast")
    pool.gradient <- attr(outPool,"gradient")
    
    ## ** export
    out <- outPool[,columns,drop=FALSE]
    attr(out,"contrast") <- pool.contrast
    attr(out,"gradient") <- pool.gradient
    return(out)

}

## * aggregate.rbindWald_lmm
##' @export
aggregate.rbindWald_lmm <- aggregate.Wald_lmm

## * aggregate.rbindWald_lmm
##' @export
aggregate.mlmm <- aggregate.Wald_lmm

## * .pool
.aggregate <- function(table, contrast, null, level, df, threshold, method,
                       Wald.Sigma, Wald.dSigma, All.Sigma, All.dSigma,
                       tol, options){

    ## ** normalize user input
    pool.method <- options$pool.method
    adj.method <- options$adj.method
    n.sample <- options$n.sampleCopula

    ## *** null
    ## evaluate or not the null hypothesis, test statistic, and
    if(length(null)==1 && is.logical(null)){
        hypo.test <- null
    }else if(is.null(null) || all(is.na(null))){
        hypo.test <- FALSE
    }else if(length(null) %in% c(1,NROW(table))){
        hypo.test <- TRUE
        table$null <- null ## update the null hypothesis
    }else{
        stop("Incorrect length for argument \'null\': should have length 1 or ",NROW(table)," (number of tests). \n")
    }
    
    ## *** level
    alpha <- 1-level
    
    ## *** df
    if(is.na(level)){
        df <- FALSE
    }

    ## ** prepare output
    n.test <- NROW(table)
    name.test <- rownames(table)
    if(!is.null(names(method))){    
        pool.name <- stats::setNames(names(method),method)
    }else{
        pool.name <- stats::setNames(method,method)
    }
    out <- as.data.frame(matrix(NA, nrow = length(pool.name), ncol = NCOL(table)+1,
                                dimnames = list(pool.name,  c(colnames(table),"method"))))
    out$model <- "all"
    out$type <- ifelse(length(unique(table$type))==1, table$type[1], "all")
    out$name <- ifelse(length(unique(table$name)), length(unique(table$name)), NA)
    out$term <- ifelse(length(unique(table$term)), length(unique(table$term)), NA)
    out$method <- method
    out$n.param <- sum(table$n.param)
    out$transformed <- ifelse(length(unique(table$type))==1, table$transformed[1], NA)
    out$tobacktransform <- FALSE

    ## ** point estimate
    if("p.rejection" %in% method){
        out[pool.name["p.rejection"],"estimate"] <- 1 - mean(stats::pt(threshold - table$statistic, df = table$df) - stats::pt(-threshold - table$statistic, df = table$df))
    }
    
    if(any(method != "p.rejection")){
        pool.contrast <- matrix(NA, nrow = sum(method != "p.rejection"), ncol = n.test, dimnames = list(setdiff(method, "p.rejection"),name.test))

        if("average" %in% method){
            pool.contrast["average",] <- 1/n.test
        }
        if("pool.se" %in% method){
            se.contrast <- 1/diag(Wald.Sigma)
            pool.contrast["pool.se",] <- se.contrast/sum(se.contrast)
        }
        if("pool.gls" %in% method || "pool.gls1" %in% method){
            if(is.invertible(Wald.Sigma, cov2cor = FALSE)){
                Wald.SigmaM1 <- solve(Wald.Sigma)
                gls.contrast  <- rowSums(Wald.SigmaM1)/sum(Wald.SigmaM1)
            }else{ ## use generalized inverse, same as MASS::ginv
                Wald.Sigma.eigen <- eigen(Wald.Sigma, symmetric = TRUE)
                index.subset <- which(abs(Wald.Sigma.eigen$values) > tol)

                ## truncate eigen value decomposition
                if(length(index.subset)==0){
                    stop("All eigenvalues of the variance-covariance matrix are close to 0 (<",tol,"). \n",sep="")
                }
                attr(out,"message.estimate") <-  paste(length(Wald.Sigma.eigen$values)-length(index.subset), " eigenvalues have been removed in the pseudo inverse")
                lambda.k <- Wald.Sigma.eigen$values[index.subset]
                q.kj <- Wald.Sigma.eigen$vector[,index.subset,drop=FALSE]

                ## (Fast) evaluation of weigths
                if(is.na(level)){
                    qbar.k <- colSums(q.kj)
                    w.k <- qbar.k^2/lambda.k
                    gls.contrast <- rowSums(sweep(q.kj, FUN = "*", MARGIN = 2, STATS = w.k/qbar.k))/sum(w.k)
                }
                ## (Slow but explicit) evaluation of the weights
                if(!is.na(level)){
                    Wald.SigmaM1 <- q.kj %*% diag(1/lambda.k) %*% t(q.kj)
                    gls.contrast  <- rowSums(Wald.SigmaM1)/sum(Wald.SigmaM1)
                }
            }

            if("pool.gls" %in% method){
                pool.contrast["pool.gls",] <- gls.contrast
            }
            if("pool.gls1" %in% method){ ## ensure no weight greater than 1
                index.contrastMAX <- which.max(abs(gls.contrast))
                xi.contrastMax <- sign(gls.contrast[index.contrastMAX])
                value.contrastMax <- gls.contrast[index.contrastMAX]

                kappaPwmax <- (1-n.test*xi.contrastMax*value.contrastMax)/(1-xi.contrastMax*n.test)
                pool.contrast["pool.gls1",] <- (1-1/kappaPwmax)/n.test + gls.contrast/kappaPwmax
            }
        }
        if("pool.rubin" %in% method){
            pool.contrast["pool.rubin",] <- 1/n.test
        }
        out[pool.name[rownames(pool.contrast)],"estimate"] <- pool.contrast %*% table$estimate
    }

    ## ** variance
    if(!is.na(level) && "pool.rubin" %in% method){
        pool.U <- mean(diag(Wald.Sigma))
        pool.B <- sum((table$estimate - mean(table$estimate))^2)/(n.test-1)
        out[pool.name["pool.rubin"],"se"] <- sqrt(pool.U + (1 + 1/n.test) * pool.B)
    }
    if(!is.na(level) && any(method != "pool.rubin")){

        grad.pool <- matrix(0, nrow = sum(method!="pool.rubin"), ncol = NCOL(All.Sigma), dimnames = list(setdiff(method,"pool.rubin"), colnames(contrast)))

        if("average" %in% method){
            grad.pool["average",] <- pool.contrast[rownames(pool.contrast)=="average",,drop=FALSE] %*% contrast
        }
        if("pool.se" %in% method){
            grad.pool.seA <- pool.contrast[rownames(pool.contrast) == "pool.se",] %*% contrast
            grad.pool.seB <- - table$estimate %*% apply(Wald.dSigma, MARGIN = 3, FUN = function(iM){diag(iM)/(table$se^4) })/sum(se.contrast)
            grad.pool.seC <- out[pool.name["pool.se"],"estimate"]/sum(se.contrast) * apply(Wald.dSigma, MARGIN = 3, FUN = function(iM){ sum(diag(iM)/table$se^4) })
            grad.pool["pool.se",] <- grad.pool.seA[1,] + grad.pool.seB[1,] + grad.pool.seC
        }

        if("pool.gls" %in% method || "pool.gls1" %in% method){

            if(is.invertible(Wald.Sigma, cov2cor = FALSE)){
                Wald.dSigmaM1 <- - apply(Wald.dSigma, MARGIN = 3, FUN = function(iM){rowSums(Wald.SigmaM1 %*% iM %*% Wald.SigmaM1)})
            }else{
                ## The Differentiation of Pseudo-Inverses and Nonlinear Least Squares Problems Whose Variables Separate.
                ## Author(s): G. H. Golub and V. Pereyra. Source: SIAM Journal on Numerical Analysis, Vol. 10, No. 2 (Apr., 1973), pp. 413-432 (formula 4.12)
                ## with simplification for symmetric matrix
                Wald.SigmaM1.SigmaM1 <- Wald.SigmaM1 %*% Wald.SigmaM1
                Wald.Im.Sigma.SigmaM1 <- diag(1, n.test) - Wald.Sigma %*% Wald.SigmaM1
                Wald.Im.SigmaM1.Sigma <- diag(1, n.test) - Wald.SigmaM1 %*% Wald.Sigma
                
                Wald.dSigmaM1 <- apply(Wald.dSigma, MARGIN = 3, FUN = function(iM){
                    rowSums( - Wald.SigmaM1 %*% iM %*% Wald.SigmaM1 + Wald.SigmaM1.SigmaM1 %*% iM %*% Wald.Im.Sigma.SigmaM1 + Wald.Im.SigmaM1.Sigma %*% iM %*% Wald.SigmaM1.SigmaM1) 
                })
            }
            ## d\psi/d\theta = d (w/sum(x) \beta) /d\theta = w/sum(w) d\beta/d\theta +  (dw/d\theta)/sum(x) \beta +  w/(dsum(x)/d\theta) \beta
            ##                                             = w/sum(w) d\beta/d\theta +  (dw/d\theta)/sum(x) \beta -  sum(dx/d\theta) \gamma/sum(x)
            grad.pool.glsA <- pool.contrast[rownames(pool.contrast) == "pool.gls",] %*% contrast
            grad.pool.glsB <- table$estimate %*% Wald.dSigmaM1 /sum(Wald.SigmaM1)
            grad.pool.glsC <- - out[pool.name["pool.gls"],"estimate"]/sum(Wald.SigmaM1) * colSums(Wald.dSigmaM1)
            grad.pool.gls <- grad.pool.glsA[1,] + grad.pool.glsB[1,] + grad.pool.glsC                

            if("pool.gls" %in% method){
                grad.pool["pool.gls",] <- grad.pool.gls
            }
            if("pool.gls1" %in% method){ ## ensure no weight greater than 1
                if(kappaPwmax > 1){
                    ## condense previous derivatives
                    grad.contrast.gls <- Wald.dSigmaM1/sum(Wald.SigmaM1) - cbind(rowSums(Wald.SigmaM1)) %*% rbind(colSums(Wald.dSigmaM1))/sum(Wald.SigmaM1)^2
                    ## range(table$estimate %*% grad.contrast.gls - grad.pool.glsB - grad.pool.glsC)

                    ## d\psi/d\theta = d (w\beta) /d\theta = w d\beta /d\theta + \beta dw/d\beta
                    ## where w =  (1-1/kappaPwmax)/n.test + gls.contrast/kappaPwmax i.e. dw =  d(gls.contrast)/kappaPwmax - gls.contrast d(kappaPwmax)/kappaPwmax^2 + d(kappaPwmax)/(kappaPwmax^2*n.test)
                    grad.pool.gls1A <- pool.contrast[rownames(pool.contrast) == "pool.gls1",] %*% contrast
                    grad.kappaPwmax <- -n.test*xi.contrastMax*grad.contrast.gls[index.contrastMAX,]/(1-xi.contrastMax*n.test)
                    grad.pool.gls1B <-  table$estimate %*%  (grad.contrast.gls/kappaPwmax - cbind(gls.contrast) %*% rbind(grad.kappaPwmax)/kappaPwmax^2)
                    grad.pool.gls1C <- sum(table$estimate) * grad.kappaPwmax/(kappaPwmax^2*n.test)

                    grad.pool["pool.gls1",] <- grad.pool.gls1A[1,] + grad.pool.gls1B[1,] + grad.pool.gls1C                
                }else{
                    grad.pool["pool.gls1",] <- grad.pool.gls
                }
            }
        }

        if("p.rejection" %in% method){
            ## Approximation: no variability of the degrees-of-freedom nor critical threshold
            grad.statistic <- contrast / table$se - (table$estimate - table$null) * apply(Wald.dSigma, MARGIN = 3, FUN = diag)  / (2 * table$se^3)
            partial.integral <- stats::dt(threshold - table$statistic, df = table$df) - stats::dt(-threshold - table$statistic, df = table$df)
            grad.pool[method == "p.rejection",] <- colMeans(grad.statistic * partial.integral)
        }
        out[pool.name[rownames(grad.pool)],"se"] <- sqrt(rowSums( (grad.pool %*% All.Sigma) * grad.pool))
    }

    ## ** null hypothesis
    if("p.rejection" %in% method && hypo.test){
        rho.linfct <- stats::cov2cor(Wald.Sigma)
        
        myMvd <- copula::mvdc(copula = copula::normalCopula(param=rho.linfct[lower.tri(rho.linfct)], dim = NROW(rho.linfct), dispstr = "un"),
                              margins = rep("t", NROW(rho.linfct)),
                              paramMargins = as.list(stats::setNames(table$df,rep("df",NROW(rho.linfct)))))
        sample.copula <- copula::rMvdc(n.sample, myMvd)

        null.integral <- do.call(cbind, lapply(1:n.test, function(iTest){
            stats::pt(threshold[iTest] - sample.copula[,iTest], df = table$df[iTest]) - stats::pt(-threshold[iTest] - sample.copula[,iTest], df = table$df[iTest])
        }))
        null.rejection <- 1 - rowMeans(null.integral)
        out[pool.name["p.rejection"],"null"] <- mean(null.rejection)
        out[pool.name["p.rejection"],"p.value"] <- mean(null.rejection>out[pool.name["p.rejection"],"estimate"])
    }
    if(any(method != "p.rejection") && hypo.test){
        out[setdiff(pool.name, pool.name["p.rejection"]),"null"] <- (pool.contrast %*% table$null)[,1]
    }

    ## ** degrees-of-freedom
    if(df == FALSE){

        out$df <- Inf

    }else{
        if("pool.rubin" %in% method){
            ## MICE's approach: https://stefvanbuuren.name/fimd/sec-whyandwhen.html
            pool.lambda <- (1+1/n.test)*pool.B / out[pool.name["pool.rubin"],"se"]^2 ## formula 2.24
            pool.nu_old <- (n.test-1)/pool.lambda^2 ## formula 2.30 (Rubin 1987b eq 3.1.6)
            pool.nu_obs <- (table$df+1)/(table$df+3)*table$df*(1-pool.lambda) ## formula 2.31 (Barnard and Rubin (1999) )
            out[pool.name["pool.rubin"],"df"] <- mean(pool.nu_old*pool.nu_obs/(pool.nu_old+pool.nu_obs)) ## formula 2.32
        }

        if(any(c("average","pool.se","pool.gls","pool.gls1","p.rejection") %in% method)){
            out[pool.name[rownames(grad.pool)],"df"] <- .df_contrast(contrast = grad.pool[,rownames(All.dSigma),drop=FALSE], vcov.param = All.Sigma, dVcov.param = All.dSigma)            
        }
    }

    ## ** export
    out$statistic <- (out$estimate - out$null)/out$se
    out$lower <- out$estimate + out$se * stats::qt(alpha/2, df = out$df)
    out$upper <- out$estimate + out$se * stats::qt(1-alpha/2, df = out$df)
    out$quantile <- stats::qt(1-alpha/2, df = out$df)
    if(any(method != "p.rejection")){
        out[method != "p.rejection","p.value"] = 2*(1-stats::pt(abs(out[method != "p.rejection","statistic"]), df = out[method != "p.rejection","df"]))
        attr(out,"contrast") <- pool.contrast
    }
    if(!is.na(level) && any(method != "pool.rubin")){
        attr(out,"gradient") <- grad.pool
    }
    return(out)
}


##----------------------------------------------------------------------
### aggregate.R ends here
