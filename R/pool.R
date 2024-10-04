### pool.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul 28 2024 (19:14) 
## Version: 
## Last-Updated: sep 30 2024 (14:26) 
##           By: Brice Ozenne
##     Update #: 203
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * pool.rbindWald_lmm
pool.rbindWald_lmm <- function(object, method, qt = NULL,
                               null, level, df, tol = 1e-10, options){

    if(is.null(options)){
        options <- LMMstar.options()
    }
    
    pool.method <- options$pool.method
    adj.method <- options$adj.method
    n.sample <- options$n.sampleCopula
        
    ## ** normalize user input

    ## *** level
    alpha <- 1-level
    if(!is.na(level) && object$args$simplify==TRUE && any(method!="pool.rubin")){
        stop("Cannot evaluate the uncertainty of pooling estimator without the full variance-covariance matrix. \n",
             "Consider setting the argument \'simplify\' to FALSE when calling anova. \n")
    }
    if(!is.na(level) && object$args$df==FALSE && any(method %in% c("pool.se","pool.gls","pool.gls1","p.rejection"))){
        stop("Cannot evaluate the uncertainty of pooling estimator without the derivative of the variance-covariance matrix. \n",
             "Consider setting the argument \'df\' to TRUE when calling lmm and anova. \n")
    }

    ## *** qt (critical quantile)
    if("p.rejection" %in% method){
        if(!is.null(qt)){
            if(is.character(qt) && length(qt)!=1){
                stop("Argument \'qt\' should be have length 1 when character. \n")
            }else if(is.numeric(qt) && (length(qt)!=1 && length(qt)!=NROW(object$univariate))){
                stop("Argument \'qt\' should either have length 1 or the number of contrasts (here ",NROW(object$univariate),") when numeric. \n")
            }
            if(!is.numeric(qt) && (!is.character(qt) || (qt %in% c("none","bonferroni","single-step","single-step2")==FALSE))){
                stop("Argument \'qt\' should either be numeric or one of \"none\", \"bonferroni\", \"single-step\", \"single-step2\". \n")
            }
            if(is.numeric(qt)){
                critical.threshold <- qt
            }else{
                critical.threshold <- stats::confint(object, method = qt, columns = "quantile", options = options)$quantile
            }
        }else{
            critical.threshold <- stats::confint(object, columns = "quantile", options = options)$quantile
        }
    }

    ## *** df
    if(is.na(level)){
        df <- FALSE
    }

    ## ** prepare output
    if(!is.null(names(method))){
        pool.name <- stats::setNames(names(method),method)
    }else{
        pool.name <- stats::setNames(method,method)## stats::setNames(poolName.rbindWald_lmm(object, method = method),method)
    }
    out <- as.data.frame(matrix(NA, nrow = length(pool.name), ncol = NCOL(object$univariate),
                                dimnames = list(pool.name,  colnames(object$univariate))))
    out$model <- "all"
    out$type <- ifelse(length(unique(object$univariate$type))==1, object$univariate$type[1], "all")
    out$name <- method
    out$term <- "pool"
    out$n.param <- sum(object$univariate$n.param)
    out$transformed <- ifelse(length(unique(object$univariate$type))==1, object$univariate$transformed[1], NA)
    out$tobacktransform <- object$univariate$tobacktransform[1]

    ## ** extract information
    ## independence between models
    independence <- object$object$independence
    
    ## estimates
    tableUni <- stats::model.tables(object, method = "none", columns = c("estimate","se","statistic","df","null"), options = options)
    n.test <- NROW(tableUni)
    name.test <- rownames(tableUni)

    ## variance
    if(any(method %in% c("pool.se","pool.gls","pool.gls1")) || ((!is.na(level) || null) && "p.rejection" %in% method)){
        Wald.Sigma <- stats::vcov(object, effects = list("Wald",c("Wald","gradient"))[[(!is.na(level))+1]],
                                  method = "none", transform.names = FALSE, options = options)
        Wald.dSigma <- attr(Wald.Sigma,"gradient")
        attr(Wald.dSigma,"gradient") <- NULL
        
    }
    if(!is.na(level) && any(method %in% c("average", "pool.se","pool.gls","pool.gls1","p.rejection"))){
        contrast <- stats::model.tables(object, effects = "contrast", transform.names = FALSE, simplify = FALSE, options = options)
        All.Sigma <- stats::vcov(object, effects = list("all",c("all","gradient"))[[df+1]],
                                 method = "none", transform.names = FALSE, options = options)
        All.dSigma <- attr(All.Sigma,"gradient")
        attr(All.Sigma,"gradient") <- NULL
    }


    ## ** point estimate
    if("p.rejection" %in% method){
        out[pool.name["p.rejection"],"estimate"] <- 1 - mean(stats::pt(critical.threshold - tableUni$statistic, df = tableUni$df) - stats::pt(-critical.threshold - tableUni$statistic, df = tableUni$df))
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

        out[pool.name[rownames(pool.contrast)],"estimate"] <- pool.contrast %*% tableUni$estimate
    }

    ## ** variance
    if(!is.na(level) && "pool.rubin" %in% method){
        pool.U <- mean(diag(Wald.Sigma))
        pool.B <- sum((tableUni$estimate - mean(tableUni$estimate))^2)/(n.test-1)
        out[pool.name["pool.rubin"],"se"] <- sqrt(pool.U + (1 + 1/n.test) * pool.B)
    }
    if(!is.na(level) && any(method != "pool.rubin")){

        grad.pool <- matrix(0, nrow = sum(method!="pool.rubin"), ncol = NCOL(All.Sigma), dimnames = list(setdiff(method,"pool.rubin"), colnames(contrast)))

        if("average" %in% method){
            grad.pool["average",] <- pool.contrast[rownames(pool.contrast)=="average",,drop=FALSE] %*% contrast
        }
        if("pool.se" %in% method){
            grad.pool.seA <- pool.contrast[rownames(pool.contrast) == "pool.se",] %*% contrast
            grad.pool.seB <- - tableUni$estimate %*% apply(Wald.dSigma, MARGIN = 3, FUN = function(iM){diag(iM)/(tableUni$se^4) })/sum(se.contrast)
            grad.pool.seC <- out[pool.name["pool.se"],"estimate"]/sum(se.contrast) * apply(Wald.dSigma, MARGIN = 3, FUN = function(iM){ sum(diag(iM)/tableUni$se^4) })
            grad.pool["pool.se",] <- grad.pool.seA[1,] + grad.pool.seB[1,] + grad.pool.seC
        }

        if("pool.gls" %in% method || "pool.gls1" %in% method){

            if(is.invertible(Wald.Sigma, cov2cor = FALSE)){
                Wald.dSigmaM1 <- - apply(Wald.dSigma, MARGIN = 3, FUN = function(iM){rowSums(Wald.SigmaM1 %*% iM %*% Wald.SigmaM1)})
            }else{
                ## The Differentiation of Pseudo-Inverses and Nonlinear Least Squares Problems Whose Variables Separate. Author(s): G. H. Golub and V. Pereyra. Source: SIAM Journal on Numerical Analysis, Vol. 10, No. 2 (Apr., 1973), pp. 413-432 (formula 4.12)
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
            grad.pool.glsB <- tableUni$estimate %*% Wald.dSigmaM1 /sum(Wald.SigmaM1)
            grad.pool.glsC <- - out[pool.name["pool.gls"],"estimate"]/sum(Wald.SigmaM1) * colSums(Wald.dSigmaM1)
            grad.pool.gls <- grad.pool.glsA[1,] + grad.pool.glsB[1,] + grad.pool.glsC                

            if("pool.gls" %in% method){
                grad.pool["pool.gls",] <- grad.pool.gls
            }
            if("pool.gls1" %in% method){ ## ensure no weight greater than 1
                if(kappaPwmax > 1){
                    ## condense previous derivatives
                    grad.contrast.gls <- Wald.dSigmaM1/sum(Wald.SigmaM1) - cbind(rowSums(Wald.SigmaM1)) %*% rbind(colSums(Wald.dSigmaM1))/sum(Wald.SigmaM1)^2
                    ## range(tableUni$estimate %*% grad.contrast.gls - grad.pool.glsB - grad.pool.glsC)

                    ## d\psi/d\theta = d (w\beta) /d\theta = w d\beta /d\theta + \beta dw/d\beta
                    ## where w =  (1-1/kappaPwmax)/n.test + gls.contrast/kappaPwmax i.e. dw =  d(gls.contrast)/kappaPwmax - gls.contrast d(kappaPwmax)/kappaPwmax^2 + d(kappaPwmax)/(kappaPwmax^2*n.test)
                    grad.pool.gls1A <- pool.contrast[rownames(pool.contrast) == "pool.gls1",] %*% contrast
                    grad.kappaPwmax <- -n.test*xi.contrastMax*grad.contrast.gls[index.contrastMAX,]/(1-xi.contrastMax*n.test)
                    grad.pool.gls1B <-  tableUni$estimate %*%  (grad.contrast.gls/kappaPwmax - cbind(gls.contrast) %*% rbind(grad.kappaPwmax)/kappaPwmax^2)
                    grad.pool.gls1C <- sum(tableUni$estimate) * grad.kappaPwmax/(kappaPwmax^2*n.test)

                    grad.pool["pool.gls1",] <- grad.pool.gls1A[1,] + grad.pool.gls1B[1,] + grad.pool.gls1C                
                }else{
                    grad.pool["pool.gls1",] <- grad.pool.gls1A[1,]
                }
            }
        }

        if("p.rejection" %in% method){
            ## Approximation: no variability of the degrees-of-freedom nor critical threshold
            grad.statistic <- contrast / tableUni$se - (tableUni$estimate - tableUni$null) * apply(Wald.dSigma, MARGIN = 3, FUN = diag)  / (2 * tableUni$se^3)
            partial.integral <- stats::dt(critical.threshold - tableUni$statistic, df = tableUni$df) - stats::dt(-critical.threshold - tableUni$statistic, df = tableUni$df)
            grad.pool[method == "p.rejection",] <- colMeans(grad.statistic * partial.integral)
        }
        out[pool.name[rownames(grad.pool)],"se"] <- sqrt(rowSums( (grad.pool %*% All.Sigma) * grad.pool))
    }

    ## ** null hypothesis
    if(null){
        if("p.rejection" %in% method){
            rho.linfct <- stats::cov2cor(Wald.Sigma)
        
            myMvd <- copula::mvdc(copula = copula::normalCopula(param=rho.linfct[lower.tri(rho.linfct)], dim = NROW(rho.linfct), dispstr = "un"),
                                  margins = rep("t", NROW(rho.linfct)),
                                  paramMargins = as.list(stats::setNames(tableUni$df,rep("df",NROW(rho.linfct)))))
            sample.copula <- copula::rMvdc(n.sample, myMvd)

            null.integral <- do.call(cbind, lapply(1:n.test, function(iTest){
                stats::pt(critical.threshold - sample.copula[,iTest], df = tableUni$df[iTest]) - stats::pt(-critical.threshold - sample.copula[,iTest], df = tableUni$df[iTest])
            }))
            out[pool.name["p.rejection"],"null"] <- mean(1 - rowMeans(null.integral))
        }
        if(any(method != "p.rejection")){
            out[setdiff(pool.name, pool.name["p.rejection"]),"null"] <- (pool.contrast %*% tableUni$null)[,1]
        }
    }

    ## ** degrees-of-freedom
    if(df == FALSE){

        out$df <- Inf

    }else{
        if("pool.rubin" %in% method){
            ## MICE's approach: https://stefvanbuuren.name/fimd/sec-whyandwhen.html
            pool.lambda <- (1+1/n.test)*pool.B / out[pool.name["pool.rubin"],"se"]^2 ## formula 2.24
            pool.nu_old <- (n.test-1)/pool.lambda^2 ## formula 2.30 (Rubin 1987b eq 3.1.6)
            pool.nu_obs <- (tableUni$df+1)/(tableUni$df+3)*tableUni$df*(1-pool.lambda) ## formula 2.31 (Barnard and Rubin (1999) )
            out[pool.name["pool.rubin"],"df"] <- mean(pool.nu_old*pool.nu_obs/(pool.nu_old+pool.nu_obs)) ## formula 2.32
        }

        if(any(c("average","pool.se","pool.gls","pool.gls1","p.rejection") %in% method)){
            out[pool.name[rownames(grad.pool)],"df"] <- .df_contrast(contrast = grad.pool[,rownames(All.dSigma),drop=FALSE], vcov.param = All.Sigma, dVcov.param = All.dSigma)
            if(object$args$simplify!=FALSE){
                if(all(out$type=="mu")){
                    attr(out,"message.df") <- "ignoring variability of the variance parameters"
                }else{
                    attr(out,"message.df") <- "ignoring variability of some parameters"
                }
            }
            ## approximate df calculation 
        }
    }

    ## ** export
    out$statistic <- (out$estimate - out$null)/out$se
    out$lower <- out$estimate + out$se * stats::qt(alpha/2, df = out$df)
    out$upper <- out$estimate + out$se * stats::qt(1-alpha/2, df = out$df)
    out$quantile <- stats::qt(1-alpha/2, df = out$df)
    out$p.value = 2*(1-stats::pt(abs(out$statistic), df = out$df))
    if(any(method != "p.rejection")){
        attr(out,"contrast") <- pool.contrast
    }
    if(!is.na(level) && any(method != "pool.rubin")){
        attr(out,"gradient") <- grad.pool
    }
    return(out)
}

## * poolName.rbindWald_lmm
poolName.rbindWald_lmm <- function(object, method){

    ## ** extract
    table <- object$univariate

    ## ** build name
    Wald.name <- rownames(table)
    sep <- object$object$sep
    pool.splitname <- strsplit(Wald.name, split = sep, fixed = TRUE)

    if(length(Wald.name)==1){
        pool.name <- paste0("<",Wald.name,">")
    }else if(all(lengths(pool.splitname)==2)){
        ## ??
        Mpool.splitname <- do.call(rbind,pool.splitname)
        pool.splitname2 <- strsplit(Mpool.splitname[,1],split="=", fixed = TRUE)
        if(all(lengths(pool.splitname2)==2)){
            Mpool.splitname2 <- do.call(rbind,pool.splitname2)
            if(length(unique(Mpool.splitname2[,1]))==1){
                pool.name <- paste0(Mpool.splitname2[1,1],"=<",paste(Mpool.splitname2[1,2],Mpool.splitname2[NROW(Mpool.splitname2),2],sep=":"),">",sep,Mpool.splitname[1,2])                   
            }else{
                pool.name <- paste0("<",paste(Mpool.splitname[1,1],Mpool.splitname[NROW(Mpool.splitname),1],sep=":"),">",sep,Mpool.splitname[1,2])                   
            }
        }else if(length(unique(Mpool.splitname[,2]))==1){
            pool.name <- paste0("<",paste(Mpool.splitname[1,1],Mpool.splitname[NROW(Mpool.splitname),1],sep=":"),">",sep,Mpool.splitname[1,2])                   
        }else{
            pool.name <- paste0("<",paste(Wald.name[1],Wald.name[length(Wald.name)],collapse=":"),">")                   
        }
    }else if(length(Wald.name)==2){
        pool.name <- paste0("<",paste(Wald.name,collapse=", "),">")
    }else{
        ## <first test, last test>
        pool.name <- paste0("<",paste(Wald.name[1],"...",Wald.name[length(Wald.name)],sep=", "),">")
    }

    ## **  export
    return(paste0(method,pool.name))
}


##----------------------------------------------------------------------
### pool.R ends here
