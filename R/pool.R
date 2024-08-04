### pool.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul 28 2024 (19:14) 
## Version: 
## Last-Updated: Aug  4 2024 (22:04) 
##           By: Brice Ozenne
##     Update #: 74
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
                               null, level, df, tol = 1e-10){

    options <- LMMstar.options()
    pool.method <- options$pool.method
    adj.method <- options$adj.method

    ## ** normalize user input

    ## *** qt (critical quantile)
    if("p.rejection" %in% method){
        if(!is.null(qt)){
            if(length(qt)!=1){
                stop("Argument \'qt\' should be have length 1. \n")
            }
            if(!is.numeric(qt) && (!is.character(qt) || (qt %in% c("none","bonferroni","single-step","single-step2")==FALSE))){
                stop("Argument \'qt\' should either be numeric or one of \"none\", \"bonferroni\", \"single-step\", \"single-step2\". \n")
            }
            if(is.numeric(qt)){
                critical.threshold <- qt
            }else{
                critical.threshold <- confint(object, method = qt, columns = "quantile")$quantile[1]
            }
        }else{
            critical.threshold <- confint(object, columns = "quantile")$quantile[1]
        }
    }

    ## ** prepare output
    if(!is.null(names(method))){
        pool.name <- names(method)
    }else{
        pool.name <- poolName.rbindWald_lmm(object, method = method)
    }
    out <- as.data.frame(matrix(NA, nrow = length(pool.name), ncol = NCOL(object$univariate),
                                dimnames = list(pool.name,  colnames(object$univariate))))
    out$model <- "all"
    out$type <- ifelse(length(unique(object$univariate$type))==1, object$univariate$type[1], "all")
    out$name <- method
    out$term <- "pool"
    out$n.param <- sum(object$univariate$n.param)
    out$transformed <- out$transformed[1]
    out$backtransformed <- out$transformed[1]

    ## ** extract information
    ## estimates
    tableUni <- model.tables(object, method = "none", columns = c("estimate","se","statistic","df"))

    ## variance
    if(!is.na(level) || any(method %in% c("pool.se","pool.gls","pool.gls1"))){

        Sigma <- vcov(object, effects = "Wald", method = "none")
        independence <- object$object$independence

    }

    ## ** point estimate
    n.test <- NROW(tableUni)
    name.test <- rownames(tableUni)
    method.contrast <- method[method != "p.rejection"]
    pool.contrast <- matrix(NA, nrow = length(method.contrast), ncol = n.test, dimnames = list(pool.name[match(method.contrast, method)],name.test))
    
    if("average" %in% method.contrast || "pool.rubin" %in% method.contrast){
        pool.contrast[method.contrast %in% c("average","pool.rubin"),] <- 1/n.test
    }
    if("pool.se" %in% method){
        iSigma.diag <- 1/diag(Sigma)
        pool.contrast[method.contrast == "pool.se",] <- iSigma.diag/sum(iSigma.diag)
    }
    if("pool.gls" %in% method.contrast || "pool.gls1" %in% method.contrast){
        if(is.invertible(Sigma, cov2cor = FALSE)){
            SigmaM1 <- solve(Sigma)
            gls.contrast  <- rowSums(SigmaM1)/sum(SigmaM1)
            message <- NULL
        }else{
            Sigma.eigen <- eigen(Sigma, symmetric = TRUE)
            index.subset <- which(abs(Sigma.eigen$values) > tol)

            ## truncate eigen value decomposition
            if(length(index.subset)==0){
                stop("All eigenvalues of the variance-covariance matrix are close to 0 (<",tol,"). \n",sep="")
            }
            message <-  length(Sigma.eigen$values)-length(index.subset)

            lambda.k <- Sigma.eigen$values[index.subset]
            q.kj <- Sigma.eigen$vector[,index.subset,drop=FALSE]

            ##  evaluate weigths
            qbar.k <- colSums(q.kj)
            w.k <- qbar.k^2/lambda.k
            gls.contrast <- out <- rowSums(sweep(q.kj, FUN = "*", MARGIN = 2, STATS = w.k/qbar.k))/sum(w.k)
        }

        if("pool.gls" %in% method.contrast){
            pool.contrast[method.contrast == "pool.gls",] <- gls.contrast
        }
    }else{
        message <- NULL
    }       
    if("pool.gls1" %in% method.contrast){ ## ensure no weight greater than 1
        index.max <- which.max(abs(gls.contrast))
        kappaPwmax <- max(c(1,(1-n.test*gls.contrast)/(1-sign(gls.contrast)*n.test))) 
        pool.contrast[method.contrast == "pool.gls1",] <- (1-1/kappaPwmax)/n.test + gls.contrast/kappaPwmax
    }
    out$estimate[match(rownames(pool.contrast), rownames(out))] <- pool.contrast %*% tableUni$estimate

    if("p.rejection" %in% method){
        out[method == "p.rejection","estimate"] <- 1 - mean(stats::pt(critical.threshold - tableUni$statistic, df = tableUni$df) - stats::pt(-critical.threshold - tableUni$statistic, df = tableUni$df))
    }
browser()
    
    ## ** variance
    if(!is.na(level) && "average" %in% method){

        out[method == "average","se"] <- sqrt(as.double( pool.contrast[method.contrast=="average",,drop=FALSE] %*% Sigma %*% t(pool.contrast[method.contrast=="average",,drop=FALSE]) ))

    }
    if(method %in% c("pool.se","pool.gls","pool.gls1")){

        
        ## *** all model parameters
        ls.theta <- lapply(coef(object, effects = "all", ordering = "by",
                                transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE),
                           function(iCoef){
                               attr(iCoef,"transform.sigma") <- transform.sigma
                               attr(iCoef,"transform.k") <- transform.k
                               attr(iCoef,"transform.rho") <- transform.rho
                               return(iCoef)
                           })

        ## GOLD STANDARD
        ## ls.estimate <- lapply(names(object$model), FUN = function(iM){ ## iM <- "A"
        ##     lava::estimate(object$model[[iM]], f = function(p){ ## p <- theta[[iM]]
        ##         iTheta <- ls.theta
        ##         iTheta[[iM]] <- p
        ##         iSigma <- vcov(object, p = iTheta, type.information = type.information, robust = robust,
        ##                        transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)[index,index,drop=FALSE]
        ##         iCoef <- coef(object, p = iTheta, effects = "contrast")
        ##         sum(iCoef/diag(iSigma))/sum(1/diag(iSigma))
        ##     }, df = df, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
        ## })
        ## pool.se <- sqrt(sum(sapply(ls.estimate,"[[","se")^2))

        ## *** variance-covariance matrix for all model parameters
        theta.Sigma <- vcov(object, effects = "all", type.information = type.information, robust = robust,
                            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)

        ## *** derivative: d(w\beta)/d\theta = dw/d\theta*\beta + w*d\beta/d\theta = dw/d\theta*\beta + w*d\beta/d\mu (as \beta only depends on the mean parameters)
        ## d\beta/d\mu
        ls.contrast <- coef(object, type = "ls.contrast") ## list of contrast matrices to extract coefficients from each model
        theta.grad <- do.call(cbind,ls.contrast) ## combine into a single matrix
        theta.grad[,theta.grad!=0] <- theta.grad[,theta.grad!=0]*pool.contrast[1,] ## times inverse variance weights
        colnames(theta.grad) <- colnames(theta.Sigma)
        ## theta.grad %*% cbind(unlist(ls.theta))

        ## dw/d\theta
        ls.gradVcov <- lapply(names(object$model), FUN = function(iM){ ## iM <- "A"

            iGrad <- numDeriv::jacobian(func = function(x){
                iSigma <- vcov(object, p = stats::setNames(list(x), iM), type.information = type.information, robust = robust,
                               transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
                if(method=="pool.se"){
                    iOut <- (1/diag(iSigma))/(sum(1/diag(iSigma)))
                }else{
                    iOut <- poolGLS(Sigma = iSigma, method = method)
                }

                return(iOut)
            }, x = ls.theta[[iM]])

            return(iGrad)

        })
        
        theta.grad <- theta.grad + table$estimate %*% do.call(cbind,ls.gradVcov)
        pool.se <- sqrt(theta.grad %*% theta.Sigma %*% t(theta.grad))
        
    }

    if(!is.na(level) && "pool.rubin" %in% method){
        pool.U <- mean(diag(Sigma))
        pool.B <- sum((tableUni$estimate - mean(tableUni$estimate))^2)/(n.test-1)
        out[method == "pool.rubin","se"] <- sqrt(pool.U + (1 + 1/n.test) * pool.B)
    }

    if(!is.na(level) && "p.rejection" %in% method){
        ## Approximation: no variability of the degrees of freedom nor critical threshold

        contrast <- model.tables(object, effects = "contrast", simplify = FALSE)[[1]]
        vcov(object, effects = c("all","gradient"))
        object$glht[[1]]$vcov
        dstat <- contrast * tableUni$se - tableUni$estimate / (2*tableUni$se^3)

        term1 <- - coef(e.lmm1)["X1"] * attr(vcov(e.lmm1, effects = c("all","gradient")), "gradient")["X1","X1",]  / (2 * vcov(e.lmm1)["X1","X1"]^(3/2))
        term2 <- stats::setNames( (names(coef(e.lmm1, effects = "all"))=="X1") / sqrt(vcov(e.lmm1)["X1","X1"]) , names(coef(e.lmm1, effects = "all")))  


        out[method == "p.rejection","se"] <- 1 - mean(stats::pt(critical.threshold - tableUni$statistic, df = tableUni$df) - stats::pt(-critical.threshold - tableUni$statistic, df = tableUni$df))
    }

    ## ** degrees of freedom
    if(!ci || df == FALSE){

        pool.df <- Inf

    }else if(method == "pool.rubin"){

        ## MICE's approach: https://stefvanbuuren.name/fimd/sec-whyandwhen.html
        pool.lambda <- (1+1/n.test)*pool.B / out$se
        pool.nu_old <- (n.test-1)/pool.lambda^2 ## formula 2.30 (Rubin 1987b eq 3.1.6)
        pool.nu_obs <- (estimate.df+1)/(estimate.df+3)*estimate.df*(1-pool.lambda) ## formula 2.31 (Barnard and Rubin (1999) )
        pool.df <- mean(pool.nu_old*pool.nu_obs/(pool.nu_old+pool.nu_obs)) ## formula 2.32

    }else if(diff(range(table$df))<0.1){

        pool.df <- mean(table$df)

    }else if(method %in% c("average","pool.fixse")){

        ## \hat{pool} = \sum_k w_k \hat{beta}_k
        ## \sigma^2_{\hat{pool}} = \Var[\hat{pool}] = \sum_k w^2_k \Var[\hat{beta}_k] = \sum_k w^2_k \sigma^2_{\hat{beta}_k} by independence
        ## \Var[\sigma^2_{\hat{pool}}] = \sum_k w^4_k \Var[\sigma^2_{\hat{beta}_k}] by independence
        ## retrieve each denominator: \Var[\sigma^2_{\hat{beta}_k}]
        lmm.ls.contrast <- coef(object, type = "ls.contrast",
                                transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
        lmm.ls.vcov <- vcov(object, effects = "all", 
                            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE, simplify = FALSE)
        lmm.ls.dVcov <- lapply(object$model, function(iM){
            attr(vcov(iM, effects = c("all","gradient"), 
                 transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE),"gradient")
        })
        lmm.df <- mapply(iC = lmm.ls.contrast, iVcov = lmm.ls.vcov, iDvcov = lmm.ls.dVcov, FUN = function(iC,iVcov,iDvcov){
            iDf <- .df_contrast(contrast = iC, vcov.param = iVcov, dVcov.param = iDvcov, return.vcov = TRUE)
            return(iDf)
        }, SIMPLIFY = FALSE)

        Vdenum.df <- sapply(lmm.df, attr, "denum")

        if(pool.contrast^4 %*% Vdenum.df <= 0){
            pool.df <- min(table$df)
        }else{
            pool.df <- max(as.double(2*pool.se^4 / (pool.contrast^4 %*% Vdenum.df)), min(table$df))
        }

    }else if(method %in% c("pool.se","pool.gls","pool.gls1")){ 

        ## ADD-HOC APPROXIMATION (ignores correlation in dVcov across models)
        theta.dVcov <- array(0, dim = rep(length(theta.grad),3), dimnames = list(colnames(theta.grad),colnames(theta.grad),colnames(theta.grad)))
        for(iBy in names(object$model)){ ## iBy <- "A"
            iDvcov <- attr(vcov(object$model[[iBy]], effects = c("all","gradient"),
                                transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE),"gradient")
            theta.dVcov[attr(theta.Sigma,"by")==iBy,attr(theta.Sigma,"by")==iBy,attr(theta.Sigma,"by")==iBy] <- iDvcov
        }
        pool.df <- .df_contrast(contrast = theta.grad, vcov.param = theta.Sigma, dVcov.param = theta.dVcov)
    }

    ## ** export
    if(length(unique(out$parameter))!=1){
        out$parameter <- NA
    }
    out$estimate <- as.double(pool.contrast %*% table$estimate)
    out$se <- as.double(pool.se)
    out$df <- as.double(pool.df)
    out$null <- as.double(pool.contrast %*% table$null)
    out$statistic <- (out$estimate - out$null)/out$se
    out$lower <- out$estimate + out$se * stats::qt(alpha/2, df = out$df)
    out$upper <- out$estimate + out$se * stats::qt(1-alpha/2, df = out$df)
    out$p.value = 2*(1-stats::pt(abs(out$statistic), df = out$df))
    rownames(out) <- pool.name
    attr(out,"contrast") <- pool.contrast
    return(out)
}

## * proportion.Wald_lmm
proportion.rbindWald_lmm <- function(object, name.method, method, qt = NULL, null, ci, df, alpha){

    


    ## ** estimate proportion
    object.ci <- confint(object, method = method, columns = c("estimate","se","df","lower","upper","statistic","null","p.value"))
    if(is.null(critical.threshold)){
        if(!is.null(attr(object.ci, "quantile"))){
            critical.threshold <- attr(object.ci, "quantile")
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

## * poolName.rbindWald_lmm
poolName.rbindWald_lmm <- function(object, method){

    ## ** extract
    table <- object$univariate

    ## ** build name
    Wald.name <- rownames(table)
    sep <- object$args$sep
    pool.splitname <- strsplit(Wald.name, split = sep, fixed = TRUE)

    if(all(lengths(pool.splitname)==2)){
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
    }else{
        ## <first test, last test>
        pool.name <- paste0("<",paste(Wald.name[1],Wald.name[length(Wald.name)],sep=", "),">")
    }

    ## **  export
    return(paste0(method,pool.name))
}


##----------------------------------------------------------------------
### pool.R ends here
