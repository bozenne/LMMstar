### poolWald.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 11 2024 (11:56) 
## Version: 
## Last-Updated: jul 15 2024 (16:35) 
##           By: Brice Ozenne
##     Update #: 93
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * poolWald.mlmm
poolWald.mlmm <- function(object, p, index, method, name.method, ci, df, alpha){

    ## ** extract information
    table <- model.tables(object, p = p, columns = c("estimate","df","null"))[index,,drop=FALSE]

    if(ci || method %in% c("pool.se","pool.gls","pool.gls1")){
        transform.sigma <-  object$args$transform.sigma
        transform.k <-  object$args$transform.k
        transform.rho <-  object$args$transform.rho
        type.information <-  object$object$type.information
        robust <-  object$args$robust

        Sigma <- vcov(object, p = p, type.information = type.information, robust = robust,
                      transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)[index,index,drop=FALSE]
        independence <- attr(object$object,"independence")

        if(ci && is.null(df)){
            if(method %in% "pool.rubin"){
                df <- TRUE
            }else if(independence & any(is.infinite(table$df)==FALSE) & method %in% c("average","pool.fixse")){
                df <- TRUE
            }else{
                df <- FALSE
            }
        }
    }

    ## ** name for the pooled estimator
    if(!is.null(name.method)){
        pool.name <- name.method
    }else{
        pool.name <- poolName.mlmm(object, index = index, method = method)
    }

    ## ** point estimate
    n.test <- NROW(table)
    if(method %in% c("average","pool.rubin")){

        pool.contrast <- matrix(1/n.test, nrow = 1, ncol = n.test)

    }else if(method %in% c("pool.fixse","pool.se")){

        iSigma.diag <- 1/diag(Sigma)
        pool.contrast <- matrix(iSigma.diag/sum(iSigma.diag), nrow = 1, ncol = n.test)

    }else if(method %in% c("pool.gls","pool.gls1")){

        poolGLS <- function(Sigma, method, tol = 1e-10){
            ## Note: without tuncation the weights can be computed as
            ## rowSums(solve(Sigma))/sum(solve(Sigma))

            ## **** truncate eigen value decomposition
            Sigma.eigen <- eigen(Sigma)
            index.subset <- which(abs(Sigma.eigen$values) > tol)

            if(length(index.subset)==0){
                stop("All eigenvalues  of the variance-covariance matrix are close to 0 (<",tol,"). \n",sep="")
            }else if(any(abs(Sigma.eigen$values) <= tol)){
                error <-  length(Sigma.eigen$values)-length(index.subset)
            }else{
                error <- NULL
            }

            lambda.k <- Sigma.eigen$values[index.subset]
            q.kj <- Sigma.eigen$vector[,index.subset,drop=FALSE]

            ## **** evaluate weigths
            qbar.k <- colSums(q.kj)
            w.k <- qbar.k^2/lambda.k
            out <- rbind(rowSums(sweep(q.kj, FUN = "*", MARGIN = 2, STATS = w.k/qbar.k))/sum(w.k))
            ## shortcut for
            ## out <- rbind(colSums(sweep(t(q.kj), FUN = "*", MARGIN = 1, STATS = w.k/qbar.k))/sum(w.k))

            if(method == "pool.gls1" && max(abs(out))>1){
                ## ensure no weight greater than 1
                index.max <- which.max(abs(out))
                kappaPwmax <- max((1-n.test*out)/(1-sign(out)*n.test))
                out <- rbind(rep((1-1/kappaPwmax)/n.test,n.test) + out[1,]/kappaPwmax)
            }

            ## **** export
            attr(out,"error") <- error
            return(out)
        }

        pool.contrast <- poolGLS(Sigma = Sigma, method = method)        
    } 

    ## ** variance
    if(!ci){

        pool.se  <- NA

    }else if(method %in% c("average","pool.fixse")){
        
        pool.se <- sqrt(as.double(pool.contrast %*% Sigma %*% t(pool.contrast)))

    }else if(method %in% c("pool.se","pool.gls","pool.gls1")){

        
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
        
    }else if(method == "pool.rubin"){
        pool.U <- mean(diag(Sigma))
        pool.B <- sum((table$estimate - mean(table$estimate))^2)/(n.test-1)
        pool.se <- sqrt(pool.U + (1 + 1/n.test) * pool.B)
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

    }else if(abs(diff(range(table$df)))<0.1){

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
            attr(vcov(iM, effects = "all", df = 2,
                 transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE),"dVcov")
        })
        lmm.df <- mapply(iC = lmm.ls.contrast, iVcov = lmm.ls.vcov, iDvcov = lmm.ls.dVcov, FUN = function(iC,iVcov,iDvcov){
            iDf <- .dfX(X.beta = iC, vcov.param = iVcov, dVcov.param = iDvcov, return.vcov = TRUE)
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
            iDvcov <- attr(vcov(object$model[[iBy]], effects = "all", df = 2,
                                transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE),"dVcov")
            theta.dVcov[attr(theta.Sigma,"by")==iBy,attr(theta.Sigma,"by")==iBy,attr(theta.Sigma,"by")==iBy] <- iDvcov
        }
        pool.df <- .dfX(X.beta = theta.grad, vcov.param = theta.Sigma, dVcov.param = theta.dVcov)
    }

    ## ** export
    out <- table[1,,drop=FALSE]
    out[setdiff(names(out),c("parameter","type","test","estimate","se","df","statistic","lower","upper","null","p.value"))] <- NA
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

## * poolName.mlmm
poolName.mlmm <- function(object, index, method){

    ## ** extract
    table <- object$univariate[index,,drop=FALSE]

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
### poolWald.R ends here
