### reparametrize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 25 2021 (11:22) 
## Version: 
## Last-Updated: sep 20 2021 (14:53) 
##           By: Brice Ozenne
##     Update #: 347
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * reparametrize
reparametrize <- function(p, type, strata, time.levels,
                          Jacobian = TRUE, dJacobian = TRUE,
                          FUN = NULL,
                          transform.sigma = NULL,
                          transform.k = NULL,
                          transform.rho = NULL,
                          transform.names = TRUE){

    n.p <- length(p)
    name.p <- names(p)
    options <- LMMstar.options()
        
    ## ** normalize user arguments
    if(!is.null(FUN)){
        if(is.function(FUN)==FALSE){
            stop("Argument \'FUN\' must either be null or a function. \n")
        }
        if(!is.null(transform.sigma)){
            warning("Argument \'transform.sigma\' is ignored when argument \'FUN\' is specified. \n")
        }
        if(!is.null(transform.k)){
            warning("Argument \'transform.k\' is ignored when argument \'FUN\' is specified. \n")
        }
        if(!is.null(transform.rho)){
            warning("Argument \'transform.rho\' is ignored when argument \'FUN\' is specified. \n")
        }
        valid.args <- c("inverse","p","strata","time.levels","type")
        if(any(names(formals(FUN)) %in% valid.args == FALSE)){
            stop("Invalid argument for \'FUN\'. \n",
                 "Invalid argument(s): \'",paste(names(formals(FUN))[names(formals(FUN)) %in% valid.args == FALSE], collapse = "\' \'"),"\' \n",
                 "Valid arguments: \'",paste(valid.args, collapse = "\' \'"),"\' \n")
        }
        if(any(valid.args %in% names(formals(FUN)) == FALSE)){
            stop("Missing argument for \'FUN\': ",paste(valid.args[valid.args %in% names(formals(FUN)) == FALSE], collapse = "\' \'"),"\' \n")
        }
        if(as.numeric(dJacobian) %in% 0:2 == FALSE){
            stop("Argument \'dJacobian\' must be 0 (FALSE), 1 (TRUE), or 2. \n")
        }
        
        ## *** transformed parameter
        p.trans <- FUN(p = p, type = type, strata = strata, time.levels =  time.levels, inverse = FALSE)
        newname <- names(p.trans)
        names(p.trans) <- names(p)
        if(length(p.trans)!=n.p){
            stop("Argument \'FUN\' must output a vector of length ",n.p," (same length as argument \'p\'). \n")
        }

        if(transform.names){
            attr(p.trans,"newname") <- newname
        }

        ## *** jacobian
        if(Jacobian){
            if(any(abs(p-FUN(p = p.trans, type = type, strata = strata, time.levels =  time.levels, inverse = TRUE))>1e-6)){
                stop("The argument \'inverse\' of the function defined in argument \'FUN\' seems to lead to incorrect value. \n",
                     "FUN(FUN(p, inverse = FALSE),inverse = TRUE) should be equal to argument \'p\'. \n")
            }

            dFUN <- function(x, type, strata, time.levels, vectorize){
                out <- numDeriv::jacobian(function(pp){FUN(p = pp, type = type, strata = strata, time.levels =  time.levels, inverse = TRUE)},x)
                if(vectorize){
                    return(as.vector(out))
                }else{
                    return(out)
                }
            }
            attr(p.trans,"Jacobian") <- dFUN(p.trans, type = type, strata = strata, time.levels, vectorize = FALSE)
            dimnames(attr(p.trans,"Jacobian")) <- list(name.p,name.p)
        }
        
        ## *** derivative of the jacobian
        if(dJacobian==1){
            ddFUN <- function(x, type, strata, time.levels, vectorize){
                out <- numDeriv::jacobian(function(pp){dFUN(x = pp, type = type, strata = strata, time.levels =  time.levels, vectorize = TRUE)},x)
                if(vectorize){ ## convert to array
                    return(out)
                }else{
                    A.out <- array(NA, dim = rep(n.p, 3), dimnames = list(name.p,name.p,name.p))
                    for(iCol in 1:NCOL(out)){
                        A.out[,,iCol] <- out[,iCol]
                    }
                    return(A.out)
                }
            }        
            attr(p.trans,"dJacobian") <- ddFUN(p.trans, type = type, strata = strata, time.levels, vectorize = FALSE)
        }else if(dJacobian==2){
            stop("Argument \'dJacobian==2\' not implemented for numerical derivatives. \n")
        }

    }else{
        init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                                x.transform.sigma = options$transform.sigma, x.transform.k = options$transform.k, x.transform.rho = options$transform.rho)
            
        ls.out <- .reparametrize(p = p, type = type, strata = strata, time.levels = time.levels,
                                 Jacobian = Jacobian, dJacobian = dJacobian, inverse = FALSE,
                                 transform.sigma = init$transform.sigma,
                                 transform.k = init$transform.k,
                                 transform.rho = init$transform.rho,
                                 transform.names = transform.names)

        p.trans <- ls.out$p
        if(transform.names){
            attr(p.trans,"newname") <- ls.out$newname
        }
        if(Jacobian){
            attr(p.trans,"Jacobian") <- ls.out$Jacobian
        }
        if(dJacobian>0){
            attr(p.trans,"dJacobian") <- ls.out$dJacobian
        }
    }
    
    ## ** export
    return(p.trans)
}

## * .reparametrize
## TODO: - extend transform.rho="cov" to UN pattern (i.e. find k corresponding to rho)
.reparametrize <- function(p, type, strata, time.levels,
                           time.k, time.rho,
                           Jacobian, dJacobian, inverse,
                           transform.sigma,
                           transform.k,
                           transform.rho,
                           transform.names){
    
    if(transform.k %in% c("sd","logsd","var","logvar") && any("k" %in% type)){
        transform.sigma  <- switch(transform.k,
                                   "var" = "square",
                                   "logvar" = "logsquare",
                                   "sd" = "none",
                                   "logsd" = "log")
        if(!is.null(attr(transform.k,"arg")) && attr(transform.k,"arg")!=transform.sigma){
            warning("Argument \'transform.sigma\' ignored when argument \'transform.k\' is \"",transform.k,"\". \n")
        }
    }
    if(transform.rho %in% "cov" && any("rho" %in% type)){
        transform.sigma <- "square"
        transform.k <- "var"
        if(!is.null(attr(transform.sigma,"arg")) && attr(transform.sigma,"arg")!=transform.sigma){
            warning("Argument \'transform.sigma\' ignored when argument \'transform.rho\' is \"cov\". \n")
        }
        if(!is.null(attr(transform.k,"arg")) && attr(transform.k,"arg")!=transform.k){
            warning("Argument \'transform.k\' ignored when argument \'transform.rho\' is \"cov\". \n")
        }
    }
    if(transform.sigma == "remove" && (Jacobian || dJacobian)){
        stop("When argument \'transform.sigma\' is set to \"none\", arguments \'Jacobian\' and \'dJacobian\' should be set to FALSE. \n")
    }
    if(transform.sigma %in% "remove" && inverse){
        stop("When argument \'transform.sigma\' is set to \"none\", argument \'inverse\' should be set to FALSE. \n")
    }

    ## ** initialize
    n.p <- length(p)
    name.p <- names(p)
    index.sigma <- which(type=="sigma")
    index.k <- which(type=="k")
    index.rho <- which(type=="rho")
    out <- list(p = p,
                newname = name.p,
                transform = (transform.sigma!="none") || (transform.k!="none") || (transform.rho!="none"),
                transform.sigma = transform.sigma,
                transform.k = transform.k,
                transform.rho = transform.rho)
    if(Jacobian){
        out$Jacobian <- matrix(0, nrow = n.p, ncol = n.p, dimnames = list(name.p,name.p))
        diag(out$Jacobian) <- 1
    }        
    if(dJacobian>0){
        out$dJacobian <- array(0, dim = rep(n.p, 3), dimnames = list(name.p,name.p,name.p))
    }

    ## *** sigma
    if(length(index.sigma)>0){
        if(transform.sigma == "one"){
            out$p[index.sigma] <- 1
            if(transform.names){
                out$newname[index.sigma] <- "reference"
            }
        }else if(transform.sigma == "log"){
            if(inverse){
                out$p[index.sigma] <- exp(p[index.sigma])
            }else{
                out$p[index.sigma] <- log(p[index.sigma])
                if(transform.names){
                    out$newname[index.sigma] <- gsub("^sigma","log(sigma)",name.p[index.sigma])
                }
                if(Jacobian){
                    for(iSigma in index.sigma){
                        out$Jacobian[iSigma,iSigma] <- p[iSigma] ## exp(x) where x=log(sigma)
                    }
                }
                if(dJacobian==1){
                    for(iSigma in index.sigma){
                        out$dJacobian[iSigma,iSigma,iSigma] <- p[iSigma]
                    }
                }else if(dJacobian==2){
                    for(iSigma in index.sigma){
                        out$dJacobian[iSigma,iSigma,iSigma] <- 1
                    }
                }
            }
        }else if(transform.sigma == "square"){
            if(inverse){
                out$p[index.sigma] <- sqrt(p[index.sigma])
            }else{
                out$p[index.sigma] <- p[index.sigma]^2
                if(transform.names){
                    out$newname[index.sigma] <- gsub("^sigma","sigma^2",name.p[index.sigma])
                }
                if(Jacobian){
                    for(iSigma in index.sigma){
                        out$Jacobian[iSigma,iSigma] <- 1/(2*p[iSigma]) ## 1/(2*sqrt(x)) where x=sigma^2
                    }
                }
                if(dJacobian == 1){
                    for(iSigma in index.sigma){
                        out$dJacobian[iSigma,iSigma,iSigma] <- -1/(4*p[iSigma]^3)
                    }
                }else if(dJacobian == 2){
                    for(iSigma in index.sigma){
                        out$dJacobian[iSigma,iSigma,iSigma] <- -1/(2*p[iSigma]^2)
                    }
                }
            }
        }else if(transform.sigma == "logsquare"){
            if(inverse){
                out$p[index.sigma] <- sqrt(exp(p[index.sigma]))
            }else{
                out$p[index.sigma] <- log(p[index.sigma]^2)
                if(transform.names){
                    out$newname[index.sigma] <- gsub("^sigma","log(sigma^2)",name.p[index.sigma])
                }
                if(Jacobian){
                    for(iSigma in index.sigma){
                        out$Jacobian[iSigma,iSigma] <- p[iSigma]/2 ## exp(x/2)/2 where x=log(sigma^2)=2log(sigma)
                    }
                }
                if(dJacobian == 1){
                    for(iSigma in index.sigma){
                        out$dJacobian[iSigma,iSigma,iSigma] <- p[iSigma]/4
                    }
                }else if(dJacobian == 2){
                    for(iSigma in index.sigma){
                        out$dJacobian[iSigma,iSigma,iSigma] <- 1/2
                    }
                }
            }
        }
    }

    ## *** k
    if(length(index.k)>0){
        if(transform.k == "log"){
            if(inverse){
                out$p[index.k] <- exp(p[index.k])
            }else{
                out$p[index.k] <- log(p[index.k])
                if(transform.names){
                    out$newname[index.k] <- gsub("^k","log(k)",name.p[index.k])
                }
                if(Jacobian){
                    for(iK in index.k){
                        out$Jacobian[iK,iK] <- p[iK] ## exp(x) where x=log(k)
                    }
                }
                if(dJacobian == 1){
                    for(iK in index.k){
                        out$dJacobian[iK,iK,iK] <- p[iK]
                    }
                }else if(dJacobian == 2){
                    for(iK in index.k){
                        out$dJacobian[iK,iK,iK] <- 1
                    }
                }

            }
        }else if(transform.k == "square"){
            if(inverse){
                out$p[index.k] <- sqrt(p[index.k])
            }else{
                out$p[index.k] <- p[index.k]^2
                if(transform.names){
                    out$newname[index.k] <- gsub("^k","k^2",name.p[index.k])
                }
                if(Jacobian){
                    for(iK in index.k){
                        out$Jacobian[iK,iK] <- 1/(2*p[iK]) ## 1/(2*sqrt(x)) where x=k^2
                    }
                }
                if(dJacobian == 1){
                    for(iK in index.k){
                        out$dJacobian[iK,iK,iK] <- -1/(4*p[iK]^3)
                    }
                }else if(dJacobian == 2){
                    for(iK in index.k){
                        out$dJacobian[iK,iK,iK] <- -1/(2*p[iK]^2)
                    }
                }
            }
        }else if(transform.k == "logsquare"){
            if(inverse){
                out$p[index.k] <- sqrt(exp(p[index.k]))
            }else{
                out$p[index.k] <- log(p[index.k]^2)
                if(transform.names){
                    out$newname[index.k] <- gsub("^k","log(k^2)",name.p[index.k])
                }
                if(Jacobian){
                    for(iK in index.k){
                        out$Jacobian[iK,iK] <- p[iK]/2 ## exp(x/2)/2 where x=log(k^2)=2log(k)
                    }
                }
                if(dJacobian == 1){
                    for(iK in index.k){
                        out$dJacobian[iK,iK,iK] <- p[iK]/4
                    }
                }else if(dJacobian == 2){
                    for(iK in index.k){
                        out$dJacobian[iK,iK,iK] <- 1/2
                    }
                }
            }
        }else if(transform.k %in% c("sd","logsd","var","logvar")){
            index.param <- which(type %in% c("sigma","k"))
            name.param <- name.p[index.param]
            type.param <- type[index.param]
            strata.param <- strata[index.param]
            Ustrata <- unique(strata.param)

            for(iStrata in Ustrata){ ## iStrata <- 1
                iIndex  <- which(strata.param==iStrata)
                iType <- type.param[iIndex]
                iName <- name.param[iIndex]
                iName.sigma <- iName[iType=="sigma"]
                iName.k <- iName[iType=="k"]

                if(transform.names){
                    out$newname[match(iName.sigma,name.p)] <- switch(transform.sigma,
                                                                     "none" = paste0(iName.sigma,":",time.levels[1]),
                                                                     "log" = paste0("log(",iName.sigma,"):",time.levels[1]),
                                                                     "square" = paste0(iName.sigma,"^2:",time.levels[1]),
                                                                     "logsquare" = paste0("log(",iName.sigma,"^2):",time.levels[1])
                                                                     )
                }
                
                if(transform.k == "sd"){
                    if(inverse){
                        out$p[iName.k] <- p[iName.k]/p[iName.sigma]
                    }else{
                        out$p[iName.k] <- p[iName.sigma]*p[iName.k]
                        if(transform.names){
                            out$newname[match(c(iName.sigma,iName.k),name.p)] <- paste0(iName.sigma,":",time.levels)
                        }
                        if(Jacobian){
                            ## \sigma1, \sigma2/\sigma1, ...   so  1 -\sigma2/\sigma1^2       
                            ##                                     0  1/sigma1
                            ## WARNING we plug-in \sigma1 and \sigma2/\sigma1
                            out$Jacobian[iName.k,iName.sigma] <- -p[iName.k]/p[iName.sigma]
                            for(iK in iName.k){
                                out$Jacobian[iK,iK] <- 1/p[iName.sigma]
                            }
                        }
                        if(dJacobian == 1){
                            ## \sigma1, \sigma2/\sigma1, ...   so  0 -\sigma2/(2\sigma1^3)  0  -1/sigma^2    
                            ##                                     0  -1/sigma1^2           0  0
                            ## WARNING we plug-in \sigma1 and \sigma2/\sigma1
                            out$dJacobian[iName.k,iName.sigma,iName.sigma] <- 2*p[iName.k]/(p[iName.sigma]^2)
                            for(iK in iName.k){
                                out$dJacobian[iK,iName.sigma,iK] <- -1/p[iName.sigma]^2
                                out$dJacobian[iK,iK,iName.sigma] <- -1/p[iName.sigma]^2
                            }

                        }else if(dJacobian == 2){
                            out$dJacobian[iName.k,iName.sigma,iName.sigma] <- p[iName.k]/p[iName.sigma]^2
                            for(iK in iName.k){
                                out$dJacobian[iK,iName.sigma,iK] <- -1/p[iName.sigma]
                                out$dJacobian[iK,iK,iName.sigma] <- -1/p[iName.sigma]^2
                            }
                        }
                    }
                }else if(transform.k == "logsd"){
                    if(inverse){
                        out$p[iName.k] <- exp(p[iName.k]-p[iName.sigma])
                    }else{
                        out$p[iName.k] <- log(p[iName.k]*p[iName.sigma])
                        if(transform.names){
                            out$newname[match(c(iName.sigma,iName.k),name.p)] <- paste0("log(",iName.sigma,"):",time.levels)
                        }
                        if(Jacobian){
                            ## \exp(\sigma1), \exp(\sigma2 - \sigma1), ...   so  \exp(\sigma1) -\exp(\sigma2 - \sigma1)
                            ##                                                      0           \exp(\sigma2 - \sigma1)
                            ## WARNING we plug-in \exp(\sigma1) and \exp(\sigma2-\sigma1)
                            out$Jacobian[iName.k,iName.sigma] <- -p[iName.k]
                            for(iK in iName.k){
                                out$Jacobian[iK,iK] <- p[iK]
                            }
                        }
                        if(dJacobian == 1){
                            ## exp(\sigma1), exp(\sigma2-\sigma1), ...   so  \exp(\sigma) exp(\sigma2-\sigma1)    0 -exp(\sigma2-\sigma1) 
                            ##                                                         0  -exp(\sigma2-\sigma1)   0 exp(\sigma2-\sigma1) 
                            ## WARNING we plug-in exp(\sigma1) and exp(\sigma2-\sigma1)
                            out$dJacobian[iName.k,iName.sigma,iName.sigma] <- p[iName.k]
                            for(iK in iName.k){
                                out$dJacobian[iK,iName.sigma,iK] <- -p[iK]
                                out$dJacobian[iK,iK,iName.sigma] <- -p[iK]
                                out$dJacobian[iK,iK,iK] <- p[iK]
                            }
                        }else if(dJacobian == 2){
                            for(iK in iName.k){
                                out$dJacobian[iK,iName.sigma,iK] <- -1
                                out$dJacobian[iK,iK,iK] <- 1
                            }                            
                        }

                    }
                }else if(transform.k == "var"){
                    if(inverse){
                        out$p[iName.k] <- sqrt(p[iName.k]/p[iName.sigma])
                    }else{
                        out$p[iName.k] <- p[iName.sigma]^2*p[iName.k]^2
                        if(transform.names){
                            out$newname[match(c(iName.sigma,iName.k),name.p)] <- paste0(iName.sigma,"^2:",time.levels)
                        }
                        if(Jacobian){
                            ## sqrt(\sigma1), sqrt(\sigma2/\sigma1), ...   so  1/(2 sqrt(\sigma1)) -\sqrt(\sigma2)/(2 \sigma1^(3/2))       
                            ##                                                      0              1 / (2*\sqrt(\sigma2)\sqrt(\sigma1))
                            ## WARNING we plug-in sqrt(\sigma1) and sqrt(\sigma2/\sigma1)
                            out$Jacobian[iName.k,iName.sigma] <- -p[iName.k]/(2*p[iName.sigma]^2)
                            for(iK in iName.k){
                                out$Jacobian[iK,iK] <- 1/(2*p[iName.sigma]^2*p[iK])
                            }
                        }
                        if(dJacobian == 1){
                            ## sqrt(\sigma1), sqrt(\sigma2/\sigma1), ...   so  -1/(4 \sigma1^(3/2)) 3\sqrt(\sigma2)/(4 \sigma1^(5/2))      0 -1/(4 \sqrt(\sigma2) \sigma1^(3/2))        
                            ##                                                      0              - 1 / (4*\sqrt(\sigma2)\sigma1^(3/2))   0 -1 / (4*\sigma2^(3/2)\sqrt(\sigma1))
                            ## WARNING we plug-in sqrt(\sigma1) and sqrt(\sigma2/\sigma1)
                            out$dJacobian[iName.k,iName.sigma,iName.sigma] <- 3*p[iName.k]/(4*p[iName.sigma]^4)
                            for(iK in iName.k){
                                out$dJacobian[iK,iName.sigma,iK] <- -1/(4*p[iName.sigma]^4*p[iK])
                                out$dJacobian[iK,iK,iName.sigma] <- -1/(4*p[iName.sigma]^4*p[iK])
                                out$dJacobian[iK,iK,iK] <- -1/(4*p[iName.sigma]^4*p[iK]^3)
                            }
                        }else if(dJacobian == 2){
                            out$dJacobian[iName.k,iName.sigma,iName.sigma] <- p[iName.k]/p[iName.sigma]^3
                            for(iK in iName.k){
                                out$dJacobian[iK,iName.sigma,iK] <- -1/(2*p[iName.sigma]^2)
                                out$dJacobian[iK,iK,iName.sigma] <- -1/(p[iName.sigma]^3*p[iK])
                                out$dJacobian[iK,iK,iK] <- -1/(2*p[iName.sigma]^2*p[iK]^2)
                            }
                        }

                    }
                }else if(transform.k == "logvar"){
                    if(inverse){
                        out$p[iName.k] <- exp(p[iName.k]/2-p[iName.sigma]/2)
                    }else{
                        out$p[iName.k] <- 2*log(p[iName.k]*p[iName.sigma])
                        if(transform.names){
                            out$newname[match(c(iName.sigma,iName.k),name.p)] <- paste0("log(",iName.sigma,"):",time.levels)
                        }
                        if(Jacobian){
                            ## \exp(\sigma1/2), \exp(\sigma2/2 - \sigma1/2), ...   so  \exp(\sigma1/2)/2 -\exp(\sigma2/2 - \sigma1/2)/2
                            ##                                                          0                 \exp(\sigma2/2 - \sigma1/2)/2
                            ## WARNING we plug-in \exp(\sigma1/2) and \exp(\sigma2/2-\sigma1/2)
                            out$Jacobian[iName.k,iName.sigma] <- -p[iName.k]/2
                            for(iK in iName.k){
                                out$Jacobian[iK,iK] <- p[iK]/2
                            }
                        }
                        if(dJacobian == 1){
                            ## \exp(\sigma1/2), \exp(\sigma2/2 - \sigma1/2), ...   so  \exp(\sigma1/2)/4  \exp(\sigma2/2 - \sigma1/2)/4    0  -\exp(\sigma2/2 - \sigma1/2)/4
                            ##                                                          0                 -\exp(\sigma2/2 - \sigma1/2)/4   0  \exp(\sigma2/2 - \sigma1/2)/4
                            ## WARNING we plug-in \exp(\sigma1/2) and \exp(\sigma2/2-\sigma1/2)
                            out$dJacobian[iName.k,iName.sigma,iName.sigma] <- p[iName.k]/4
                            for(iK in iName.k){
                                out$dJacobian[iK,iName.sigma,iK] <- -p[iK]/4
                                out$dJacobian[iK,iK,iName.sigma] <- -p[iK]/4
                                out$dJacobian[iK,iK,iK] <- p[iK]/4
                            }
                        }else if(dJacobian == 2){
                            for(iK in iName.k){
                                out$dJacobian[iK,iName.sigma,iK] <- -1/2
                                out$dJacobian[iK,iK,iK] <- 1/2
                            }
                        }

                    }
                }
            }
        }
    }
            
    ## *** rho
    if(length(index.rho)>0){
        if(transform.rho == "atanh"){
            if(inverse){
                out$p[index.rho] <- tanh(p[index.rho])
            }else{
                out$p[index.rho] <- atanh(p[index.rho])
                if(transform.names){
                    out$newname[index.rho] <- paste0("atanh(",name.p[index.rho],")")
                }
                if(Jacobian){
                    for(iRho in index.rho){
                        out$Jacobian[iRho,iRho] <- (1-p[iRho]^2) ## (1-tanh(x)^2) where x=atanh(rho)
                    }
                }
                if(dJacobian == 1){
                    for(iRho in index.rho){
                        out$dJacobian[iRho,iRho,iRho] <- -2*p[iRho]*(1-p[iRho]^2)
                    }
                }else if(dJacobian == 2){
                    for(iRho in index.rho){
                        out$dJacobian[iRho,iRho,iRho] <- -2*p[iRho]
                    }
                }
            }
        }else if(transform.rho == "cov"){
            index.param <- which(type %in% c("rho","sigma","k"))
            name.param <- name.p[index.param]
            type.param <- type[index.param]
            strata.param <- strata[index.param]
            Ustrata <- unique(strata.param)

            for(iStrata in Ustrata){ ## iStrata <- 1
                iIndex  <- which(strata.param==iStrata)
                iType <- type.param[iIndex]
                iName <- name.param[iIndex]
                iName.sigma <- iName[iType=="sigma"]
                iName.k <- iName[iType=="k"]
                iName.rho <- iName[iType=="rho"]

                if(length(index.k)==0){
                    if(inverse){
                        out$p[iName.rho] <- p[iName.rho]/p[iName.sigma]
                    }else{
                        out$p[iName.rho] <- p[iName.rho]*p[iName.sigma]^2
                        if(transform.names){
                            out$newname[match(iName.sigma,name.p)] <- paste0(iName.sigma,"^2")
                            out$newname[match(iName.rho,name.p)] <- gsub(pattern = "^Rho|^cor",replacement = "cov",iName.rho,")")
                        }
                        if(Jacobian){
                            ##  \rho/\sigma leads to -\rho/\sigma^2 and 1/\sigma
                            ## WARNING we plug-in \sigma^2 and \rho*\sigma^2,
                            out$Jacobian[iName.rho,iName.sigma] <- -p[iName.rho]/p[iName.sigma]^2
                            for(iRho in iName.rho){
                                out$Jacobian[iRho,iRho] <- 1/p[iName.sigma]^2
                            }
                        }
                        if(dJacobian == 1){
                            ## sqrt(\sigma) \rho/\sigma, ...   so  -1/(4\sigma^(3/2))  2\rho/\sigma^3   0    -1/\sigma^2
                            ##                                     0                   -1/\sigma^2      0     0
                            ## WARNING we plug-in \sigma^2 and \rho*\sigma^2,
                            out$dJacobian[iName.rho,iName.sigma,iName.sigma] <- 2*p[iName.rho]/p[iName.sigma]^(4)
                            for(iRho in iName.rho){
                                out$dJacobian[iRho,iName.sigma,iRho] <- -1/p[iName.sigma]^4
                                out$dJacobian[iRho,iRho,iName.sigma] <- -1/p[iName.sigma]^4
                            }
                        }else if(dJacobian == 2){
                            out$dJacobian[iName.rho,iName.sigma,iName.sigma] <- 2*p[iName.rho]/p[iName.sigma]^3
                            for(iRho in iName.rho){
                                out$dJacobian[iRho,iName.sigma,iRho] <- -1/p[iName.sigma]^2
                                out$dJacobian[iRho,iRho,iName.sigma] <- -2/p[iName.sigma]^3
                            }
                        }
                    }
                }else if(is.null(time.rho)){
                    stop("Missing \'time.rho\' values. \n")
                }else{
                    if(inverse){
                        out$p[iName.rho] <- sapply(iName.rho, function(iRho){ ## iRho <- iName.rho[2]
                            iTime.k <- time.k[strata[names(time.k)] == strata[iRho]]
                            iName.sigma1 <- c(names(which(iTime.k==time.rho[1,iRho])),iName.sigma)[1]
                            iName.sigma2 <- c(names(which(iTime.k==time.rho[2,iRho])),iName.sigma)[1]
                            return(p[iRho]/sqrt(p[iName.sigma1]*p[iName.sigma2]))
                        })
                    }else{
                        out$p[iName.rho] <- sapply(iName.rho, function(iRho){ ## iRho <- iName.rho[1]
                            iTime.k <- time.k[strata[names(time.k)] == strata[iRho]]
                            iName.sigma1 <- names(which(iTime.k==time.rho[1,iRho]))
                            iName.sigma2 <- names(which(iTime.k==time.rho[2,iRho]))
                            if(length(iName.sigma1)>0){
                                iSigma1 <- p[iName.sigma1]
                            }else{
                                iSigma1 <- 1
                            }
                            if(length(iName.sigma2)>0){
                                iSigma2 <- p[iName.sigma2]
                            }else{
                                iSigma2 <- 1
                            }
                            return(p[iRho]*iSigma1*iSigma2*p[iName.sigma]^2)
                        })
                        if(transform.names){
                            out$newname[match(iName.rho,name.p)] <- gsub(pattern = "^Rho|^cor",replacement = "cov",iName.rho,")")
                        }
                        if(Jacobian){
                            stop("Not implemented in presence of multiple variance parameters. \n")
                            ## \rho/(\sigma1 \sigma2)  leads to  -\rho/(\sigma1^2 \sigma2) and -\rho/(\sigma1 \sigma2^2) and 1/(\sigma1 \sigma2)
                            ## WARNING we plug-in k1^2 sigma^2, k2^2 sigma^2, k1^2 k2^2 sigma^2 rho
                            ## factor.tempo <- matrix(1, nrow = 2, ncol = length(iName.rho),
                            ##                        dimnames = list(NULL, iName.rho))
                            ## index.tempo <- match(time.rho[,iName.rho,drop=FALSE], time.k[,iName.k])
                            ## factor.tempo[which(!is.na(index.tempo))] <- p[names(time.k)[stats::na.omit(index.tempo)]]

                            ## out$Jacobian[iName.rho,iName.sigma] <- 
                            ## for(iRho in iName.rho){ ## iRho <- iName.rho[1]
                            ##     iK <- names(time.k[iName.k])[time.k[iName.k] %in% time.rho[,iRho]]
                            ##     out$Jacobian[iRho,iName.sigma] <- - p[iRho]/(2*p[iName.sigma]^2)
                            ##     out$Jacobian[iRho,iK] <- - p[iRho]/(2*p[iName.sigma]^2*p[iK]^2)
                            ##     out$Jacobian[iRho,iRho] <- 1/(p[iName.sigma]^2*prod(factor.tempo[,iRho]))
                            ## }
                            
                            
                                
                            ## ##  \rho/\sigma leads to -\rho/\sigma^2 and 1/\sigma
                            ## ## WARNING we plug-in \sigma^2 and \rho*\sigma^2,
                            ## out$Jacobian[iName.rho,iName.sigma] <- -p[iName.rho]/p[iName.sigma]^2
                            ## for(iRho in iName.rho){
                            ##     out$Jacobian[iRho,iRho] <- 1/p[iName.sigma]^2
                            ## }
                        }
                        if(dJacobian == 1){
                            ## only off-diagonal elements (sigma,k,rho)
                            stop("Not implemented in presence of multiple variance parameters. \n")
                        }else if(dJacobian == 2){
                            ## diagonal elements
                            ## only off-diagonal elements (sigma,k,rho)
                            stop("Not implemented in presence of multiple variance parameters. \n")
                        }

                    }
                }
            }
        }
    }

    ## ** export
    return(out)
}

## * .init_transform
.init_transform <- function(transform.sigma, transform.k, transform.rho, 
                            x.transform.sigma, x.transform.k, x.transform.rho){

    ## ** transform
    if(is.null(transform.sigma)){
        transform.sigma.save <- transform.sigma
        transform.sigma <- x.transform.sigma
        attr(transform.sigma,"arg") <- transform.sigma.save
    }else{
        if(identical(transform.sigma,"")){
            transform.sigma <- "none"
        }
        transform.sigma <- match.arg(transform.sigma, c("none","one","log","square","logsquare","remove"))
        attr(transform.sigma,"arg") <- NULL
    }
    if(is.null(transform.k)){
        transform.k.save <- transform.k
        transform.k <- x.transform.k
        attr(transform.k,"arg") <- transform.k.save
    }else{
        if(identical(transform.k,"")){
            transform.k <- "none"
        }
        transform.k <- match.arg(transform.k, c("none","log","square","logsquare","sd","logsd","var","logvar","remove"))
        attr(transform.k,"arg") <- NULL
    }
    if(is.null(transform.rho)){
        transform.rho.save <- transform.rho
        transform.rho <- x.transform.rho
        attr(transform.rho,"arg") <- transform.rho.save
    }else{
        if(identical(transform.rho,"")){
            transform.rho <- "none"
        }
        transform.rho <- match.arg(transform.rho, c("none","atanh", "cov","remove"))
        attr(transform.rho,"arg") <- NULL
    }

    ## ** x.transform
    if(!is.null(x.transform.sigma) && !is.null(x.transform.k) && !is.null(x.transform.rho)){
        if(is.function(transform.sigma) || is.function(transform.k) || is.function(transform.rho)){
            test.notransform <- FALSE
        }else{
            test.notransform <- (transform.sigma==x.transform.sigma) && (transform.k==x.transform.k) && (transform.rho==x.transform.rho)
        }
    }else{
        test.notransform <- NULL
    }

    ## ** export
    return(list(transform.sigma = transform.sigma,
                transform.k = transform.k,
                transform.rho = transform.rho,
                test.notransform = test.notransform))
}

##----------------------------------------------------------------------
### reparametrize.R ends here
