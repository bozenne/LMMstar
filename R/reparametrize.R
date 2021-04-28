### reparametrize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 25 2021 (11:22) 
## Version: 
## Last-Updated: Apr 26 2021 (23:19) 
##           By: Brice Ozenne
##     Update #: 200
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
    options <- LMMstar.options()

    n.p <- length(p)
    name.p <- names(p)

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
        if(dJacobian){
            if(any(abs(p-FUN(p = p.trans, type = type, strata = strata, time.levels =  time.levels, inverse = TRUE))>1e-6)){
                stop("The argument \'inverse\' of the function defined in argument \'FUN\' seems to lead to incorrect value. \n",
                     "FUN(FUN(p, inverse = FALSE),inverse = TRUE) should be equal to argument \'p\'. \n")
            }

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
        }

    }else{
        if(is.null(transform.sigma)){
            transform.sigma <- options$transform.sigma
        }else{
            transform.sigma <- match.arg(transform.sigma, c("none","log","square","logsquare"))
        }
        if(is.null(transform.k)){
            transform.k <- options$transform.k
        }else{
            transform.k <- match.arg(transform.k, c("none","log","square","logsquare","sd","logsd","var","logvar"))
        }
        if(transform.k %in% c("sd","logsd","var","logvar") && any(type=="k") && transform.sigma != "none"){
            warning("Argument \'transform.sigma\' ignored when argument \'transform.k\' set to \"sd\", \"logsd\", \"var\", or \"logvar\". \n")
        }
        if(is.null(transform.rho)){
            transform.rho <- options$transform.rho
        }else{
            transform.rho <- match.arg(transform.rho, c("none","atanh", "cov"))
        }
        if(transform.rho %in% c("cov") && any(type=="rho") && transform.sigma != "none"){
            warning("Argument \'transform.sigma\' ignored when argument \'transform.rho\' set to \"cov\". \n")
        }

        p.trans <- .reparametrize(p = p, type = type, strata = strata, time.levels = time.levels,
                                  Jacobian = Jacobian, dJacobian = dJacobian, inverse = FALSE,
                                  transform.sigma = transform.sigma,
                                  transform.k = transform.k,
                                  transform.rho = transform.rho,
                                  transform.names = transform.names)
    }
    
    ## ** export
    return(p.trans)
}

## * .reparametrize
## TODO: - extend transform.rho="cov" to UN pattern (i.e. find k corresponding to rho)
.reparametrize <- function(p, type, strata, time.levels,
                           Jacobian, dJacobian, inverse,
                           transform.sigma,
                           transform.k,
                           transform.rho,
                           transform.names){
    
    ## ** initialize
    n.p <- length(p)
    name.p <- names(p)
    index.sigma <- which(type=="sigma")
    index.k <- which(type=="k")
    index.rho <- which(type=="rho")

    p.trans <- p
    if(transform.names){
        attr(p.trans,"newname") <- name.p
    }
    if(Jacobian){
        attr(p.trans, "Jacobian") <- matrix(0, nrow = n.p, ncol = n.p, dimnames = list(name.p,name.p))
        diag(attr(p.trans, "Jacobian")) <- 1
    }        
    if(dJacobian){
        attr(p.trans, "dJacobian") <- array(0, dim = rep(n.p, 3), dimnames = list(name.p,name.p,name.p))
    }

    ## *** sigma
    ntest.k <- (transform.k %in% c("sd","logsd","var","logvar") == FALSE) || (length(index.k) == 0)
    ntest.rho <- (transform.rho %in% c("cov") == FALSE) || (length(index.rho) == 0)
    if(ntest.k && ntest.rho){
        if(transform.sigma == "log"){
            if(inverse){
                p.trans[index.sigma] <- exp(p[index.sigma])
            }else{
                p.trans[index.sigma] <- log(p[index.sigma])
                if(transform.names){
                    attr(p.trans,"newname")[index.sigma] <- paste0("log(",name.p[index.sigma],")")
                }
                if(Jacobian){
                    attr(p.trans, "Jacobian")[index.sigma,index.sigma] <- p[index.sigma] ## exp(x) where x=log(sigma)
                }
                if(dJacobian){
                    attr(p.trans, "dJacobian")[index.sigma,index.sigma,index.sigma] <- p[index.sigma]
                }
            }
        }else if(transform.sigma == "square"){
            if(inverse){
                p.trans[index.sigma] <- sqrt(p[index.sigma])
            }else{
                p.trans[index.sigma] <- p[index.sigma]^2
                if(transform.names){
                    attr(p.trans,"newname")[index.sigma] <- paste0(name.p[index.sigma],"^2")
                }
                if(Jacobian){
                    attr(p.trans, "Jacobian")[index.sigma,index.sigma] <- 1/(2*p[index.sigma]) ## 1/(2*sqrt(x)) where x=sigma^2
                }
                if(dJacobian){
                    attr(p.trans, "dJacobian")[index.sigma,index.sigma,index.sigma] <- -1/(4*p[index.sigma]^3)
                }
            }
        }else if(transform.sigma == "logsquare"){
            if(inverse){
                p.trans[index.sigma] <- sqrt(exp(p[index.sigma]))
            }else{
                p.trans[index.sigma] <- log(p[index.sigma]^2)
                if(transform.names){
                    attr(p.trans,"newname")[index.sigma] <- paste0("log(",name.p[index.sigma],"^2)")
                }
                if(Jacobian){
                    attr(p.trans, "Jacobian")[index.sigma,index.sigma] <- p[index.sigma]/2 ## exp(x/2)/2 where x=log(sigma^2)=2log(sigma)
                }
                if(dJacobian){
                    attr(p.trans, "dJacobian")[index.sigma,index.sigma,index.sigma] <- p[index.sigma]/4
                }
            }
        }
    }

    ## *** k
    if(length(index.k)){
        if(transform.k == "log"){
            if(inverse){
                p.trans[index.k] <- exp(p[index.k])
            }else{
                p.trans[index.k] <- log(p[index.k])
                if(transform.names){
                    attr(p.trans,"newname")[index.k] <- paste0("log(",name.p[index.k],")")
                }
                if(Jacobian){
                    attr(p.trans, "Jacobian")[index.k,index.k] <- p[index.k] ## exp(x) where x=log(k)
                }
                if(dJacobian){
                    attr(p.trans, "dJacobian")[index.k,index.k,index.k] <- p[index.k]
                }

            }
        }else if(transform.k == "square"){
            if(inverse){
                p.trans[index.k] <- sqrt(p[index.k])
            }else{
                p.trans[index.k] <- p[index.k]^2
                if(transform.names){
                    attr(p.trans,"newname")[index.k] <- paste0(name.p[index.k],"^2")
                }
                if(Jacobian){
                    attr(p.trans, "Jacobian")[index.k,index.k] <- 1/(2*p[index.k]) ## 1/(2*sqrt(x)) where x=k^2
                }
                if(dJacobian){
                    attr(p.trans, "dJacobian")[index.k,index.k,index.k] <- -1/(4*p[index.k]^3)
                }
            }
        }else if(transform.k == "logsquare"){
            if(inverse){
                p.trans[index.k] <- sqrt(exp(p[index.k]))
            }else{
                p.trans[index.k] <- log(p[index.k]^2)
                if(transform.names){
                    attr(p.trans,"newname")[index.k] <- paste0("log(",name.p[index.k],"^2)")
                }
                if(Jacobian){
                    attr(p.trans, "Jacobian")[index.k,index.k] <- p[index.k]/2 ## exp(x/2)/2 where x=log(k^2)=2log(k)
                }
                if(dJacobian){
                    attr(p.trans, "dJacobian")[index.k,index.k,index.k] <- p[index.k]/4
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
            
                if(transform.k == "sd"){
                    if(inverse){
                        p.trans[iName.k] <- p[iName.k]/p[iName.sigma]
                    }else{
                        p.trans[iName.k] <- p[iName.sigma]*p[iName.k]
                        if(transform.names){
                            attr(p.trans, "newname")[match(c(iName.sigma,iName.k),name.p)] <- paste0(iName.sigma,":",time.levels)
                        }
                        if(Jacobian){
                            ## \sigma1, \sigma2/\sigma1, ...   so  1 -\sigma2/\sigma1^2       
                            ##                                     0  1/sigma1
                            ## WARNING we plug-in \sigma1 and \sigma2/\sigma1
                            attr(p.trans, "Jacobian")[iName.k,iName.sigma] <- -p[iName.k]/p[iName.sigma]
                            for(iK in iName.k){
                                attr(p.trans, "Jacobian")[iK,iK] <- 1/p[iName.sigma]
                            }
                        }
                        if(dJacobian){
                            ## \sigma1, \sigma2/\sigma1, ...   so  0 -\sigma2/(2\sigma1^3)  0  -1/sigma^2    
                            ##                                     0  -1/sigma1^2           0  0
                            ## WARNING we plug-in \sigma1 and \sigma2/\sigma1
                            attr(p.trans, "dJacobian")[iName.k,iName.sigma,iName.sigma] <- 2*p[iName.k]/(p[iName.sigma]^2)
                            for(iK in iName.k){
                                attr(p.trans, "dJacobian")[iK,iName.sigma,iK] <- -1/p[iName.sigma]^2
                                attr(p.trans, "dJacobian")[iK,iK,iName.sigma] <- -1/p[iName.sigma]^2
                            }

                        }
                    }
                }else if(transform.k == "logsd"){
                    if(inverse){
                        p.trans[iName.sigma] <- exp(p[iName.sigma])
                        p.trans[iName.k] <- exp(p[iName.k]-p[iName.sigma])
                    }else{
                        p.trans[iName.sigma] <- log(p[iName.sigma])
                        p.trans[iName.k] <- log(p[iName.k]*p[iName.sigma])
                        if(transform.names){
                            attr(p.trans, "newname")[match(c(iName.sigma,iName.k),name.p)] <- paste0("log(",iName.sigma,"):",time.levels)
                        }
                        if(Jacobian){
                            ## \exp(\sigma1), \exp(\sigma2 - \sigma1), ...   so  \exp(\sigma1) -\exp(\sigma2 - \sigma1)
                            ##                                                      0           \exp(\sigma2 - \sigma1)
                            ## WARNING we plug-in \exp(\sigma1) and \exp(\sigma2-\sigma1)
                            attr(p.trans, "Jacobian")[iName.sigma,iName.sigma] <- p[iName.sigma]
                            attr(p.trans, "Jacobian")[iName.k,iName.sigma] <- -p[iName.k]
                            for(iK in iName.k){
                                attr(p.trans, "Jacobian")[iK,iK] <- p[iK]
                            }
                        }
                        if(dJacobian){
                            ## exp(\sigma1), exp(\sigma2-\sigma1), ...   so  \exp(\sigma) exp(\sigma2-\sigma1)    0 -exp(\sigma2-\sigma1) 
                            ##                                                         0  -exp(\sigma2-\sigma1)   0 exp(\sigma2-\sigma1) 
                            ## WARNING we plug-in exp(\sigma1) and exp(\sigma2-\sigma1)
                            attr(p.trans, "dJacobian")[iName.sigma,iName.sigma,iName.sigma] <- p[iName.sigma]
                            attr(p.trans, "dJacobian")[iName.k,iName.sigma,iName.sigma] <- p[iName.k]
                            for(iK in iName.k){
                                attr(p.trans, "dJacobian")[iK,iName.sigma,iK] <- -p[iK]
                                attr(p.trans, "dJacobian")[iK,iK,iName.sigma] <- -p[iK]
                                attr(p.trans, "dJacobian")[iK,iK,iK] <- p[iK]
                            }
                        }

                    }
                }else if(transform.k == "var"){
                    if(inverse){
                        p.trans[iName.sigma] <- sqrt(p[iName.sigma])
                        p.trans[iName.k] <- sqrt(p[iName.k]/p[iName.sigma])
                    }else{
                        p.trans[iName.sigma] <- p[iName.sigma]^2
                        p.trans[iName.k] <- p[iName.sigma]^2*p[iName.k]^2
                        if(transform.names){
                            attr(p.trans, "newname")[match(c(iName.sigma,iName.k),name.p)] <- paste0(iName.sigma,"^2:",time.levels)
                        }
                        if(Jacobian){
                            ## sqrt(\sigma1), sqrt(\sigma2/\sigma1), ...   so  1/(2 sqrt(\sigma1)) -\sqrt(\sigma2)/(2 \sigma1^(3/2))       
                            ##                                                      0              1 / (2*\sqrt(\sigma2)\sqrt(\sigma1))
                            ## WARNING we plug-in sqrt(\sigma1) and sqrt(\sigma2/\sigma1)
                            attr(p.trans, "Jacobian")[iName.sigma,iName.sigma] <- 1/(2*p[iName.sigma])
                            attr(p.trans, "Jacobian")[iName.k,iName.sigma] <- -p[iName.k]/(2*p[iName.sigma]^2)
                            for(iK in iName.k){
                                attr(p.trans, "Jacobian")[iK,iK] <- 1/(2*p[iName.sigma]^2*p[iK])
                            }
                        }
                        if(dJacobian){
                            ## sqrt(\sigma1), sqrt(\sigma2/\sigma1), ...   so  -1/(4 \sigma1^(3/2)) 3\sqrt(\sigma2)/(4 \sigma1^(5/2))      0 -1/(4 \sqrt(\sigma2) \sigma1^(3/2))        
                            ##                                                      0              - 1 / (4*\sqrt(\sigma2)\sigma1^(3/2))   0 -1 / (4*\sigma2^(3/2)\sqrt(\sigma1))
                            ## WARNING we plug-in sqrt(\sigma1) and sqrt(\sigma2/\sigma1)
                            attr(p.trans, "dJacobian")[iName.sigma,iName.sigma,iName.sigma] <- -1/(4*p[iName.sigma]^3)
                            attr(p.trans, "dJacobian")[iName.k,iName.sigma,iName.sigma] <- 3*p[iName.k]/(4*p[iName.sigma]^4)
                            for(iK in iName.k){
                                attr(p.trans, "dJacobian")[iK,iName.sigma,iK] <- -1/(4*p[iName.sigma]^4*p[iK])
                                attr(p.trans, "dJacobian")[iK,iK,iName.sigma] <- -1/(4*p[iName.sigma]^4*p[iK])
                                attr(p.trans, "dJacobian")[iK,iK,iK] <- -1/(4*p[iName.sigma]^4*p[iK]^3)
                            }
                        }

                    }
                }else if(transform.k == "logvar"){
                    if(inverse){
                        p.trans[iName.sigma] <- exp(p[iName.sigma]/2)
                        p.trans[iName.k] <- exp(p[iName.k]/2-p[iName.sigma]/2)
                    }else{
                        p.trans[iName.sigma] <- 2*log(p[iName.sigma])
                        p.trans[iName.k] <- 2*log(p[iName.k]*p[iName.sigma])
                        if(transform.names){
                            attr(p.trans, "newname")[match(c(iName.sigma,iName.k),name.p)] <- paste0("log(",iName.sigma,"):",time.levels)
                        }
                        if(Jacobian){
                            ## \exp(\sigma1/2), \exp(\sigma2/2 - \sigma1/2), ...   so  \exp(\sigma1/2)/2 -\exp(\sigma2/2 - \sigma1/2)/2
                            ##                                                          0                 \exp(\sigma2/2 - \sigma1/2)/2
                            ## WARNING we plug-in \exp(\sigma1/2) and \exp(\sigma2/2-\sigma1/2)
                            attr(p.trans, "Jacobian")[iName.sigma,iName.sigma] <- p[iName.sigma]/2
                            attr(p.trans, "Jacobian")[iName.k,iName.sigma] <- -p[iName.k]/2
                            for(iK in iName.k){
                                attr(p.trans, "Jacobian")[iK,iK] <- p[iK]/2
                            }
                        }
                        if(dJacobian){
                            ## \exp(\sigma1/2), \exp(\sigma2/2 - \sigma1/2), ...   so  \exp(\sigma1/2)/4  \exp(\sigma2/2 - \sigma1/2)/4    0  -\exp(\sigma2/2 - \sigma1/2)/4
                            ##                                                          0                 -\exp(\sigma2/2 - \sigma1/2)/4   0  \exp(\sigma2/2 - \sigma1/2)/4
                            ## WARNING we plug-in \exp(\sigma1/2) and \exp(\sigma2/2-\sigma1/2)
                            attr(p.trans, "dJacobian")[iName.sigma,iName.sigma,iName.sigma] <- p[iName.sigma]/4
                            attr(p.trans, "dJacobian")[iName.k,iName.sigma,iName.sigma] <- p[iName.k]/4
                            for(iK in iName.k){
                                attr(p.trans, "dJacobian")[iK,iName.sigma,iK] <- -p[iK]/4
                                attr(p.trans, "dJacobian")[iK,iK,iName.sigma] <- -p[iK]/4
                                attr(p.trans, "dJacobian")[iK,iK,iK] <- p[iK]/4
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
                p.trans[index.rho] <- tanh(p[index.rho])
            }else{
                p.trans[index.rho] <- atanh(p[index.rho])
                if(transform.names){
                    attr(p.trans,"newname")[index.rho] <- paste0("atanh(",name.p[index.rho],")")
                }
                if(Jacobian){
                    attr(p.trans, "Jacobian")[index.rho,index.rho] <- (1-p[index.rho]^2) ## (1-tanh(x)^2) where x=atanh(rho)
                }
                if(dJacobian){
                    attr(p.trans, "dJacobian")[index.rho,index.rho,index.rho] <- -2*p[index.rho]*(1-p[index.rho]^2)
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
                        p.trans[iName.sigma] <- sqrt(p[iName.sigma])
                        p.trans[iName.rho] <- p[iName.rho]/p[iName.sigma]
                    }else{
                        p.trans[iName.sigma] <- p[iName.sigma]^2
                        p.trans[iName.rho] <- p[iName.rho]*p[iName.sigma]^2
                        if(transform.names){
                            attr(p.trans,"newname")[iName.sigma] <- paste0(name.p[iName.sigma],"^2")
                            attr(p.trans,"newname")[iName.rho] <- gsub(pattern = "^Rho|^cor",replacement = "cov",name.p[iName.rho],")")
                        }
                        if(Jacobian){
                            ## sqrt(\sigma) \rho/\sigma, ...   so  1/(2\sqrt(\sigma))  -\rho/\sigma^2
                            ##                                     0                    1/\sigma
                            ## WARNING we plug-in \sigma^2 and \rho*\sigma^2,
                            attr(p.trans, "Jacobian")[iName.sigma,iName.sigma] <- 1/(2*p[iName.sigma])
                            attr(p.trans, "Jacobian")[iName.rho,iName.sigma] <- -p[iName.rho]/p[iName.sigma]^2
                            for(iRho in iName.rho){
                                attr(p.trans, "Jacobian")[iRho,iRho] <- 1/p[iName.sigma]^2
                            }
                        }
                        if(dJacobian){
                            ## sqrt(\sigma) \rho/\sigma, ...   so  -1/(4\sigma^(3/2))  2\rho/\sigma^3   0    -1/\sigma^2
                            ##                                     0                   -1/\sigma^2      0     0
                            ## WARNING we plug-in \sigma^2 and \rho*\sigma^2,
                            attr(p.trans, "dJacobian")[iName.sigma,iName.sigma,iName.sigma] <- -1/(4*p[iName.sigma]^3)
                            attr(p.trans, "dJacobian")[iName.rho,iName.sigma,iName.sigma] <- 2*p[iName.rho]/p[iName.sigma]^(4)
                            for(iRho in iName.rho){
                                attr(p.trans, "dJacobian")[iRho,iName.sigma,iRho] <- -1/p[iName.sigma]^4
                                attr(p.trans, "dJacobian")[iRho,iRho,iName.sigma] <- -1/p[iName.sigma]^4
                            }
                        }
                    }
                }else{
                    stop("Not yet implemented in presence of multiple variance parameters. \n")
                    browser()
                }
            }
        }
    }

    return(p.trans)
}

##----------------------------------------------------------------------
### reparametrize.R ends here
