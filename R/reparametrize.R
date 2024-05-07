### reparametrize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 25 2021 (11:22) 
## Version: 
## Last-Updated: maj  7 2024 (11:27) 
##           By: Brice Ozenne
##     Update #: 791
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * reparametrize
reparametrize <- function(p, type, level, sigma, k.x, k.y, 
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
        valid.args <- c("inverse","p","type","sigma","k.x","k.y")

        if(any(names(formals(FUN)) %in% c(valid.args,"...") == FALSE)){
            stop("Invalid argument for \'FUN\'. \n",
                 "Invalid argument(s): \'",paste(names(formals(FUN))[names(formals(FUN)) %in% valid.args == FALSE], collapse = "\' \'"),"\' \n",
                 "Valid arguments: \'",paste(valid.args, collapse = "\' \'"),"\' \n")
        }
        if(any(valid.args %in% names(formals(FUN)) == FALSE) && "..." %in% names(formals(FUN)) == FALSE){
            stop("Missing argument for \'FUN\': ",paste(valid.args[valid.args %in% names(formals(FUN)) == FALSE], collapse = "\' \'"),"\' \n")
        }
        if(as.numeric(dJacobian) %in% 0:2 == FALSE){
            stop("Argument \'dJacobian\' must be 0 (FALSE), 1 (TRUE), or 2. \n")
        }
        
        ## *** transformed parameter
        p.trans <- FUN(p = p, type = type, sigma = sigma, k.x = k.x, k.y = k.y, inverse = FALSE)
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
            if(any(abs(p-FUN(p = p.trans, type = type, sigma = sigma, k.x = k.x, k.y = k.y, inverse = TRUE))>1e-6)){
                stop("The argument \'inverse\' of the function defined in argument \'FUN\' seems to lead to incorrect value. \n",
                     "FUN(FUN(p, inverse = FALSE),inverse = TRUE) should be equal to argument \'p\'. \n")
            }

            dFUN <- function(x, type, sigma, k.x, k.y, vectorize){
                out <- numDeriv::jacobian(function(pp){FUN(p = pp, type = type, sigma = sigma, k.x = k.x, k.y = k.y, inverse = TRUE)},x)
                if(vectorize){
                    return(as.vector(out))
                }else{
                    return(out)
                }
            }
            attr(p.trans,"Jacobian") <- dFUN(p.trans, type = type, sigma = sigma, k.x = k.x, k.y = k.y, vectorize = FALSE)
            dimnames(attr(p.trans,"Jacobian")) <- list(name.p,name.p)
        }
        
        ## *** derivative of the jacobian
        if(dJacobian==1){
            ddFUN <- function(x, type, sigma, k.x, k.y, vectorize){
                out <- numDeriv::jacobian(function(pp){dFUN(x = pp, type = type, sigma = sigma, k.x = k.x, k.y = k.y, vectorize = TRUE)},x)
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
            attr(p.trans,"dJacobian") <- ddFUN(p.trans, type = type, sigma = sigma, k.x = k.x, k.y = k.y, vectorize = FALSE)
        }else if(dJacobian==2){
            stop("Argument \'dJacobian==2\' not implemented for numerical derivatives. \n")
        }

    }else{
        init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                                x.transform.sigma = options$transform.sigma, x.transform.k = options$transform.k, x.transform.rho = options$transform.rho)
        attr(init$transform.sigma,"arg") <- transform.sigma
        attr(init$transform.k,"arg") <- transform.k
        attr(init$transform.rho,"arg") <- transform.rho

        ls.out <- .reparametrize(p = p, type = type, level = level, sigma = sigma, k.x = k.x, k.y = k.y,
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
.reparametrize <- function(p, type, level,
                           sigma, k.x, k.y,
                           Jacobian, dJacobian, inverse,
                           transform.sigma,
                           transform.k,
                           transform.rho,
                           transform.names){

    arg.transform.sigma <- attr(transform.sigma,"arg")
    arg.transform.k <- attr(transform.k,"arg")

    if(transform.rho %in% "cov" && any("rho" %in% type)){
        transform.sigma <- "square"
        transform.k <- "var"
        if(!is.null(arg.transform.sigma) && arg.transform.sigma!=transform.sigma){
            warning("Argument \'transform.sigma\' ignored when argument \'transform.rho\' is \"cov\". \n")
        }
        if(!is.null(arg.transform.k) && arg.transform.k!=transform.k){
            warning("Argument \'transform.k\' ignored when argument \'transform.rho\' is \"cov\". \n")
        }
    }
    if(transform.k %in% c("sd","logsd","var","logvar") && any("k" %in% type)){
        transform.sigma <- switch(transform.k,
                                  "sd" = "none",
                                  "logsd" = "log",
                                  "var" = "square",
                                  "logvar" = "logsquare")

        if(!is.null(arg.transform.sigma) && arg.transform.sigma!=transform.sigma){
            warning(paste0("Argument \'transform.sigma\' ignored when argument \'transform.k\' is \"",transform.k,"\". \n"))
        }
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
    
 
    ## *** sigma (or k with log/square/logsquare)
    if(length(index.sigma)>0){

        transform.sigma <- switch(transform.sigma,
                                  "none" = .reparametrize.none,
                                  "one" = .reparametrize.one,
                                  "log" = .reparametrize.log,
                                  "square" = .reparametrize.square,
                                  "logsquare" = .reparametrize.logsquare)

        if(!transform.names){
            level.sigma <- NULL
        }else if(transform.k %in% c("sd","logsd","var","logvar")){
            level.sigma <- sapply(level,"[[",1)
        }else{
            level.sigma <- stats::setNames(rep(NA,n.p),name.p)
            pos.column <- which(grepl(":",name.p,fixed=TRUE))
            if(length(pos.column)>1){ ## add strata
                level.sigma[pos.column] <- sapply(strsplit(name.p[pos.column],":", fixed = TRUE), function(iString){
                    paste(c(":",iString[-1]),collapse="")
                })
            }
        }
        out <- transform.sigma(p = p, out = out,
                               index = index.sigma,
                               level = level.sigma, type = "sigma",
                               inverse = inverse, transform.names = transform.names, Jacobian = Jacobian, dJacobian = dJacobian)
               
    }

    ## *** k
    if(length(index.k)>0){

        transform.k <- switch(transform.k,
                              "none" = .reparametrize.none,
                              "log" = .reparametrize.log,
                              "square" = .reparametrize.square,
                              "logsquare" = .reparametrize.logsquare,
                              "sd" = .reparametrize.sd,
                              "logsd" = .reparametrize.logsd,
                              "var" = .reparametrize.var,
                              "logvar" = .reparametrize.logvar)

        if(!transform.names){
            level.k <- NULL
        }else{
            level.k <- sapply(level,"[[",1)
        }
        out <- transform.k(p = p, out = out,
                           index = index.k, indexSigma = match(sigma[index.k], name.p),
                           level = level.k, type = "k",
                           inverse = inverse, transform.names = transform.names, Jacobian = Jacobian, dJacobian = dJacobian)
               

    }

    ## *** rho
    if(length(index.rho)>0){
        transform.rho <- switch(transform.rho,
                                "none" = .reparametrize.none,
                                "atanh" = .reparametrize.atanh,
                                "cov" = .reparametrize.cov)

        if(!transform.names){
            level.rho <- NULL
        }else{
            level.rho <- sapply(level,"[[",1)
        }
        out <- transform.rho(p = p, out = out,
                             index = index.rho, indexSigma = match(sigma[index.rho], name.p), indexKx = match(k.x[index.rho], name.p), indexKy = match(k.y[index.rho], name.p),
                             level = level.rho, type = "rho",
                             inverse = inverse, transform.names = transform.names, Jacobian = Jacobian, dJacobian = dJacobian)
    }     

    ## ** export
    return(out)
}


## * .reparametrize.none
.reparametrize.none <- function(p, out, index, indexSigma, indexKx, indexKy, level, inverse, type,
                                transform.names, Jacobian, dJacobian){

    if(!inverse && transform.names && all(!is.na(level[index]))){
        out$newname[index] <- paste0(type, level[index])
    }

    return(out)

}

## * .reparametrize.one
.reparametrize.one <- function(p, out, index, indexSigma, indexKx, indexKy, level, inverse, type,
                               transform.names, Jacobian, dJacobian){

    out$p[index] <- 1

    if(transform.names){
        out$newname[index] <- "reference"
    }

    if(Jacobian){
        diag(out$Jacobian)[index] <- 0
    }

    return(out)

}

## * .reparametrize.log
.reparametrize.log <- function(p, out, index, indexSigma, indexKx, indexKy, level, inverse, type,
                               transform.names, Jacobian, dJacobian){

    if(inverse){

        ## f^{-1}:y -> x=exp(y)
        out$p[index] <- exp(p[index])

    }else{
           
        ## f:x -> y=log(x)
        out$p[index] <- log(p[index])
            
        if(transform.names){
            if(all(is.na(level[index]))){
                out$newname[index] <- paste0("log(",type,")")
            }else{
                out$newname[index] <- paste0("log(",type,")", level[index])
            }
        }
            
        if(Jacobian){
            ## df^{-1}:y -> exp(y)=x
            diag(out$Jacobian)[index] <- p[index] 
        }

        if(dJacobian==1){ ## orignal scale after tranformation
            for(iIndex in index){
                ## d2f^{-1}:y -> exp(y)=x
                out$dJacobian[iIndex,iIndex,iIndex] <- p[iIndex]
            }
        }else if(dJacobian==2){ ## newscale after transformation
            for(iIndex in index){
                ## d/dx df^{-1}:y -> 1
                out$dJacobian[iIndex,iIndex,iIndex] <- 1
            }
        }
    }

    return(out)

}

## * .reparametrize.square
.reparametrize.square <- function(p, out, index, indexSigma, indexKx, indexKy, level, inverse, type,
                                  transform.names, Jacobian, dJacobian){
    if(inverse){

        ## f^{-1}:y -> x=sqrt(y)
        out$p[index] <- sqrt(p[index])

    }else{

        ## f:x -> y=x^2
        out$p[index] <- p[index]^2
                
        if(transform.names){
            if(all(is.na(level[index]))){
                out$newname[index] <- paste0(type, "^2")
            }else{
                out$newname[index] <- paste0(type, "^2", level[index])
            }
        }

        if(Jacobian){
            ## df^{-1}:y -> 1/(2*sqrt(y))=1/(2x)
            diag(out$Jacobian)[index] <- 1/(2*p[index]) 
        }

        if(dJacobian == 1){
            for(iIndex in index){
                ## d2f^{-1}:y -> -1/(4*y^(3/2))=-1/(4x^3)
                out$dJacobian[iIndex,iIndex,iIndex] <- -1/(4*p[iIndex]^3)
            }
        }else if(dJacobian == 2){
            for(iIndex in index){
                ## d/dx df^{-1}:y -> -1/(2x^2)
                out$dJacobian[iIndex,iIndex,iIndex] <- -1/(2*p[iIndex]^2)
            }
        }
    }

    return(out)
}

## * .reparametrize.logsquare
.reparametrize.logsquare <- function(p, out, index, indexSigma, indexKx, indexKy, level, inverse, type,
                                     transform.names, Jacobian, dJacobian){
    if(inverse){

        ## f^{-1}:y -> x=sqrt(exp(y))=exp(y/2)
        out$p[index] <- exp(p[index]/2)

    }else{

        ## f:x -> y=log(x^2)=2log(x)
        out$p[index] <- log(p[index]^2)

        if(transform.names){
            if(all(is.na(level[index]))){
                out$newname[index] <- paste0("log(",type,"^2)")
            }else{
                out$newname[index] <- paste0("log(",type,"^2)", level[index])
            }
        }

        if(Jacobian){
            ## df^{-1}:y -> x=d exp(y/2)= exp(y/2)/2 = x/2
            diag(out$Jacobian)[index] <- p[index]/2 
        }

        if(dJacobian == 1){
            for(iIndex in index){
                ## d2f^{-1}:y -> exp(y/2)/4 = x/4
                out$dJacobian[iIndex,iIndex,iIndex] <- p[iIndex]/4
            }
        }else if(dJacobian == 2){
            for(iIndex in index){
                ## d/dx df^{-1}:y -> 1/2
                out$dJacobian[iIndex,iIndex,iIndex] <- 1/2
            }
        }
    }

    return(out)
}

## * .reparametrize.sd
.reparametrize.sd <- function(p, out, index, indexSigma, indexKx, indexKy, level, inverse, type,
                              transform.names, Jacobian, dJacobian){

    if(type!="k"){
        stop("Only implemented for type \'k\' parameters. \n")
    }
    
    n.index <- length(index)

    if(inverse){
        ## f^{-1}:(y1,y2) -> (x1,x2)=(y1,y2/y1)
        out$p[index] <- p[index]/p[indexSigma]

    }else{

        ## f:(x1,x2) -> (y1,y2)=(x1,x1*x2)
        out$p[index] <- p[index]*p[indexSigma]
         
        if(transform.names){
            out$newname[index] <- paste0("sigma", level[index])
        }

        if(Jacobian){
            ## df^{-1}:(y1,y2) -> (1,-y2/y1^2)=  (1,-x2/x1)            /y1 
            ##                 -> (0,1/y1)    =  (0,1/x1)              /y2
            ## Here we are only interested in the second column, i.e. (-x2/x1,1/x1)
            for(iK in 1:n.index){
                iIndex.k <- index[iK]
                iIndex.sigma <- indexSigma[iK]
                out$Jacobian[iIndex.k,iIndex.sigma] <- -p[iIndex.k]/p[iIndex.sigma]
            }
            diag(out$Jacobian)[index] <- 1/p[indexSigma]
        }
    }
    
    if(dJacobian == 1){
        for(iK in 1:n.index){
            ## df^{-1}:(y1,y2) -> (-y2/y1^2,1/y1)     
            ## d2f^{-1}:(y1,y2) -> (2 y2/y1^3,-1/y1^2) = (2x2/x1^2,-1/x1^1)    
            ##                  -> (-1/y1^2,0)         = (-1/x^2,0)
            iIndex.k <- index[iK]
            iIndex.sigma <- indexSigma[iK]
            out$dJacobian[iIndex.k,iIndex.sigma,iIndex.sigma] <- 2*p[iIndex.k]/(p[iIndex.sigma]^2)
            out$dJacobian[iIndex.k,iIndex.k,iIndex.sigma] <- -1/p[iIndex.sigma]^2
            out$dJacobian[iIndex.k,iIndex.sigma,iIndex.k] <- -1/p[iIndex.sigma]^2
        }
    }else if(dJacobian == 2){
        for(iK in 1:n.index){
            ## df^{-1}:(y1,y2) -> (-x2/x1,1/x1)     
            ## d/dx df^{-1}:(y1,y2) -> (x2/x1^2,-1/x1^2)
            ##                      -> (-1/x1,0)
            iIndex.k <- index[iK]
            iIndex.sigma <- indexSigma[iK]
            out$dJacobian[iIndex.k,iIndex.sigma,iIndex.sigma] <- p[iIndex.k]/p[iIndex.sigma]^2
            out$dJacobian[iIndex.k,iIndex.k,iIndex.sigma] <- -1/p[iIndex.sigma]^2
            out$dJacobian[iIndex.k,iIndex.sigma,iIndex.k] <- -1/p[iIndex.sigma]
        }
    }

    return(out)
}

## * .reparametrize.logsd
.reparametrize.logsd <- function(p, out, index, indexSigma, indexKx, indexKy, level, inverse, type,
                                 transform.names, Jacobian, dJacobian){

    if(type!="k"){
        stop("Only implemented for type \'k\' parameters. \n")
    }
    
    n.index <- length(index)

    if(inverse){

        ## f^{-1}:(y1,y2) -> (x1,x2)=(exp(y1),exp(y2-y1))
        out$p[index] <- exp(p[index]-p[indexSigma])

    }else{

        ## f:(x1,x2) -> (y1,y2)=(log(x1),log(x1*x2))
        out$p[index] <- log(p[index]*p[indexSigma])

        if(transform.names){
            out$newname[index] <- paste0("log(sigma)", level[index])
        }

        if(Jacobian){
            ## df^{-1}:(y1,y2) -> (exp(y1),-exp(y2-y1)) = (x1, -x2)            /y1
            ##                 -> (0,exp(y2-y1)) = (0, x2)                     /y2
            ## Here we are only interested in the second column, i.e. (-x2,x2)
            for(iK in 1:n.index){
                iIndex.k <- index[iK]
                iIndex.sigma <- indexSigma[iK]
                out$Jacobian[iIndex.k,iIndex.sigma] <- -p[iIndex.k]
            }
            diag(out$Jacobian)[index] <- p[index]
        }
        if(dJacobian == 1){
            for(iK in 1:n.index){
                ## df^{-1}:(y1,y2) -> (-exp(y2-y1),exp(y2-y1))
                ## d2f^{-1}:(y1,y2) -> (exp(y2-y1),-exp(y2-y1)) = (x2,-x2)
                ##                  -> (-exp(y2-y1),exp(y2-y1)) = (-x2,x2)
                iIndex.k <- index[iK]
                iIndex.sigma <- indexSigma[iK]
                out$dJacobian[iIndex.k,iIndex.sigma,iIndex.sigma] <- p[iIndex.k]
                out$dJacobian[iIndex.k,iIndex.k,iIndex.sigma] <- -p[iIndex.k]
                out$dJacobian[iIndex.k,iIndex.sigma,iIndex.k] <- -p[iIndex.k]
                out$dJacobian[iIndex.k,iIndex.k,iIndex.k] <- p[iIndex.k]
            }
        }else if(dJacobian == 2){
            for(iK in 1:n.index){
                ## df^{-1}:(y1,y2) -> (-x2,x2)
                ## d/dx df^{-1}:(y1,y2) -> (0,0)
                ##                      -> (-1,1)
                iIndex.k <- index[iK]
                iIndex.sigma <- indexSigma[iK]
                out$dJacobian[iIndex.k,iIndex.sigma,iIndex.k] <- -1
                out$dJacobian[iIndex.k,iIndex.k,iIndex.k] <- 1
            }                            
        }

    }

    return(out)
}

## * .reparametrize.var
.reparametrize.var <- function(p, out, index, indexSigma, indexKx, indexKy, level, inverse, type,
                               transform.names, Jacobian, dJacobian){

    if(type!="k"){
        stop("Only implemented for type \'k\' parameters. \n")
    }
    
    n.index <- length(index)

    if(inverse){

        ## f^{-1}:(y1,y2) -> (x1,x2)=(sqrt(y1),sqrt(y2/y1))
        out$p[index] <- sqrt(p[index]/p[indexSigma])

    }else{

        ## f:(x1,x2) -> (y1,y2)=(x1^2,x1^2*x2^2)
        out$p[index] <- p[indexSigma]^2*p[index]^2

        if(transform.names){
            out$newname[index] <- paste0("sigma^2", level[index])
        }

        if(Jacobian){
            ## df^{-1}:(y1,y2) -> (1/(2*sqrt(y1)),-sqrt(y2)/(2*y1^(3/2))=  (1/(2*x1),-x2/(2*x1^2))            /y1 
            ##                 -> (0,1/(2*sqrt(y2*y1)))              =  (0,1/(2*x2*x1^2))                     /y2
            ## Here we are only interested in the second column, i.e. (-x2/(2*x1^2),1/(2*x1^2*x2))
            for(iK in 1:n.index){
                iIndex.k <- index[iK]
                iIndex.sigma <- indexSigma[iK]
                out$Jacobian[iIndex.k,iIndex.sigma] <- -p[iIndex.k]/(2*p[iIndex.sigma]^2)
            }
            diag(out$Jacobian)[index] <- 1/(2*p[indexSigma]^2*p[index])
            
        }
        if(dJacobian == 1){
            ## sqrt(\sigma1), sqrt(\sigma2/\sigma1), ...   so  -1/(4 \sigma1^(3/2)) 3\sqrt(\sigma2)/(4 \sigma1^(5/2))      0 -1/(4 \sqrt(\sigma2) \sigma1^(3/2))        
            ##                                                      0              - 1 / (4*\sqrt(\sigma2)\sigma1^(3/2))   0 -1 / (4*\sigma2^(3/2)\sqrt(\sigma1))
            ## WARNING we plug-in sqrt(\sigma1) and sqrt(\sigma2/\sigma1)
            for(iK in 1:n.index){
                ## df^{-1}:(y1,y2) -> (-sqrt(y2)/(2*y1^(3/2)),1/(2*sqrt(y2*y1)))
                ## d2f^{-1}:(y1,y2) -> (3sqrt(y2)/(4*y1^(5/2)),-1/(4*sqrt(y2)*y1^(3/2)))) = (3 x2/(4x1^4), -1/(4*x2*x1^4))
                ##                  -> (-1/(4*sqrt(y2)*y1^(3/2)),-1/(4*y2^(3/2)*y1))) = (-1/(4*x2*x1^4), -1/(4*x2^3*x1^4))
                iIndex.k <- index[iK]
                iIndex.sigma <- indexSigma[iK]
                out$dJacobian[iIndex.k,iIndex.sigma,iIndex.sigma] <- 3*p[iIndex.k]/(4*p[iIndex.sigma]^4)
                out$dJacobian[iIndex.k,iIndex.sigma,iIndex.k] <- -1/(4*p[iIndex.sigma]^4*p[iIndex.k])
                out$dJacobian[iIndex.k,iIndex.k,iIndex.sigma] <- -1/(4*p[iIndex.sigma]^4*p[iIndex.k])
                out$dJacobian[iIndex.k,iIndex.k,iIndex.k] <- -1/(4*p[iIndex.sigma]^4*p[iIndex.k]^3)
            }
        }else if(dJacobian == 2){
            for(iK in 1:n.index){
                ## df^{-1}:(y1,y2) -> (-x2/(2*x1^2),1/(2*x1^2*x2))
                ## d/dx df^{-1}:(y1,y2) -> (x2/(4*x1^3),-1/(4*x1^3*x2))
                ##                      -> (-1/(2*x1^2),-1/(2*x1^2*x2^2))
                iIndex.k <- index[iK]
                iIndex.sigma <- indexSigma[iK]
                out$dJacobian[iIndex.k,iIndex.sigma,iIndex.sigma] <- p[iIndex.k]/(4*p[iIndex.sigma]^3)
                out$dJacobian[iIndex.k,iIndex.k,iIndex.sigma] <- -1/(4*p[iIndex.sigma]^3*p[iIndex.k])
                out$dJacobian[iIndex.k,iIndex.sigma,iIndex.k] <- -1/(2*p[iIndex.sigma]^2)
                out$dJacobian[iIndex.k,iIndex.k,iIndex.k] <- -1/(2*p[iIndex.sigma]^2*p[iIndex.k]^2)
            }
        }
    }

    return(out)
}

## * .reparametrize.logvar
.reparametrize.logvar <- function(p, out, index, indexSigma, indexKx, indexKy, level, inverse, type,
                                  transform.names, Jacobian, dJacobian){
    if(type!="k"){
        stop("Only implemented for type \'k\' parameters. \n")
    }
    
    n.index <- length(index)

    if(inverse){

        ## f^{-1}:(y1,y2) -> (x1,x2)=(sqrt(exp(y1)),sqrt(exp(y2-y1)))=(exp(y1/2),exp(y2/2-y1/2))
        out$p[index] <- exp(p[index]/2-p[indexSigma]/2)

    }else{

        ## f:(x1,x2) -> (y1,y2)=(log(x1^2),log(x1^2*x2^2))=(2*log(x1),2*log(x1*x2))
        out$p[index] <- 2*log(p[index]*p[indexSigma])

        if(transform.names){
            out$newname[index] <- paste0("log(sigma^2)", level[index])
        }
        
        if(Jacobian){
            ## df^{-1}:(y1,y2) -> (exp(y1/2)/2,-exp(y2/2-y1/2)/2) = (x1/2, -x2/2)            /y1
            ##                 -> (0,exp(y2/2-y1/2)/2) = (0, x2*2)                          /y2
            ## Here we are only interested in the second column, i.e. (-x2/2/,x2/2)
            for(iK in 1:n.index){
                iIndex.k <- index[iK]
                iIndex.sigma <- indexSigma[iK]
                out$Jacobian[iIndex.k,iIndex.sigma] <- -p[iIndex.k]/2
            }
            diag(out$Jacobian)[index] <- p[index]/2
        }
        if(dJacobian == 1){
            ## df^{-1}:(y1,y2) -> (-exp(y2/2-y1/2)/2,exp(y2/2-y1/2)/2)
            ## d2f^{-1}:(y1,y2) -> (exp(y2/2-y1/2)/4,-exp(y2/2-y1/2)/4) = (x2/4,-x2/4)
            ##                  -> (-exp(y2/2-y1/2)/4,exp(y2/2-y1/2)/4) = (-x2/4,x2/4)
            for(iK in 1:n.index){
                iIndex.k <- index[iK]
                iIndex.sigma <- indexSigma[iK]
                out$dJacobian[iIndex.k,iIndex.sigma,iIndex.sigma] <- p[iIndex.k]/4
                out$dJacobian[iIndex.k,iIndex.sigma,iIndex.k] <- -p[iIndex.k]/4
                out$dJacobian[iIndex.k,iIndex.k,iIndex.sigma] <- -p[iIndex.k]/4
                out$dJacobian[iIndex.k,iIndex.k,iIndex.k] <- p[iIndex.k]/4
            }
        }else if(dJacobian == 2){
            for(iK in 1:n.index){
                ## df^{-1}:(y1,y2) -> (-x2/2/,x2/2)
                ## d/dx df^{-1}:(y1,y2) -> (0,0)
                ##                      -> (-1/2,1/2)
                iIndex.k <- index[iK]
                iIndex.sigma <- indexSigma[iK]
                out$dJacobian[iIndex.k,iIndex.sigma,iIndex.k] <- -1/2
                out$dJacobian[iIndex.k,iIndex.k,iIndex.k] <- 1/2
            }            
        }

    }
    return(out)                
}

## * .reparametrize.atanh
.reparametrize.atanh <- function(p, out, index, indexSigma, indexKx, indexKy, level, inverse, type,
                                 transform.names, Jacobian, dJacobian){
    if(type!="rho"){
        stop("Only implemented for type \'rho\' parameters. \n")
    }
    
    n.index <- length(index)

    if(inverse){

        ## f^{-1}:y -> x=tanh(y)
        out$p[index] <- tanh(p[index])

    }else{

        ## f:x -> y=atanh(x)
        out$p[index] <- atanh(p[index])

        if(transform.names){
            out$newname[index] <- paste0("atanh(rho)", level[index])
        }

        if(Jacobian){
            ## df^{-1}:y -> d tanh(y) = 1-tanh(y)^2 = 1-x^2
            diag(out$Jacobian)[index] <- (1-p[index]^2) ## (1-tanh(x)^2) where x=atanh(rho)
        }
        if(dJacobian == 1){
            for(iRho in index){
                ## d2f^{-1}:y -> d tanh(y) = - 2 tanh(y)(1-tanh(y)^2) = -2 x (1-x^2)
                out$dJacobian[iRho,iRho,iRho] <- -2*p[iRho]*(1-p[iRho]^2)
            }
        }else if(dJacobian == 2){
            for(iRho in index){
                ## d/dx df^{-1}:y -> -2x
                out$dJacobian[iRho,iRho,iRho] <- -2*p[iRho]
            }
        }
    }    
    
    return(out)
}

## * .reparametrize.cov
.reparametrize.cov <- function(p, out, index, indexSigma, indexKx, indexKy, level, inverse, type,
                               transform.names, Jacobian, dJacobian){

    if(type!="rho"){
        stop("Only implemented for type \'rho\' parameters. \n")
    }
    
    n.index <- length(index)

    Sigma <- p[indexSigma]
    Kx <- p[indexKx]
    Ky <- p[indexKy]
        
    if(inverse){

        ## fa^{-1}:(y1,y4) -> (x1,x4)=(sqrt(y1),y4/y1)
        ## fb^{-1}:(y1,y2,y3) -> (x1,x2,x3)=(sqrt(y1),sqrt(y2/y1),y4/sqrt(y1*y2))
        ## fc^{-1}:(y1,y2,y3) -> (x1,x2,x3)=(sqrt(y1),sqrt(y2/y1),y4/y2)
        ## fd^{-1}:(y1,y2,y3,y4) -> (x1,x2,x3,x4)=(sqrt(y1),sqrt(y2/y1),sqrt(y3/y1),y4/sqrt(y2*y3))
        Kx[is.na(Kx)] <- Sigma[is.na(Kx)]
        Ky[is.na(Ky)] <- Sigma[is.na(Ky)]
        out$p[index] <- p[index]/sqrt(Kx*Ky)

    }else{

        ## fa:(x1,x4) -> (y1,y4)=(x1^2,x1^2*x4)
        ## fb:(x1,x2,x4) -> (y1,y2,y4)=(x1^2,x1^2*x2^2,x1^2*x2*x4)
        ## fc:(x1,x2,x4) -> (y1,y2,y4)=(x1^2,x1^2*x2^2,x1^2*x2^2*x4)
        ## fd:(x1,x2,x3,x4) -> (y1,y2,y3,y4)=(x1^2,x1^2*x2^2,x1^2*x3^2,x1^2*x2*x3*x4)
        Kx[is.na(Kx)] <- 1
        Ky[is.na(Ky)] <- 1
        out$p[index] <- p[index]*Sigma^2*Kx*Ky

        if(transform.names){
            out$newname[index] <- paste0("cov", level[index])
        }

        if(Jacobian){
            ## dfa^{-1}:(y1,y4) -> (1/(2*sqrt(y1)),-y4/y1^2)=  (1/(2*x1),-x4/x1^2)          /y1 
            ##                  -> (0,1/y1)                 =  (0,1/x1^2)                   /y4
            ## Here we are only interested in the second column, i.e. (-x4/(x1^2),1/(x1^2))

            ## dfb^{-1}:(y1,y2,y4) -> -y4/(2*sqrt(y2)*y1^(3/2)) =  -x4 /(2*y1) = -x4 / (2*x1^2)      /y1 
            ##                     -> -y4/(2*sqrt(y1)*y2^(3/2)) =  -x4 /(2*y2) = -x4 / (2*x1^2*x2^2) /y2
            ##                     -> 1/sqrt(y1*y2)             =  1/(x1^2*x2)                       /y4

            ## dfc^{-1}:(y1,y2,y4) -> 0                                        /y1 
            ##                     -> -y4/y2^2 = -x4/y2 = -x4/(x1^2*x2^2)      /y2
            ##                     -> 1/y2              =  1/(x1^2*x2^2)       /y4

            ## dfd^{-1}:(y1,y2,y3,y4) -> 0                                                            /y1 
            ##                        -> -y4/(2*sqrt(y3)*y2^(3/2)) = -x4/(2*y2) = -x4/(2*x1^2*x2^2)   /y2
            ##                        -> -y4/(2*sqrt(y2)*y3^(3/2)) = -x4/(2*y3) = -x4/(2*x1^2*x3^2)   /y3
            ##                        -> 1/sqrt(y2*y3)     =  1/(x1^2*x2*x3)                          /y4

            for(iRho in 1:n.index){ ## iRho <- 1
                iIndex.sigma <- indexSigma[iRho]
                iIndex.kx <- indexKx[iRho]
                iIndex.ky <- indexKy[iRho]
                iIndex.rho <- index[iRho]
                if(is.na(iIndex.kx) && is.na(iIndex.ky)){
                    out$Jacobian[iIndex.rho,iIndex.sigma] <- -p[iIndex.rho]/p[iIndex.sigma]^2 ## /y1
                }else if(!is.na(iIndex.kx) && is.na(iIndex.ky)){
                    out$Jacobian[iIndex.rho,iIndex.sigma] <- -p[iIndex.rho]/(2*p[iIndex.sigma]^2) ## /y1                    
                    out$Jacobian[iIndex.rho,iIndex.kx] <- -p[iIndex.rho]/(2*p[iIndex.sigma]^2*p[iIndex.kx]^2) ## /y2                    
                }else if(is.na(iIndex.kx) && !is.na(iIndex.ky)){
                    out$Jacobian[iIndex.rho,iIndex.sigma] <- -p[iIndex.rho]/(2*p[iIndex.sigma]^2) ## /y1                    
                    out$Jacobian[iIndex.rho,iIndex.ky] <- -p[iIndex.rho]/(2*p[iIndex.sigma]^2*p[iIndex.ky]^2) ## /y2                   
                }else{ ## both non-NA
                    ## 0 for /y1
                    out$Jacobian[iIndex.rho,iIndex.kx] <- out$Jacobian[iIndex.rho,iIndex.sigma]-p[iIndex.rho]/(2*p[iIndex.sigma]^2*p[iIndex.kx]^2) ## /y3
                    out$Jacobian[iIndex.rho,iIndex.ky] <- out$Jacobian[iIndex.rho,iIndex.sigma]-p[iIndex.rho]/(2*p[iIndex.sigma]^2*p[iIndex.ky]^2) ## /y4                   
                }
            }
            diag(out$Jacobian)[index] <- 1/(Sigma^2*Kx*Ky) ## /y4 ok
            
        }
        if(dJacobian == 1){
            ##    fa^{-1}:(y1,y4) -> (x1,x4)=(sqrt(y1),y4/y1)
            ##    dfa^{-1}:(y1,y4) -> (-y4/y1^2,1/y1)
            ## so d2fa^{-1}:(y1,y4) -> (2 y4/y1^3,-1/y1^2)   =  (2 x4/x1^4,-1/x1^4)     // y1
            ##                      -> (-1/y1^2,0)           = (-1/x1^4,0)             // y4

            ##    fb^{-1}:(y1,y2,y3) -> (x1,x2,x3)=(sqrt(y1),sqrt(y2/y1),y4/sqrt(y1*y2))
            ##    dfb^{-1}:(y1,y2,y4) -> ( -y4/(2*sqrt(y2)*y1^(3/2)),-y4/(2*sqrt(y1)*y2^(3/2)),1/sqrt(y1*y2) )
            ## so d2fb^{-1}:(y1,y2,y4) -> (3*y4/(4*sqrt(y2)*y1^(5/2)), y4/(4*y1^(3/2)*y2^(3/2)), -1/(2*y1^(3/2)*sqrt(y2))) = (3*x4/(4*x1^4)   ,x4/(4*x1^4*x2^2) ,-1/(2*x1^4*x2)) ## y1
            ##                         -> (y4/(4*y1^(3/2)*y2^(3/2)), 3*y4/(4*sqrt(y1)*y2^(5/2)), -1/(2*sqrt(y1)*y2^(3/2))) = (x4/(4*x1^4*x2^2),3*x4/(4*x1^4*x2^4),-1/(2*x1^4*x2^3)) ## y2
            ##                         -> (-1/(2*sqrt(y2)*y1^(3/2)),-1/(2*sqrt(y1)*y2^(3/2)) , 0)                          = (-1/(2*x1^4*x2)  , -1/(2*x1^4*x2^3), 0)     ## y4
        
            ## fd^{-1}:(y1,y2,y3,y4) -> (x1,x2,x3,x4)=(sqrt(y1),sqrt(y2/y1),sqrt(y3/y1),y4/sqrt(y2*y3))
            ## dfd^{-1}:(y1,y2,y3,y4) -> (0, -y4/(2*sqrt(y3)*y2^(3/2)), -y4/(2*sqrt(y2)*y3^(3/2)), 1/sqrt(y2*y3))
            ## so d2fd^{-1}:(y1,y2,y3,y4) -> (0, 0, 0, 0)
            ##                            -> (0, 3*y4/(4*sqrt(y3)*y2^(5/2)), y4/(4*y2^(3/2)*y3^(3/2)), -1/(2y2^(3/2)*sqrt(y3))) = (0,3*x4/(4*x2^4*x1^4),x4/(4*x2^2*x3^2*x1^4),-1/(2*x2^3*x3*x1^4)))
            ##                            -> (0, y4/(4*y3^(3/2)*y2^(3/2)), 3y4/(4*sqrt(y2)*y3^(5/2)), -1/(2sqrt(y2)*y3^(3/2))) = (0,x4/(4*x3^2^x2^2*x1^4),3*x4/(4*y3^4*x1^4),-1/(2*x2*x3^3*x1^4))
            ##                            -> (0, -1/(2*sqrt(y3)*y2^(3/2)), -1/(2*sqrt(y2)*y3^(3/2)), 0) = (0,-1/(2*x3*x2^3*x1^4),-1/(2*x2*x3^3*x1^4),0)

            for(iRho in 1:n.index){ ## iRho <- 1
                iIndex.sigma <- indexSigma[iRho]
                iIndex.kx <- indexKx[iRho]
                iIndex.ky <- indexKy[iRho]
                iIndex.rho <- index[iRho]
                if(is.na(iIndex.kx) && is.na(iIndex.ky)){

                    out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.sigma] <- 2*p[iIndex.rho]/p[iIndex.sigma]^4
                    out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.rho] <- -1/p[iIndex.sigma]^4
                    out$dJacobian[iIndex.rho,iIndex.rho,iIndex.sigma] <- -1/p[iIndex.sigma]^4

                }else if(!is.na(iIndex.kx) && is.na(iIndex.ky)){
                    
                    out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.sigma] <- 3*p[iIndex.rho]/(4*p[iIndex.sigma]^4)               
                    out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.kx] <- out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.kx] + p[iIndex.rho]/(4*p[iIndex.sigma]^4*p[iIndex.kx]^2)               
                    out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.rho] <- -1/(2*p[iIndex.sigma]^4*p[iIndex.kx])               

                    out$dJacobian[iIndex.rho,iIndex.kx,iIndex.sigma] <- out$dJacobian[iIndex.rho,iIndex.kx,iIndex.sigma] + p[iIndex.rho]/(4*p[iIndex.sigma]^4*p[iIndex.kx]^2)
                    out$dJacobian[iIndex.rho,iIndex.kx,iIndex.kx] <- out$dJacobian[iIndex.rho,iIndex.kx,iIndex.kx] + 3*p[iIndex.rho]/(4*p[iIndex.sigma]^4*p[iIndex.kx]^4)
                    out$dJacobian[iIndex.rho,iIndex.kx,iIndex.rho] <- out$dJacobian[iIndex.rho,iIndex.kx,iIndex.rho] - 1/(2*p[iIndex.sigma]^4*p[iIndex.kx]^3)

                    out$dJacobian[iIndex.rho,iIndex.rho,iIndex.sigma] <- -1/(2*p[iIndex.sigma]^4*p[iIndex.kx])
                    out$dJacobian[iIndex.rho,iIndex.rho,iIndex.kx] <- out$dJacobian[iIndex.rho,iIndex.rho,iIndex.kx] - 1/(2*p[iIndex.sigma]^4*p[iIndex.kx]^3)


                }else if(is.na(iIndex.kx) && !is.na(iIndex.ky)){

                    out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.sigma] <- 3*p[iIndex.rho]/(4*p[iIndex.sigma]^4)               
                    out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.ky] <- out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.ky] + p[iIndex.rho]/(4*p[iIndex.sigma]^4*p[iIndex.ky]^2)               
                    out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.rho] <- -1/(2*p[iIndex.sigma]^4*p[iIndex.ky])               

                    out$dJacobian[iIndex.rho,iIndex.ky,iIndex.sigma] <- out$dJacobian[iIndex.rho,iIndex.ky,iIndex.sigma] + p[iIndex.rho]/(4*p[iIndex.sigma]^4*p[iIndex.ky]^2)
                    out$dJacobian[iIndex.rho,iIndex.ky,iIndex.ky] <- out$dJacobian[iIndex.rho,iIndex.ky,iIndex.ky] + 3*p[iIndex.rho]/(4*p[iIndex.sigma]^4*p[iIndex.ky]^4)
                    out$dJacobian[iIndex.rho,iIndex.ky,iIndex.rho] <- out$dJacobian[iIndex.rho,iIndex.ky,iIndex.rho] - 1/(2*p[iIndex.sigma]^4*p[iIndex.ky]^3)

                    out$dJacobian[iIndex.rho,iIndex.rho,iIndex.sigma] <- -1/(2*p[iIndex.sigma]^4*p[iIndex.ky])
                    out$dJacobian[iIndex.rho,iIndex.rho,iIndex.ky] <- out$dJacobian[iIndex.rho,iIndex.rho,iIndex.ky] - 1/(2*p[iIndex.sigma]^4*p[iIndex.ky]^3)

                }else{ ## both non-NA
                    ## 0 for /y1
                    out$dJacobian[iIndex.rho,iIndex.kx,iIndex.kx] <- out$dJacobian[iIndex.rho,iIndex.kx,iIndex.kx]+3*p[iIndex.rho]/(4*p[iIndex.sigma]^4*p[iIndex.kx]^4)
                    out$dJacobian[iIndex.rho,iIndex.kx,iIndex.ky] <- out$dJacobian[iIndex.rho,iIndex.kx,iIndex.ky]+p[iIndex.rho]/(4*p[iIndex.sigma]^4*p[iIndex.kx]^2*p[iIndex.ky]^2)
                    out$dJacobian[iIndex.rho,iIndex.kx,iIndex.rho] <- out$dJacobian[iIndex.rho,iIndex.kx,iIndex.rho]-1/(2*p[iIndex.sigma]^4*p[iIndex.kx]^3*p[iIndex.ky])

                    out$dJacobian[iIndex.rho,iIndex.ky,iIndex.kx] <- out$dJacobian[iIndex.rho,iIndex.ky,iIndex.kx]+p[iIndex.rho]/(4*p[iIndex.sigma]^4*p[iIndex.kx]^2*p[iIndex.ky]^2)
                    out$dJacobian[iIndex.rho,iIndex.ky,iIndex.ky] <- out$dJacobian[iIndex.rho,iIndex.ky,iIndex.ky]+3*p[iIndex.rho]/(4*p[iIndex.sigma]^4*p[iIndex.ky]^4)
                    out$dJacobian[iIndex.rho,iIndex.ky,iIndex.rho] <- out$dJacobian[iIndex.rho,iIndex.ky,iIndex.rho]-1/(2*p[iIndex.sigma]^4*p[iIndex.kx]*p[iIndex.ky]^3)

                    out$dJacobian[iIndex.rho,iIndex.rho,iIndex.kx] <- out$dJacobian[iIndex.rho,iIndex.rho,iIndex.kx]-1/(2*p[iIndex.sigma]^4*p[iIndex.kx]^3*p[iIndex.ky]) 
                    out$dJacobian[iIndex.rho,iIndex.rho,iIndex.ky] <- out$dJacobian[iIndex.rho,iIndex.rho,iIndex.ky]-1/(2*p[iIndex.sigma]^4*p[iIndex.kx]*p[iIndex.ky]^3) 
                }
                
            }
            
            
        }else if(dJacobian == 2){

            for(iRho in 1:n.index){
                ## dfa^{-1}:(y1,y4) -> (-x4/(x1^2),1/(x1^2))
                ## d/dx dfa^{-1}:(y1,y4) -> (2*x4/(x1^3),-2/(x1^3))
                ##                       -> (-1/(x1^2),0)

                ## dfb^{-1}:(y1,y2,y4) -> (-x4 / (2*x1^2), -x4 / (2*x1^2*x2^2), 1/(x1^2*x2))
                ## d/dx dfb^{-1}:(y1,y2,y4) -> (x4 / (x1^3),   x4 / (x1^3*x2^2),   -2/(x1^3*x2))
                ##                          -> (0,             x4 / (x1^2*x2^3),   -1/(x1^2*x2^2))
                ##                          -> (-1 / (2*x1^2), -1 / (2*x1^2*x2^2),  0)
                
                ## dfd^{-1}:(y1,y2,y3,y4) -> (0, -x4/(2*x1^2*x2^2), -x4/(2*x1^2*x3^2), 1/(x1^2*x2*x3))
                ## d/dx dfd^{-1}:(y1,y2,y3,y4) -> (0, x4/(x1^3*x2^2),   x4/(x1^3*x3^2),  -2/(x1^3*x2*x3))
                ##                             -> (0, x4/(x1^2*x2^3),   0,               -1/(x1^2*x2^2*x3))
                ##                             -> (0, 0,                x4/(x1^2*x3^3),  -1/(x1^2*x2*x3^2))
                ##                             -> (0, -1/(2*x1^2*x2^2), -1/(2*x1^2*x3^2), 0)
                for(iRho in 1:n.index){ ## iRho <- 1
                    iIndex.sigma <- indexSigma[iRho]
                    iIndex.kx <- indexKx[iRho]
                    iIndex.ky <- indexKy[iRho]
                    iIndex.rho <- index[iRho]
                    if(is.na(iIndex.kx) && is.na(iIndex.ky)){
                        out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.sigma] <- 2*p[iIndex.rho]/(p[iIndex.sigma]^3)
                        out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.rho] <- -1/p[iIndex.sigma]^2
                    }else if(!is.na(iIndex.kx) && is.na(iIndex.ky)){
                        out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.sigma] <- p[iIndex.rho]/(p[iIndex.sigma]^3)
                        out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.rho] <- -1/(2*p[iIndex.sigma]^2)
                    
                        out$dJacobian[iIndex.rho,iIndex.kx,iIndex.sigma] <- out$dJacobian[iIndex.rho,iIndex.kx,iIndex.sigma]+p[iIndex.rho]/(p[iIndex.sigma]^3*p[iIndex.kx]^2)
                        out$dJacobian[iIndex.rho,iIndex.kx,iIndex.kx] <- out$dJacobian[iIndex.rho,iIndex.kx,iIndex.kx]+p[iIndex.rho]/(p[iIndex.sigma]^2*p[iIndex.kx]^3)
                        out$dJacobian[iIndex.rho,iIndex.kx,iIndex.rho] <- out$dJacobian[iIndex.rho,iIndex.kx,iIndex.rho]-1/(2*p[iIndex.sigma]^2*p[iIndex.kx]^2)
                    }else if(is.na(iIndex.kx) && !is.na(iIndex.ky)){
                        out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.sigma] <- p[iIndex.rho]/(p[iIndex.sigma]^3)
                        out$dJacobian[iIndex.rho,iIndex.sigma,iIndex.rho] <- -1/(p[iIndex.sigma]^2)
                   
                        out$dJacobian[iIndex.rho,iIndex.ky,iIndex.sigma] <- out$dJacobian[iIndex.rho,iIndex.ky,iIndex.sigma]+p[iIndex.rho]/(p[iIndex.sigma]^3*p[iIndex.ky]^2)
                        out$dJacobian[iIndex.rho,iIndex.ky,iIndex.ky] <- out$dJacobian[iIndex.rho,iIndex.ky,iIndex.ky]+p[iIndex.rho]/(p[iIndex.sigma]^2*p[iIndex.ky]^3)
                        out$dJacobian[iIndex.rho,iIndex.ky,iIndex.rho] <- out$dJacobian[iIndex.rho,iIndex.ky,iIndex.rho]-1/(2*p[iIndex.sigma]^2*p[iIndex.ky]^2)
                    }else{ ## both non-NA
                        out$dJacobian[iIndex.rho,iIndex.kx,iIndex.sigma] <- out$dJacobian[iIndex.rho,iIndex.kx,iIndex.sigma]+p[iIndex.rho]/(p[iIndex.sigma]^3*p[iIndex.kx]^2) 
                        out$dJacobian[iIndex.rho,iIndex.kx,iIndex.ky] <- out$dJacobian[iIndex.rho,iIndex.kx,iIndex.ky]+p[iIndex.rho]/(p[iIndex.sigma]^2*p[iIndex.kx]^3) 
                        out$dJacobian[iIndex.rho,iIndex.kx,iIndex.rho] <- out$dJacobian[iIndex.rho,iIndex.kx,iIndex.rho]-1/(2*p[iIndex.sigma]^2*p[iIndex.kx]^2) 
                        
                        out$dJacobian[iIndex.rho,iIndex.ky,iIndex.sigma] <- out$dJacobian[iIndex.rho,iIndex.ky,iIndex.sigma]+p[iIndex.rho]/(p[iIndex.sigma]^3*p[iIndex.ky]^2)
                        out$dJacobian[iIndex.rho,iIndex.ky,iIndex.ky] <- out$dJacobian[iIndex.rho,iIndex.ky,iIndex.ky]+p[iIndex.rho]/(p[iIndex.sigma]^2*p[iIndex.ky]^3)
                        out$dJacobian[iIndex.rho,iIndex.ky,iIndex.rho] <- out$dJacobian[iIndex.rho,iIndex.ky,iIndex.rho]-1/(2*p[iIndex.sigma]^2*p[iIndex.ky]^2)
                    }
                    
                    out$dJacobian[iIndex.rho,iIndex.rho,iIndex.sigma] <- -2/(p[iIndex.sigma]^3*p[iIndex.kx]*p[iIndex.ky])
                    if(!is.na(iIndex.kx)){
                        out$dJacobian[iIndex.rho,iIndex.rho,iIndex.kx] <- out$dJacobian[iIndex.rho,iIndex.rho,iIndex.kx]-1/(p[iIndex.sigma]^2*p[iIndex.kx]^2*p[iIndex.ky])
                    }
                    if(!is.na(iIndex.ky)){
                        out$dJacobian[iIndex.rho,iIndex.rho,iIndex.ky] <- out$dJacobian[iIndex.rho,iIndex.rho,iIndex.ky]-1/(p[iIndex.sigma]^2*p[iIndex.kx]*p[iIndex.ky]^2)
                    }
                    

                    
                }
            }
        }
    }
  
    return(out)
}

## * .init_transform
##' @description Initalize the transformations for the variance and correlation coefficients.
##' Can also reparametrise the input p to no transformation (used by estimate when the user ask for certain transformations)
##' @param p [numeric vector] value for the model parameters 
##' @param transform.sigma,transform.k,transform.rho [character] user input for the transformations
##' @param x.transform.sigma,x.transform.k,x.transform.rho [character] transformations used when fitting the object
##' @param normalize.p [logical] should the model parameter be re-parametrized to no transformation
##' @noRd
.init_transform <- function(p, transform.sigma, transform.k, transform.rho, 
                            x.transform.sigma, x.transform.k, x.transform.rho,
                            table.param){

    ## ** normalize input
    ## several way to say no transform
    ## do not use identical(,) because transform.sigma/k/rho may contain an attribute
    if(length(transform.sigma==1) && ((transform.sigma == "") || (transform.sigma == "no") || (transform.sigma == FALSE))){
        transform.sigma <- "none"
    }
    if(length(transform.k==1) && ((transform.k == "") || (transform.k == "no") || (transform.k == FALSE))){
        transform.k <- "none"
    }
    if(length(transform.rho==1) && ((transform.rho == "") || (transform.rho == "no") || (transform.rho == FALSE))){
        transform.rho <- "none"
    }

    ## get attributes
    p.transform.sigma <- attr(p,"transform.sigma")
    p.transform.k <- attr(p,"transform.k")
    p.transform.rho <- attr(p,"transform.rho")
        
    ## ** transform
    ## transformation attribute given to p by the estimate function overrule object transformation
    if(is.null(transform.rho)){
        if(!is.null(p.transform.rho)){
            transform.rho <- p.transform.rho
        }else{
            transform.rho <- x.transform.rho
        }
        attr(transform.rho,"arg") <- NULL
    }else{
        transform.rho.save <- transform.rho
        attr(transform.rho.save,"arg") <- NULL
        transform.rho <- match.arg(transform.rho, c("none","atanh", "cov"))
        attr(transform.rho,"arg") <- transform.rho.save
    }

    if(is.null(transform.k)){
        if(!is.null(p.transform.k)){
            transform.k <- p.transform.k
        }else if(transform.rho == "cov"){
            transform.k <- "var"
        }else{
            transform.k <- x.transform.k
        }
        attr(transform.k,"arg") <- NULL
    }else{
        transform.k.save <- transform.k
        attr(transform.k.save,"arg") <- NULL
        transform.k <- match.arg(transform.k, c("none","log","square","logsquare","sd","logsd","var","logvar"))
        attr(transform.k,"arg") <- transform.k.save
    }
    
    if(is.null(transform.sigma)){
        if(!is.null(p.transform.sigma)){
            transform.sigma <- p.transform.sigma
        }else if(transform.k %in% c("sd","logsd","var","logvar")){
            transform.sigma <- switch(transform.k,
                                      "sd" = "none",
                                      "logsd" = "log",
                                      "var" = "square",
                                      "logvar" = "logsquare")
        }else{
            transform.sigma <- x.transform.sigma
        }
        attr(transform.sigma,"arg") <- NULL 
    }else{
        transform.sigma.save <- transform.sigma
        attr(transform.sigma.save,"arg") <- NULL
        transform.sigma <- match.arg(transform.sigma, c("none","one","log","square","logsquare"))
        attr(transform.sigma,"arg") <- transform.sigma.save
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

    ## ** reparametrisation
    if(!is.null(p)){
        param.name <- table.param$name
        if(any(duplicated(names(p)))){
            stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
        }
        if(any(param.name %in% names(p) == FALSE)){
            stop("Incorrect argument \'p\': missing parameter(s) \"",paste(param.name[param.name %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
        }

        if((!is.null(p.transform.sigma) && p.transform.sigma!="none") || (!is.null(p.transform.k) && p.transform.k!="none") || (!is.null(p.transform.rho) && p.transform.rho!="none")){
            indexVcov <- which(table.param$type %in% c("sigma","k","rho"))
            paramVcov.name <- table.param$name[indexVcov]

            p <- c(p[table.param$name[table.param$type=="mu"]],
                   .reparametrize(p = p[paramVcov.name], type = table.param$type[indexVcov], level = table.param$level[indexVcov], 
                                  sigma = table.param$sigma[indexVcov], k.x = table.param$sigma[indexVcov], k.y = table.param$sigma[indexVcov],
                                  Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE, 
                                  transform.sigma = p.transform.sigma,
                                  transform.k = p.transform.k,
                                  transform.rho = p.transform.rho,
                                  transform.names = FALSE)$p)
        }else{
            p <- p[param.name]
        }        
    }
    
    ## ** export
    return(list(p = p,
                transform.sigma = transform.sigma,
                transform.k = transform.k,
                transform.rho = transform.rho,
                test.notransform = test.notransform))
}
##----------------------------------------------------------------------
### reparametrize.R ends here
