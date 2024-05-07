### sampleRem.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (14:23) 
## Version: 
## Last-Updated: maj  7 2024 (10:13) 
##           By: Brice Ozenne
##     Update #: 86
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * sampleRem (documentation)
##' @title Sample Longuitudinal Data
##' @description Sample longuitudinal data with covariates
##'
##' @param n [integer,>0] sample size
##' @param n.times [integer,>0] number of visits (i.e. measurements per individual).
##' @param mu [numeric vector] expected measurement value at each visit (when all covariates are fixed to 0). Must have length \code{n.times}.
##' @param sigma [numeric vector,>0] standard error of the measurements at each visit (when all covariates are fixed to 0). Must have length \code{n.times}.
##' @param lambda [numeric vector] covariance between the measurement at each visit and the individual latent variable. Must have length \code{n.times}.
##' @param beta [numeric vector of length 10] regression coefficient between the covariates and the latent variable. 
##' @param gamma [numeric matrix with n.times rows and 10 columns] regression coefficient specific to each timepoint (i.e. interaction with time). 
##' @param format [character] Return the data in the wide format (\code{"wide"}) or long format (\code{"long"}).
##' Can also be \code{"wide+"} or \code{"long+"} to export as attributes the function arguments and the latent variable model used to generate the data.
##' @param latent [logical] Should the latent variable be output?
##'
##' @details The generative model is a latent variable model where each outcome (\eqn{Y_j}) load on the latent variable (\eqn{\eta}) with a coefficient lambda:
##' \deqn{Y_j = \mu_j + \lambda_j*\eta + \sigma_j\epsilon_j}
##' The latent variable is related to the covariates (\eqn{X_1,\ldots,X_(10)}):
##' \deqn{\eta = \alpha + \beta_1 X_1 + ... + \beta_{10} X_{10} + \xi}
##' \eqn{\epsilon_j} and \eqn{\xi} are independent random variables with standard normal distribution.
##'
##' @return a data.frame
##' 
##' @keywords datagen

## * sampleRem (examples)
##' @examples
##' #### generate data in the wide format ####
##' set.seed(10)
##' dW <- sampleRem(100, n.times = 3, format = "wide")
##' head(dW)
##'
##' #### generate data in the long format ####
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' head(dL)

## * sampleRem (code)
##' @export
sampleRem <- function(n, n.times,
                      mu = 1:n.times,
                      sigma = rep(1,n.times),
                      lambda = rep(1,n.times),
                      beta = c(2,1,0,0,0,1,1,0,0,0),
                      gamma = matrix(0, nrow = n.times, ncol = 10),
                      format = "wide",
                      latent = FALSE){

    name.Y <- paste0("Y",1:n.times)

    ## ** check arguments
    if(identical(format,"wide+")){
        format <- "wide"
        exportArgs <- TRUE
    }else if(identical(format,"long+")){
        format <- "long"
        exportArgs <- TRUE
    }else{
        exportArgs <- FALSE
    }
    format <- match.arg(format, c("wide","long"))
    
    if(length(mu)!=n.times){
        stop("Argument \'mu\' must have length argument \'n.times\' \n")
    }
    if(length(sigma)!=n.times){
        stop("Argument \'sigma\' must have length argument \'n.times\' \n")
    }
    if(length(lambda)!=n.times){
        stop("Argument \'lambda\' must have length argument \'n.times\' \n")
    }
    ## ** generative model
    m <- lava::lvm()
    m <- lava::`latent<-`(m, value = ~eta)
    
    ## covariates
    m <- lava::`distribution<-`(m, "X1", value = lava::binomial.lvm(size = 1, p = 0.5))
    m <- lava::`distribution<-`(m, "X2", value = lava::binomial.lvm(size = 1, p = 0.1))
    m <- lava::`distribution<-`(m, "X3", value = lava::binomial.lvm(size = 2, p = 0.5))
    m <- lava::`distribution<-`(m, "X4", value = lava::binomial.lvm(size = 3, p = 0.5))
    m <- lava::`distribution<-`(m, "X5", value = lava::poisson.lvm(lambda = 1))

    m <- lava::`distribution<-`(m, "X6", value = lava::gaussian.lvm())
    m <- lava::`distribution<-`(m, "X7", value = lava::gaussian.lvm(sd = 3))
    m <- lava::`distribution<-`(m, "X8", value = lava::gaussian.lvm(mean = 2, sd = 2))
    m <- lava::`distribution<-`(m, "X9", value = lava::Gamma.lvm(shape = 1, rate = 1))
    m <- lava::`distribution<-`(m, "X10", value = lava::beta.lvm())

    m <- lava::`regression<-`(m, eta ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, value = beta)

    ## outcome
    for(iT in 1:n.times){ ## iT <- 1
        m <- lava::`distribution<-`(m, name.Y[iT], value = lava::gaussian.lvm(mean = mu[iT], sd = sigma[iT]))
        m <- lava::`regression<-`(m, eta ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, value = beta)
        m <- lava::`regression<-`(m, stats::as.formula(paste0(name.Y[iT],"~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10")), value = gamma[iT,])
        m <- lava::`regression<-`(m, stats::as.formula(paste0(name.Y[iT],"~eta")), value =  lambda[iT])
    }
    
    ## ** generate data
    d <- lava::sim(m, n = n, latent = latent)
    d <- cbind(id = 1:n, d)
    
    ## ** reshape
    if(format == "long"){
        d <- stats::reshape(d, direction  = "long",
                            idvar = "id",
                            varying = name.Y,
                            v.names = "Y",
                            timevar = "visit")
        d$visit <- as.factor(d$visit)
        d <- d[order(d$id,d$visit),c("id","visit","Y",paste0("X",1:10),if(latent){"eta"})]
        rownames(d) <- NULL
    }
    
    ## ** export
    if(exportArgs){
        attr(d,"call") <- match.call()
        attr(d,"mu") <- mu
        attr(d,"sigma") <- sigma
        attr(d,"lambda") <- lambda
        attr(d,"beta") <- beta
        attr(d,"gamma") <- gamma
        attr(d,"lvm") <- m
    }
    return(d)
}

######################################################################
### sampleRem.R ends here
