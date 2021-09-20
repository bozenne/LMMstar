### test-reparametrize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 25 2021 (11:54) 
## Version: 
## Last-Updated: sep 20 2021 (16:53) 
##           By: Brice Ozenne
##     Update #: 30
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

if(FALSE){
    library(testthat)
    library(numDeriv)

    library(LMMstar)
}

context("Check reparametrize")

## * no transformation
p <- c("sigma" = 1.5, "k" = 2)
type <- c("sigma", "k")
strata <- c(1, 1)
time.level <- c(1, 2)

GS.none <- reparametrize(p = p, type = type, strata = strata, time.level = time.level,
                    FUN = function(p, type, strata, time.levels, inverse){
                        p
                    })

test.none <- reparametrize(p = p, type = type, strata = strata, time.level = time.level,
                      transform.sigma = "none",
                      transform.k = "none",
                      transform.rho = "none")

expect_equal(GS.none, test.none, tol = 1e-6)

## * transformation involving a single parameter 
p <- c("sigma" = 1.5, "k" = 2, "rho" =  0.5)
type <- c("sigma", "k", "rho")
strata <- c(1, 1, 1)
time.level <- c(1, 2)

## ** log, square, atanh
GS.sp1 <- reparametrize(p = p, type = type, strata = strata, time.level = time.level,
                        FUN = function(p, type, strata, time.levels, inverse){
                          
                          if(inverse){
                              p[type=="sigma"] <- exp(p[type=="sigma"])
                              p[type=="k"] <- sqrt(p[type=="k"])
                              p[type=="rho"] <- tanh(p[type=="rho"])
                          }else{
                              p[type=="sigma"] <- log(p[type=="sigma"])
                              p[type=="k"] <- p[type=="k"]^2
                              p[type=="rho"] <- atanh(p[type=="rho"])
                          }

                            return(p)
                        }, transform.names = FALSE)
test.sp1 <- reparametrize(p = p, type = type, strata = strata, time.level = time.level,                          
                          transform.sigma = "log",
                          transform.k = "square",
                          transform.rho = "atanh",
                          transform.names = FALSE)
##     sigma         k       rho 
## 0.4054651 0.6931472 0.5493061 
## attr(,"Jacobian")
##       sigma k  rho
## sigma   1.5 0 0.00
## k       0.0 2 0.00
## rho     0.0 0 0.75

expect_equal(GS.sp1, test.sp1, tol = 1e-5)

testB.sp1 <- .reparametrize(p = setNames(as.double(test.sp1),names(p)), type = type, strata = strata, time.level = time.level,
               Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
               transform.sigma = "log",
               transform.k = "square",
               transform.rho = "atanh",
               transform.names = FALSE)
expect_equal(testB.sp1$p, p)

## ** logsquare
GS.sp2 <- reparametrize(p = p, type = type, strata = strata, time.level = time.level,
                        FUN = function(p, type, strata, time.levels, inverse){
                          
                          if(inverse){
                              p[type=="sigma"] <- exp(p[type=="sigma"]/2)
                              p[type=="k"] <- exp(p[type=="k"]/2)
                              p[type=="rho"] <- tanh(p[type=="rho"])
                          }else{
                              p[type=="sigma"] <- log(p[type=="sigma"]^2)
                              p[type=="k"] <- log(p[type=="k"]^2)
                              p[type=="rho"] <- atanh(p[type=="rho"])
                          }

                          return(p)
                      }, transform.names = FALSE)


test.sp2 <- reparametrize(p = p, type = type, strata = strata, time.level = time.level,
                      transform.sigma = "logsquare",
                      transform.k = "logsquare",
                      transform.rho = "atanh",
                      transform.names = FALSE)

expect_equal(GS.sp2, test.sp2, tol = 1e-5)

testB.sp2 <- .reparametrize(p = setNames(as.double(test.sp2),names(p)), type = type, strata = strata, time.level = time.level,
               Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
               transform.sigma = "logsquare",
               transform.k = "logsquare",
               transform.rho = "atanh",
               transform.names = FALSE)
expect_equal(testB.sp2$p, p)


## * transformation involving multiple parameters
## ** sd
GS.mp1 <- reparametrize(p = p, type = type, strata = strata, time.level = time.level,
                        FUN = function(p, type, strata, time.levels, inverse){
                            sigma <- p[type=="sigma"]
                            k <- p[type=="k"]
                            rho <- p[type=="rho"]
                            if(inverse){
                                return(setNames(c(sigma,k/sigma,rho),names(p)))
                            }else{
                                return(setNames(c(sigma*c(1,k),rho), names(p)))
                            }
                        }, transform.names = FALSE)
test.mp1 <- reparametrize(p = p, type = type, strata = strata, time.levels = time.level,
                          transform.sigma = "none",
                          transform.k = "sd",
                          transform.rho = "none",
                          transform.names = FALSE)

expect_equal(test.mp1, GS.mp1, tol = 1e-5)

testB.mp1 <- .reparametrize(p = setNames(as.double(test.mp1),names(p)), type = type, strata = strata, time.level = time.level,
               Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
               transform.sigma = "none",
               transform.k = "sd",
               transform.rho = "none",
               transform.names = FALSE)
expect_equal(testB.mp1$p, p)

## ** logsd
GS.mp2 <- reparametrize(p = p, type = type, strata = strata, time.level = time.level,
                        FUN = function(p, type, strata, time.levels, inverse){
                            sigma <- p[type=="sigma"]
                            k <- p[type=="k"]
                            rho <- p[type=="rho"]
                            if(inverse){
                                return(setNames(c(exp(sigma),exp(k-sigma),rho),names(p)))
                            }else{
                                return(setNames(c(log(sigma*c(1,k)),rho), names(p)))
                            }
                        }, transform.names = FALSE)

test.mp2 <- reparametrize(p = p, type = type, strata = strata, time.levels = time.level,
                          transform.sigma = "none",
                          transform.k = "logsd",
                          transform.rho = "none",
                          transform.names = FALSE)

expect_equal(test.mp2, GS.mp2, tol = 1e-5)

testB.mp2 <- .reparametrize(p = setNames(as.double(test.mp2),names(p)), type = type, strata = strata, time.level = time.level,
                       Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                       transform.sigma = "none",
                       transform.k = "logsd",
                       transform.rho = "none",
                       transform.names = FALSE)
expect_equal(testB.mp2$p, p)

## ** var / cov
GS.mp3 <- reparametrize(p = p, type = type, strata = strata, time.level = time.level,
                        FUN = function(p, type, strata, time.levels, inverse){
                            sigma <- p[type=="sigma"]
                            k <- p[type=="k"]
                            rho <- p[type=="rho"]
                            if(inverse){
                                return(setNames(c(sqrt(sigma),sqrt(k/sigma),rho),names(p)))
                            }else{
                                return(setNames(c(sigma^2*c(1,k^2),rho), names(p)))
                            }
                        }, transform.names = FALSE)

test.mp3 <- reparametrize(p = p, type = type, strata = strata, time.levels = time.level,
                          transform.sigma = "none",
                          transform.k = "var",
                          transform.rho = "none", transform.names = FALSE)

expect_equal(test.mp1, GS.mp1, tol = 1e-5)

testB.mp3 <- .reparametrize(p = setNames(as.double(test.mp3),names(p)), type = type, strata = strata, time.level = time.level,
               Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
               transform.sigma = "none",
               transform.k = "var",
               transform.rho = "none",
               transform.names = FALSE)
expect_equal(testB.mp3$p, p)

## ** logvar / logcov
GS.mp4 <- reparametrize(p = p, type = type, strata = strata, time.level = time.level,
                        FUN = function(p, type, strata, time.levels, inverse){
                            sigma <- p[type=="sigma"]
                            k <- p[type=="k"]
                            rho <- p[type=="rho"]
                            if(inverse){
                                return(setNames(c(exp(sigma/2),exp(k/2-sigma/2),rho),names(p)))
                            }else{
                                return(setNames(c(log(sigma^2*c(1,k^2)),rho), names(p)))
                            }
                        }, transform.names = FALSE)

test.mp4 <- reparametrize(p = p, type = type, strata = strata, time.levels = time.level,
                          transform.sigma = "none",
                          transform.k = "logvar",
                          transform.rho = "none", transform.names = FALSE)

expect_equal(test.mp4, GS.mp4, tol = 1e-5)

testB.mp4 <- .reparametrize(p = setNames(as.double(test.mp4),names(p)), type = type, strata = strata, time.level = time.level,
                            Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                            transform.sigma = "none",
                            transform.k = "logvar",
                            transform.rho = "none",
                            transform.names = FALSE)
expect_equal(testB.mp4$p, p)



## ** cov (CS)
p.bis <- c("sigma" = 1.5, "rho" =  0.5)
type.bis <- c("sigma", "rho")
strata.bis <- c(1, 1)

GS.mp5 <- reparametrize(p = p.bis, type = type.bis, strata = strata.bis, time.level = time.level,
                        FUN = function(p, type, strata, time.levels, inverse){
                            sigma <- p[type=="sigma"]
                            rho <- p[type=="rho"]
                            if(inverse){
                                return(setNames(c(sqrt(sigma),rho/sigma),names(p)))
                            }else{
                                return(setNames(c(sigma^2,sigma^2*rho), names(p)))
                            }
                        }, transform.names = FALSE)

test.mp5 <- reparametrize(p = p.bis, type = type.bis, strata = strata.bis,
                          transform.sigma = "none",
                          transform.k = "none",
                          transform.rho = "cov",
                          transform.names = FALSE)

expect_equal(test.mp5, GS.mp5, tol = 1e-5)

testB.mp5 <- .reparametrize(p = setNames(as.double(test.mp5),names(p.bis)), type = type.bis, strata = strata.bis, time.level = time.level,
                            Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                            transform.sigma = "none",
                            transform.k = "none",
                            transform.rho = "cov",
                            transform.names = FALSE)
expect_equal(testB.mp5$p, p.bis)

##----------------------------------------------------------------------
### test-reparametrize.R ends here
