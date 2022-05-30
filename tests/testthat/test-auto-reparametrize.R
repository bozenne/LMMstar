### test-auto-reparametrize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 25 2021 (11:54) 
## Version: 
## Last-Updated: maj 30 2022 (09:25) 
##           By: Brice Ozenne
##     Update #: 46
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
sigma <- c(NA, "sigma")
k.x <- k.y <- c(NA, NA)

test_that("no transformation", {
GS.none <- reparametrize(p = p, type = type,  
                         FUN = function(p, ...){
                             p
                         }, transform.names = FALSE)

test.none <- reparametrize(p = p, type = type, 
                           transform.sigma = "none",
                           transform.k = "none",
                           transform.rho = "none",
                           transform.names = FALSE)

expect_equal(GS.none, test.none, tol = 1e-6)
})

## * transformation involving a single parameter 
p <- c("sigma" = 1.5, "k" = 2, "rho" =  0.5)
type <- c("sigma", "k", "rho")
sigma <- c(NA, "sigma", "sigma")
k.x <- c(NA, NA, NA)
k.y <- c(NA, NA, "k")
level <- c("A","B","(A,B)")

## ** log, square, atanh
test_that("log, square, atanh transformation", {
GS.sp1 <- reparametrize(p = p, type = type, 
                        FUN = function(p, type, inverse, ...){
                         
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

test.sp1 <- reparametrize(p = p, type = type, 
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

testB.sp1 <- .reparametrize(p = setNames(as.double(test.sp1),names(p)), type = type, 
                            Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                            transform.sigma = "log",
                            transform.k = "square",
                            transform.rho = "atanh",
                            transform.names = FALSE)
expect_equal(testB.sp1$p, p, tol = 1e-5)
})
## ** logsquare, logsquare, none
test_that("logsquare, logsquare, none transformation", {
GS.sp2 <- reparametrize(p = p, type = type, 
                        FUN = function(p, type, inverse, ...){
                          
                          if(inverse){
                              p[type=="sigma"] <- exp(p[type=="sigma"]/2)
                              p[type=="k"] <- exp(p[type=="k"]/2)
                          }else{
                              p[type=="sigma"] <- log(p[type=="sigma"]^2)
                              p[type=="k"] <- log(p[type=="k"]^2)
                          }

                          return(p)
                      }, transform.names = FALSE)


test.sp2 <- reparametrize(p = p, type = type, 
                          transform.sigma = "logsquare",
                          transform.k = "logsquare",
                          transform.rho = "none",
                          transform.names = FALSE)

expect_equal(GS.sp2, test.sp2, tol = 1e-5)

testB.sp2 <- .reparametrize(p = setNames(as.double(test.sp2),names(p)), type = type, 
                            Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                            transform.sigma = "logsquare",
                            transform.k = "logsquare",
                            transform.rho = "none",
                            transform.names = FALSE)
expect_equal(testB.sp2$p, p, tol = 1e-5)
})

## * transformation involving multiple parameters
## ** sd
test_that("sd transformation", {
GS.mp1 <- reparametrize(p = p, type = type, 
                        FUN = function(p, type, inverse, ...){
                            sigma <- p[type=="sigma"]
                            k <- p[type=="k"]
                            rho <- p[type=="rho"]
                            if(inverse){
                                return(setNames(c(sigma,k/sigma,rho),names(p)))
                            }else{
                                return(setNames(c(sigma*c(1,k),rho), names(p)))
                            }
                        }, transform.names = FALSE)

test.mp1 <- reparametrize(p = p, type = type, sigma = sigma, 
                          transform.k = "sd",
                          transform.rho = "none",
                          transform.names = FALSE)

expect_equal(test.mp1, GS.mp1, tol = 1e-5)

testB.mp1 <- .reparametrize(p = setNames(as.double(test.mp1),names(p)), type = type, sigma = sigma,
               Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
               transform.sigma = "none",
               transform.k = "sd",
               transform.rho = "none",
               transform.names = FALSE)
expect_equal(testB.mp1$p, p, tol = 1e-5)
})

## ** logsd
test_that("logsd transformation", {
GS.mp2 <- reparametrize(p = p, type = type, 
                        FUN = function(p, type, inverse, ...){
                            sigma <- p[type=="sigma"]
                            k <- p[type=="k"]
                            rho <- p[type=="rho"]
                            if(inverse){
                                return(setNames(c(exp(sigma),exp(k-sigma),rho),names(p)))
                            }else{
                                return(setNames(c(log(sigma*c(1,k)),rho), names(p)))
                            }
                        }, transform.names = FALSE)

test.mp2 <- reparametrize(p = p, type = type, sigma = sigma, 
                          transform.k = "logsd",
                          transform.rho = "none",
                          transform.names = FALSE)
expect_equal(test.mp2, GS.mp2, tol = 1e-5)

testB.mp2 <- .reparametrize(p = setNames(as.double(test.mp2),names(p)), type = type, sigma = sigma, 
                            Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                            transform.sigma = "none",
                            transform.k = "logsd",
                            transform.rho = "none",
                            transform.names = FALSE)
expect_equal(testB.mp2$p, p, tol = 1e-5)
})

## ** var
test_that("var transformation", {
GS.mp3 <- reparametrize(p = p, type = type,
                        FUN = function(p, type, inverse, ...){
                            sigma <- p[type=="sigma"]
                            k <- p[type=="k"]
                            rho <- p[type=="rho"]
                            if(inverse){
                                return(setNames(c(sqrt(sigma),sqrt(k/sigma),rho),names(p)))
                            }else{
                                return(setNames(c(sigma^2*c(1,k^2),rho), names(p)))
                            }
                        }, transform.names = FALSE)

test.mp3 <- reparametrize(p = p, type = type, sigma = sigma, 
                          transform.k = "var",
                          transform.rho = "none", transform.names = FALSE)

expect_equal(test.mp3, GS.mp3, tol = 1e-5)

testB.mp3 <- .reparametrize(p = setNames(as.double(test.mp3),names(p)), type = type, sigma = sigma,
                            Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                            transform.sigma = "none",
                            transform.k = "var",
                            transform.rho = "none",
                            transform.names = FALSE)
expect_equal(testB.mp3$p, p, tol = 1e-5)
})

## ** logvar
test_that("logvar transformation", {
GS.mp4 <- reparametrize(p = p, type = type, 
                        FUN = function(p, type, inverse, ...){
                            sigma <- p[type=="sigma"]
                            k <- p[type=="k"]
                            rho <- p[type=="rho"]
                            if(inverse){
                                return(setNames(c(exp(sigma/2),exp(k/2-sigma/2),rho),names(p)))
                            }else{
                                return(setNames(c(log(sigma^2*c(1,k^2)),rho), names(p)))
                            }
                        }, transform.names = FALSE)

test.mp4 <- reparametrize(p = p, type = type, sigma = sigma,
                          transform.k = "logvar",
                          transform.rho = "none", transform.names = FALSE)

expect_equal(test.mp4, GS.mp4, tol = 1e-5)

testB.mp4 <- .reparametrize(p = setNames(as.double(test.mp4),names(p)), type = type, sigma = sigma, 
                            Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                            transform.sigma = "none",
                            transform.k = "logvar",
                            transform.rho = "none",
                            transform.names = FALSE)
expect_equal(testB.mp4$p, p, tol = 1e-5)
})

## ** cov (CS)
test_that("cov transformation (CS)", {
GS.mp5 <- reparametrize(p = p[c(1,3)], type = type[c(1,3)], 
                        FUN = function(p, type, inverse, ...){
                            sigma <- p[type=="sigma"]
                            rho <- p[type=="rho"]
                            if(inverse){
                                return(setNames(c(sqrt(sigma),rho/sigma),names(p)))
                            }else{
                                return(setNames(c(sigma^2,sigma^2*rho), names(p)))
                            }
                        }, transform.names = FALSE)

test.mp5 <- reparametrize(p = p[c(1,3)], type = type[c(1,3)], sigma = sigma[c(1,3)], k.x = k.x[c(1,3)], k.y = k.x[c(1,3)],
                          transform.rho = "cov",
                          transform.names = FALSE)

expect_equal(test.mp5, GS.mp5, tol = 1e-5)

testB.mp5 <- .reparametrize(p = setNames(as.double(test.mp5),names(p[c(1,3)])), type[c(1,3)], sigma = sigma[c(1,3)], k.x = k.x[c(1,3)], k.y = k.x[c(1,3)],
                            Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                            transform.sigma = "none",
                            transform.k = "none",
                            transform.rho = "cov",
                            transform.names = FALSE)
expect_equal(testB.mp5$p, p[c(1,3)], tol = 1e-5)
})

## ** cov (UN)
test_that("cov transformation (UN)", {
GS.mp6 <- reparametrize(p = p, type = type, 
                        FUN = function(p, type, inverse, ...){
                            sigma <- p[type=="sigma"]
                            k <- p[type=="k"]
                            rho <- p[type=="rho"]
                            if(inverse){
                                return(setNames(c(sqrt(sigma),sqrt(k/sigma),rho/sqrt(sigma*k)),names(p)))
                            }else{
                                return(setNames(c(sigma^2*c(1,k^2),sigma^2*k*rho), names(p)))
                            }
                        }, transform.names = FALSE)

test.mp6 <- reparametrize(p = p, type = type, sigma = sigma, k.x = k.x, k.y = k.y,
                          transform.rho = "cov",
                          transform.names = FALSE)


expect_equal(test.mp6, GS.mp6, tol = 1e-5)

testB.mp6 <- .reparametrize(p = setNames(as.double(test.mp6),names(p)), type, sigma = sigma, k.x = k.x, k.y = k.y,
                            Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                            transform.sigma = "none",
                            transform.k = "none",
                            transform.rho = "cov",
                            transform.names = FALSE)
expect_equal(testB.mp6$p, p, tol = 1e-5)


p2 <- c("sigma" = 1.5, "k1" = 2, "k2" = 0.9, "rho01" =  0.5, "rho02" =  0.25, "rho12" =  0.6)
type2 <- c("sigma", "k", "k", "rho", "rho","rho")
sigma2 <- c(NA, "sigma", "sigma", "sigma", "sigma", "sigma")
k2.x <- c(NA, NA, NA, NA, NA, "k1")
k2.y <- c(NA, NA, NA, "k1", "k2", "k2")
level2 <- c("A","B","C","A:B","A:C","B:C")

GS.mp7 <- reparametrize(p = p2, type = type2, sigma = sigma2, k.x = k2.x, k.y = k2.y,
                        FUN = function(p, type, sigma, k.x, k.y, inverse, ...){
                            k.x <- p[k.x[type=="rho"]]
                            k.y <- p[k.y[type=="rho"]]
                            sigma <- p[type=="sigma"]
                            k <- p[type=="k"]
                            rho <- p[type=="rho"]
                            if(inverse){
                                k.x[is.na(k.x)] <- sigma
                                k.y[is.na(k.y)] <- sigma
                                return(setNames(c(sqrt(sigma),sqrt(k/sigma),rho/sqrt(k.x*k.y)),names(p)))
                            }else{
                                k.x[is.na(k.x)] <- 1
                                k.y[is.na(k.y)] <- 1
                                return(setNames(c(sigma^2*c(1,k^2),rho*sigma^2*k.x*k.y), names(p)))
                            }
                        }, transform.names = FALSE)

test.mp7 <- reparametrize(p = p2, type = type2, sigma = sigma2, k.x = k2.x, k.y = k2.y,
                          transform.rho = "cov",
                          transform.names = FALSE)

expect_equal(test.mp7, GS.mp7, tol = 1e-5)

attr(test.mp7,"Jacobian")-attr(GS.mp7,"Jacobian")
rdiff <- (attr(test.mp7,"dJacobian")-attr(GS.mp7,"dJacobian"))/abs(attr(test.mp7,"dJacobian"))
range(rdiff[!is.na(rdiff) & !is.infinite(rdiff)])
})

##----------------------------------------------------------------------
### test-auto-reparametrize.R ends here
