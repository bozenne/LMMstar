### test-lmm-examples.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (10:13) 
## Version: 
## Last-Updated: May  4 2021 (23:47) 
##           By: Brice Ozenne
##     Update #: 64
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
    library(lava)

    library(LMMstar)
}

context("Check lmm on simple examples")


## * Linear regression
## ** simulate data
n <- 5e1
p <- 3
X.name <- paste0("X",1:p)
link.lvm <- paste0("Y~",X.name)
formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))

m <- lvm(formula.lvm)
distribution(m,~Id) <- Sequence.lvm(a = 1, b = n)
distribution(m,~Gender) <- binomial.lvm()
set.seed(10)
d <- lava::sim(m,n)
d$id <- paste0("id",1:NROW(d))
d$time <- "t1"


## ** single variance parameter (ML)
e.lmm <- lmm(Y ~ X1 + X2 + X3, variance = ~time|id, structure = "CS", data = d, debug = 2,
             method = "ML")
expect_warning(lmm(Y ~ X1 + X2 + X3, variance = ~time|id, structure = "UN", data = d))
e.gls <- gls(Y ~ X1 + X2 + X3, data = d, method = "ML")
e.lava <- estimate(lvm(Y~X1+X2+X3),data = d)

## *** coef
expect_equal(coef(e.lmm, effects = "mean"), coef(e.gls), tol = 1e-6)
expect_equal(unname(coef(e.lmm, transform.sigma = "square")), unname(coef(e.lava)), tol = 1e-6)


## *** logLikelihood
expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)
## coef(e.lmm, transform.sigma = "log")
## coef(e.lmm, transform.sigma = "square")
## coef(e.lmm, transform.sigma = "logsquare")

## no transformation
newp <- coef(e.lmm)+1
newp.lava <- coef(e.lava) + 1 ; newp.lava["Y~~Y"] <- (sqrt(newp.lava["Y~~Y"]-1)+1)^2
expect_equal(logLik(e.lmm, p = newp), as.double(logLik(e.lava, p = newp.lava)), tol = 1e-6)

## *** score
expect_true(all(abs(score(e.lmm)) < 1e-6))
## score(e.lmm, transform.sigma = "log")
## score(e.lmm, transform.sigma = "square")
## score(e.lmm, transform.sigma = "logsquare")

## no transformation
newp <- coef(e.lmm, transform.sigma = "none")+1
GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- score(e.lmm, p = newp, transform.sigma = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p)}, x = newp.log)
test <- score(e.lmm, p = newp, transform.sigma = "log")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2
GS <- jacobian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(e.lmm, p = p)}, x = newp.2)
GS0 <- score(e.lava, p = newp.2)
test <- score(e.lmm, p = newp, transform.sigma = "square")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)
expect_equal(as.double(GS0), as.double(GS), tol = 1e-6)

## *** information
test <- information(e.lmm, p = coef(e.lmm, transform.sigma = "none"), transform.sigma = "none")
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, transform.sigma = "none"))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(e.lmm, p = coef(e.lmm, transform.sigma = "none"), transform.sigma = "log")
GS <- -hessian(func = function(p){p["sigma"]<-exp(p["sigma"]);logLik(e.lmm, p = p)}, x = coef(e.lmm, transform.sigma = "log", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(e.lmm, p = coef(e.lmm, transform.sigma = "none"), transform.sigma = "square") 
GS <- -hessian(func = function(p){p["sigma"]<-sqrt(p["sigma"]);logLik(e.lmm, p = p)}, x = coef(e.lmm, transform.sigma = "square", transform.names = FALSE))
GS0 <- -hessian(func = function(p){logLik(e.lava, p = p)}, x = coef(e.lava))
expect_equal(as.double(test["sigma^2","sigma^2"]), as.double(nobs(e.lmm)[1]/(2*coef(e.lmm, effect = "variance", transform.sigma = "square")^2)), tol = 1e-6)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)
expect_equal(as.double(GS0), as.double(GS), tol = 1e-6)

## no transformation 
newp <- coef(e.lmm, transform.sigma = "none")+1
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- information(e.lmm, p = newp, transform.sigma = "none")
expect_equal(as.double(test[1:4,1:4]), as.double(GS[1:4,1:4]), tol = 1e-6) ## does not match as some terms do not cancel

## *** variance-covariance
test <- vcov(e.lmm, p = coef(e.lmm, transform.sigma = "none"), transform.sigma = "none")
GS <- solve(-hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, transform.sigma = "none")))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- vcov(e.lmm, p = coef(e.lmm, transform.sigma = "none"), transform.sigma = "square", df = TRUE) 
GS <- solve(-hessian(func = function(p){p["sigma"] <- sqrt(p["sigma"]);logLik(e.lmm, p = p)}, x = coef(e.lmm, transform.sigma = "square", transform.names = FALSE)))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

GS <- c(rep(nobs(e.lmm)[1],sum(e.lmm$param$type=="mu")), nobs(e.lmm)[1]/4)
expect_equal(unname(attr(test,"df")), unname(GS))

## ** single variance parameter (REML)
e.lmm <- lmm(Y ~ X1 + X2 + X3, variance = ~time|id, structure = "CS", data = d, debug = 2,
             method = "REML")
e.gls <- gls(Y ~ X1 + X2 + X3, data = d, method = "REML")

## *** coef
expect_equal(coef(e.lmm, effects = "mean"), coef(e.gls), tol = 1e-6)

## *** logLikelihood
expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)

## *** score
expect_true(all(abs(score(e.lmm)) < 1e-6))

## no transformation
newp <- coef(e.lmm)+1
GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- score(e.lmm, p = newp, transform.sigma = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## transformation
newp.log <- newp; newp.log["sigma"] <- log(newp.log["sigma"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p)}, x = newp.log)
test <- score(e.lmm, p = newp, transform.sigma = "log")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## *** information
test1 <- information(e.lmm, p = coef(e.lmm), transform.sigma = "none", type = "observed")
test2 <- information(e.lmm, p = coef(e.lmm), transform.sigma = "none", type = "expected")
## difference: in one case (information)  : n / sigma^2
##             in the other case (hessian): n / sigma^2 - 3 * sum(residuals^2)/sigma^4
## which are equal when sum(residuals^2)/n = sigma^2
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, transform.sigma = "none"))
expect_equal(as.double(test1), as.double(GS), tol = 1e-6)

test1 <- information(e.lmm, p = coef(e.lmm), transform.sigma = "log", type = "observed")
test2 <- information(e.lmm, p = coef(e.lmm), transform.sigma = "log", type = "expected")
GS <- -hessian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p)}, x = coef(e.lmm, transform.sigma = "log", transform.names = FALSE))
expect_equal(as.double(test1), as.double(GS), tol = 1e-6)

test1 <- information(e.lmm, p = coef(e.lmm), transform.sigma = "square", type = "observed")
test2 <- information(e.lmm, p = coef(e.lmm), transform.sigma = "square", type = "expected")
GS <- -hessian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(e.lmm, p = p)}, x = coef(e.lmm, transform.sigma = "square", transform.names = FALSE))
expect_equal(as.double(test1), as.double(GS), tol = 1e-6)
expect_equal(as.double(test1["sigma^2","sigma^2"]), as.double((nobs(e.lmm)[1]-length(coef(e.lmm, effect = "mean")))/(2*coef(e.lmm, effect = "variance", transform = 2)^2)), tol = 1e-6)

## no transformation 
newp <- coef(e.lmm, transform = FALSE)+1
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- information(e.lmm, p = newp, transform = FALSE)
expect_equal(as.double(test[1:4,1:4]), as.double(GS[1:4,1:4]), tol = 1e-6) ## does not match as some terms do not cancel

## *** variance-covariance
test <- attr(vcov(e.lmm, df = TRUE, transform = 2),"df")
GS <- c(rep(nobs(e.lmm)[1]-sum(e.lmm$param$type=="mu"),sum(e.lmm$param$type=="mu")), (nobs(e.lmm)[1]-sum(e.lmm$param$type=="mu"))/4)
expect_equal(unname(test), unname(GS))

## ** multiple variance parameters (ML)
e.lmm <- lmm(Y ~ X1 + X2 + X3, variance = ~Gender|id, structure = "UN", data = d, debug = 2,
             method = "ML")
e.gls <- gls(Y ~ X1 + X2 + X3, data = d, weights = varIdent(form=~1|Gender), method = "ML")

## *** coef
expect_equal(coef(e.lmm, effects = "mean"), coef(e.gls), tol = 1e-6)
## coef(e.lmm, transform = 2)
## sigma(e.gls)^2

## *** logLikelihood
expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)
## coef(e.lmm, transform.sigma = "log", transform.k = "log")
## coef(e.lmm, transform.sigma = "square", transform.k = "square")
## coef(e.lmm, transform.sigma = "logsquare", transform.k = "logsquare")

## *** score
expect_true(all(abs(score(e.lmm, transform.sigma = "none", transform.k = "none")) < 1e-4))

## no transformation
newp <- coef(e.lmm, transform.sigma = "none", transform.k = "none")+1
GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- score(e.lmm, p = newp, transform.sigma = "none", transform.k = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"]); newp.log["k.0"] <- log(newp["k.0"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); p["k.0"] <- exp(p["k.0"]); logLik(e.lmm, p = p)}, x = newp.log)
test <- score(e.lmm, p = newp, transform.sigma = "log", transform.k = "log")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2; newp.2["k.0"] <- newp["k.0"]^2*newp["sigma"]^2
GS <- jacobian(func = function(p){p[c("sigma","k.0")] <- c(sqrt(p["sigma"]),sqrt(p["k.0"]/p["sigma"])); logLik(e.lmm, p = p, transform = FALSE)}, x = newp.2)
test <- score(e.lmm, p = newp, transform = 2)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## TO KEEP: DEBUG TRANSFORMATION
## .transformDeriv(transform = 2, sigma = newp["sigma"], k = newp["k.0"], rho = NULL, pattern.param = attr(e.lmm$design$X.var,"Upattern.param")[[2]])
## FCT_TRANS <- function(p){
##     c(sqrt(p["sigma"]), sqrt(p["k.0"]/p["sigma"]))
## }
## jacobian(FCT_TRANS, newp.2[c("sigma","k.0")])
## c(1/(2*sqrt(newp.2["sigma"])),-sqrt(newp.2["k.0"])/(2*newp.2["sigma"]^(3/2)), 1/(2*sqrt(newp.2["sigma"])*sqrt(newp.2["k.0"])))

## *** information
test0 <- information(e.lmm, p = coef(e.lmm, transform = FALSE), transform = FALSE, type.information = "expected") ## using sigma
test <- information(e.lmm, p = coef(e.lmm, transform = FALSE), transform = FALSE, type.information = "observed") ## using sigma
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, transform = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test0 <- information(e.lmm, p = coef(e.lmm, transform.sigma = "none"), transform.sigma = "log", type.information = "expected") ## using log(sigma) and log(k)
test <- information(e.lmm, p = coef(e.lmm, transform.sigma = "none"), transform.sigma = "log", type.information = "observed") ## using log(sigma) and log(k)
GS <- -hessian(func = function(p){p["sigma"]<-exp(p["sigma"]);p["k.0"]<-exp(p["k.0"]);logLik(e.lmm, p = p, transform = FALSE)}, x = coef(e.lmm, transform = TRUE, transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test0 <- information(e.lmm, p = coef(e.lmm, transform = FALSE), transform = 2, type.information = "expected") ## using log(sigma) and log(k)
test <- information(e.lmm, p = coef(e.lmm, transform = FALSE), transform = 2, type.information = "observed") ## using log(sigma) and log(k)
GS <- -hessian(func = function(p){p["sigma"]<-exp(p["sigma"]);p["k.0"]<-exp(p["k.0"]);logLik(e.lmm, p = p, transform = FALSE)}, x = coef(e.lmm, transform = TRUE))
expect_equal(as.double(test[1:4,1:4]), as.double(GS[1:4,1:4]), tol = 1e-6)
expect_equal(as.double(test[5:6,5:6]), as.double(GS[5:6,5:6]), tol = 1e-6)

## *** debug
FCT_OMEGA <- function(p, transform, vectorize){
    rho <- 0.5
    if(transform){
        out <- matrix(c(p[1],rho*sqrt(p[1]*p[2]), rho*sqrt(p[1]*p[2]),p[2]),2,2)
    }else{
        out <- p[1]^2*matrix(c(1,rho*p[2], rho*p[2],p[2]^2),2,2)
    }
    if(vectorize){
        return(as.vector(out))
    }else{
        out
    }
}
FCT_dOMEGA <- function(p, transform, vectorize){
    rho <- 0.5
    if(transform){
        out <- unname(cbind(c(1,rho*sqrt(p[2]/p[1])/2, rho*sqrt(p[2]/p[1])/2,0),
                            c(0,rho*sqrt(p[1]/p[2])/2, rho*sqrt(p[1]/p[2])/2,1)))
    }else{
        out <- unname(cbind(2*p[1]*c(1,rho*p[2], rho*p[2],p[2]^2),
                            p[1]^2*c(0,rho, rho,2*p[2])))
    }
    if(vectorize){
        return(as.vector(out))
    }else{
        out
    }
}
FCT_TRANS <- function(p){
    c(sqrt(p[1]), sqrt(p[2]/p[1]))
}
FCT_dTRANS <- function(p, vectorize){
    out <- matrix(c(1/(2*sqrt(p[1])), -sqrt(p[2])/(2*p[1]^{3/2}), 0, 1/(2*sqrt(p[2]*p[1]))),2,2) ## c(1/(2*p[1]), -p[2]/(2*p[1]^{2}), 0, p[1]/(2*p[2]))
    if(vectorize){
        return(as.vector(out))
    }else{
        out
    }
}
FCT_ddTRANS <- function(p){
    ## out <- cbind(c(-1/(4*p[1]^(3/2)), 3*sqrt(p[2])/(2*p[1]^{5/2}), 0, -1/(4*sqrt(p[2])*p[1]^(3/2))),
    ##              c(0, -1/(4*sqrt(p[2])*p[1]^{3/2}), 0, -1/(4*sqrt(p[1])*p[2]^(3/2))))
    out <- cbind(c(1/(4*p[1]^2), -p[2]/(2*p[1]^{2}), 0, p[1]/(2*p[2])),
                 c(1/(2*sqrt(p[1])), -sqrt(p[2])/(2*p[1]^{3/2}), 0, 1/(2*sqrt(p[2]*p[1]))))
    out
}
FCT_OMEGA(coef(e.lmm, effects = "variance"), transform = FALSE, vectorize=FALSE)
FCT_OMEGA(coef(e.lmm, transform.sigma = "square", effects = "variance"), transform.sigma = "square", vectorize=FALSE)

FCT_dOMEGA(coef(e.lmm, effects = "variance"), transform = FALSE, vectorize=FALSE)
GS0 <- jacobian(function(p){FCT_OMEGA(p, transform = FALSE, vectorize = TRUE)}, coef(e.lmm, effects = "variance"))
FCT_dOMEGA(coef(e.lmm, transform = 2, effects = "variance"), transform = TRUE, vectorize=FALSE)
GS1 <- jacobian(function(p){FCT_OMEGA(p, transform = TRUE, vectorize = TRUE)}, coef(e.lmm, transform = 2, effects = "variance"))
FCT_dTRANS(coef(e.lmm, transform = 2, effects = "variance"), vectorize = FALSE)
jac <- jacobian(FCT_TRANS, coef(e.lmm, transform = 2, effects = "variance"))

GS0 %*% jac - GS1
GS0[1,,drop=FALSE] %*% jac[,1,drop=FALSE] - GS1[1,1]

dGS0 <- jacobian(function(p){FCT_dOMEGA(p, transform = FALSE, vectorize = TRUE)}, coef(e.lmm, effects = "variance"))
dGS1 <- jacobian(function(p){FCT_dOMEGA(p, transform = TRUE, vectorize = TRUE)}, coef(e.lmm, transform = 2, effects = "variance"))
djac <- jacobian(function(p){FCT_dTRANS(p, vectorize = TRUE)}, coef(e.lmm, transform = 2, effects = "variance"))

dGS00 <- jacobian(function(p){FCT_dOMEGA(p, transform = FALSE, vectorize = FALSE)[1,,drop=FALSE]}, coef(e.lmm, effects = "variance"))
djac00 <- jacobian(function(p){FCT_dTRANS(p, vectorize = FALSE)[,1,drop=FALSE]}, coef(e.lmm, transform = 2, effects = "variance"))
dGS10 <- jacobian(function(p){FCT_dOMEGA(p, transform = TRUE, vectorize = FALSE)[1,1]}, coef(e.lmm, transform = 2, effects = "variance"))



t(dGS00 %*% jac[,1,drop=FALSE]) + GS0[1,,drop=FALSE] %*% FCT_ddTRANS(coef(e.lmm, transform = 0, effects = "variance")) - dGS10
2/(2*sqrt(coef(e.lmm, transform = 2, effects = "variance")[1])) - 2*coef(e.lmm, effects = "variance")[1] * 1/(2*coef(e.lmm, transform = 2, effects = "variance")[1])

2*sqrt(coef(e.lmm, transform = 2, effects = "variance")[1])
1/(2*sqrt(coef(e.lmm, transform = 2, effects = "variance")[1]))

2
-1/(2*coef(e.lmm, transform = 2, effects = "variance")[1])

x <- 2; y  <- 5; a <- sqrt(x); b <- sqrt(y/x)
jacobian(FCT_TRANS, c(x,y))
1/(2*a); -b/(2*a^2); 1/(2*a^2*b)
1/(2*sqrt(x)); -sqrt(y)/(2*x^(3/2)); 1/(2*sqrt(x)*sqrt(y))

FCT_TRANS <- function(p){
    p[c("sigma", "k.0")]  <- c(sqrt(p["sigma"]), sqrt(p["k.0"]/p["sigma"]))
    return(p)
}



scoreGS <- -jacobian(func = function(p){logLik(e.lmm, p = p, transform = FALSE)}, x = coef(e.lmm))
scoreGS.trans <- -jacobian(func = function(p){logLik(e.lmm, p = FCT_TRANS(p), transform = FALSE)}, x = coef(e.lmm, transform = 2, transform.names = FALSE))
## GS.trans <- -jacobian(func = function(p){p[c("sigma","k.0")]<- c(sqrt(p["sigma"]), sqrt(p["k.0"]/p["sigma"]));logLik(e.lmm, p = p, transform = FALSE)}, x = coef(e.lmm, transform = 2, transform.names = FALSE))

matrix(jacobian(FCT_3,point) %*% FCT_1(point), 6, 6)
matrix(FCT_1(point) %*% t(jacobian(FCT_3,point)), 6, 6)

array(jacobian(FCT_3,point) %*% FCT_1(point), c(6, 6, 6))

infoGS.trans + t(FCT_2(point)) %*% jacobian(func = function(p){FCT_1(p)}, x = point) 


infoGS <- -hessian(func = function(p){logLik(e.lmm, p = p, transform = FALSE)}, x = coef(e.lmm))
infoGS.trans <- -hessian(func = function(p){logLik(e.lmm, p = FCT_TRANS(p), transform = FALSE)}, x = coef(e.lmm, transform = 2, transform.names = FALSE))

FCT_1 <- function(p){score(e.lmm, p = FCT_TRANS(p))}
FCT_2 <- function(p){jacobian(func = FCT_TRANS, x = p)}
FCT_3 <- function(p){as.double(jacobian(func = FCT_TRANS, x = p))}
point <- coef(e.lmm, transform = 2, transform.names = FALSE)
infoGS.trans2 <- -jacobian(func = function(p){FCT_1(p) %*% FCT_2(p)}, x = point)

infoGS.trans2 <- -jacobian(func = function(p){FCT_1(p)}, x = point) %*% t(FCT_2(point))

-3.774442e+00 * 0.7909782
-3.774442e+00 * 0.7829036
[5,]    0    0    0    0  0.6272681 0.0000000
[6,]    0    0    0    0 -0.7829036 0.7909782



round(infoGS.trans - infoGS.trans2,2)

term1 <- -jacobian(func = function(p){score(e.lmm, p = FCT_TRANS(p))}, x = coef(e.lmm, transform = 2, transform.names = FALSE)) %*% jacobian(func = FCT_TRANS, x = coef(e.lmm, transform = 2, transform.names = FALSE))
term2 <- score(e.lmm, p = coef(e.lmm, transform = 0, transform.names = FALSE))


infoGS.trans3 <-  -jacobian(func = function(p){ %*% jacobian(func = FCT_TRANS, x = p)}, x = coef(e.lmm, transform = 2, transform.names = FALSE))


M.d2Trans <- jacobian(func = function(p){as.double(jacobian(func = FCT_TRANS, x = p))}, x = coef(e.lmm))
A.d2Trans <- array(NA, dim = rep(length(coef(e.lmm)),3), dimnames = list(names(coef(e.lmm)),names(coef(e.lmm)),names(coef(e.lmm))))
for(iP in 1:length(coef(e.lmm))){
    A.d2Trans[,,iP] <- M.d2Trans[,iP]
}
GS[,1] %*% A.d2Trans[1,1,]
GS-GS.trans

infoGS %*% jacobian(func = FCT_TRANS, x = coef(e.lmm)) 

test <- information(e.lmm, p = coef(e.lmm, transform = FALSE), transform = 2) ## using sigma^2
GS <- -hessian(func = function(p){p["sigma"]<-sqrt(p["sigma"]);logLik(e.lmm, p = p, transform = FALSE)}, x = coef(e.lmm, transform = 2))
GS0 <- -hessian(func = function(p){logLik(e.lava, p = p)}, x = coef(e.lava))
expect_equal(as.double(test["sigma","sigma"]), as.double(nobs(e.lmm)[1]/(2*coef(e.lmm, effect = "variance", transform = 2)^2)), tol = 1e-6)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)
expect_equal(as.double(GS0), as.double(GS), tol = 1e-6)

test <- vcov(e.lmm, p = coef(e.lmm, transform = FALSE), transform = FALSE) ## using sigma
GS <- solve(-hessian(func = function(p){logLik(e.lmm, p = p, transform = FALSE)}, x = coef(e.lmm, transform = FALSE)))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## no transformation 
newp <- coef(e.lmm,transform = FALSE)+1
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- information(e.lmm, p = newp, transform = FALSE)
expect_equal(as.double(test[1:4,1:4]), as.double(GS[1:4,1:4]), tol = 1e-6) ## does not match as some terms do not cancel

## *** variance-covariance

## * Mixed model
## ** simulate data
m <- lvm(c(Y1,Y2,Y3,Y4) ~ age + gender)
categorical(m, labels = c("male","female")) <- ~gender
transform(m, id~gender) <- function(x){1:NROW(x)}
distribution(m, ~age) <- gaussian.lvm(mean = 50, sd = 10)

set.seed(10)
dW <- lava::sim(m, 1e2)

## move to the long format
name.varying <- paste0("Y",1:4)
dL <- reshape(dW, direction  = "long",
              idvar = c("id","age","gender"),
              varying = name.varying,
              v.names = "Y",
              timevar = "visit")
rownames(dL) <- NULL
dL$visit <- factor(dL$visit,
                   levels = 1:length(name.varying),
                   labels = name.varying)
 
## ** fit lmm
eCS.lmm <- lmm(Y ~ visit + age + gender, variance = ~visit|id, structure = "CS", data = dL, debug = 2, method = "ML")
eUN.lmm <- lmm(Y ~ visit + age + gender, variance = ~visit|id, structure = "UN", data = dL, debug = 2)

eCSs.lmm <- lmm(Y ~ visit*gender + age*gender, variance = gender~visit|id, structure = "CS", data = dL, debug = 2)
eUNs.lmm <- lmm(Y ~ visit*gender + age*gender, variance = gender~visit|id, structure = "UN", data = dL, debug = 2)

## ** coef method
coef(eCS.lmm)
coef(eCS.lmm, type = "gls", strata = "1")

coef(eCSs.lmm)
coef(eCSs.lmm, type = "gls", strata = "male")

coef(eUN.lmm)
coef(eUN.lmm, type = "gls", strata = "1")

## ** formula method
formula(eCS.lmm)

## ** getVarCov method
getVarCov(eCS.lmm)
getVarCov(eCS.lmm, type = "gls")

getVarCov(eCSs.lmm)
getVarCov(eCSs.lmm, type = "gls")

getVarCov(eUN.lmm)
getVarCov(eUN.lmm, type = "gls")

## ** model.matrix
model.matrix(eCS.lmm)

## ** model.matrix
nobs(eCS.lmm)

## ** model.matrix
sd(residuals(eCS.lmm))
sd(residuals(eCS.lmm, type.residual = "pearson"))
sd(residuals(eCS.lmm, type.residual = "normalized"))

residuals(eCS.lmm, format = "long")

## ** score method
score(eCS.lmm)
score(eCS.lmm, data = dL, p = coef(eCS.lmm, effects = c("mean","variance")))
score(eCS.lmm, type = "gls", strata = "1")

## ** summary method
summary(eCS.lmm)

## ** vcov method
vcov(eCS.lmm)
vcov(eCS.lmm, data = dL, p = coef(eCS.lmm, effects = c("mean","variance")))
vcov(eCS.lmm, effects = c("mean","variance"))
vcov(eCS.lmm, type = "gls", strata = "1")




##----------------------------------------------------------------------
### test-lmm-examples.R ends here
