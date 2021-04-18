### test-lmm-examples.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (10:13) 
## Version: 
## Last-Updated: Apr 17 2021 (00:16) 
##           By: Brice Ozenne
##     Update #: 28
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
    library(lava)
}

context("Check lmm on a simple example")


## * Linear regression
## ** simulate data
n <- 5e1
p <- 3
X.name <- paste0("X",1:p)
link.lvm <- paste0("Y~",X.name)
formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))

m <- lvm(formula.lvm)
distribution(m,~Id) <- Sequence.lvm(0)
set.seed(10)
d <- lava::sim(m,n)
d$id <- paste0("id",1:NROW(d))
d$time <- "t1"


## ** single variance parameter (ML)
e.lmm <- lmm(Y ~ X1 + X2 + X3, variance = ~time|id, structure = "CS", data = d, debug = 2,
             method = "ML")
e.gls <- gls(Y ~ X1 + X2 + X3, data = d, method = "ML")
e.lava <- estimate(lvm(Y~X1+X2+X3),data = d)

## *** coef
expect_equal(coef(e.lmm, effects = "mean"), coef(e.gls), tol = 1e-6)
expect_equal(unname(c(coef(e.lmm, effects = "mean"), coef(e.lmm, effects = "var")^2)), unname(coef(e.lava)), tol = 1e-6)

## *** logLikelihood
expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)

## no transformation
newp <- coef(e.lmm)+1
newp.lava <- coef(e.lava) + 1 ; newp.lava["Y~~Y"] <- (sqrt(newp.lava["Y~~Y"]-1)+1)^2
expect_equal(logLik(e.lmm, p = newp), as.double(logLik(e.lava, p = newp.lava)), tol = 1e-6)

## *** score
expect_true(all(abs(score(e.lmm)) < 1e-6))

## no transformation
newp <- coef(e.lmm)+1
GS <- jacobian(func = function(p){logLik(e.lmm, p = p, transform = FALSE)}, x = newp)
test <- score(e.lmm, p = newp, transform = FALSE)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)
## transformation
newp.log <- newp; newp.log["sigma"] <- log(newp.log["sigma"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p, transform = TRUE)}, x = newp.log)
test <- score(e.lmm, p = newp, transform = TRUE)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## *** variance-covariance
test <- information(e.lmm, p = coef(e.lmm), transform = FALSE) ## using sigma
GS <- -hessian(func = function(p){logLik(e.lmm, p = p, transform = FALSE)}, x = coef(e.lmm))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- vcov(e.lmm, p = coef(e.lmm), transform = FALSE) ## using sigma
GS <- solve(-hessian(func = function(p){logLik(e.lmm, p = p, transform = FALSE)}, x = coef(e.lmm)))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## no transformation 
newp <- coef(e.lmm)+1
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- information(e.lmm, p = newp, transform = FALSE)
expect_equal(as.double(test[1:4,1:4]), as.double(GS[1:4,1:4]), tol = 1e-6) ## does not match as some terms do not cancel

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
GS <- jacobian(func = function(p){logLik(e.lmm, p = p, transform = FALSE)}, x = newp)
test <- score(e.lmm, p = newp, transform = FALSE)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)
## transformation
newp.log <- newp; newp.log["sigma"] <- log(newp.log["sigma"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p, transform = TRUE)}, x = newp.log)
test <- score(e.lmm, p = newp, transform = TRUE)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## *** variance-covariance
test1 <- information(e.lmm, p = coef(e.lmm), transform = FALSE, type = "observed") ## using sigma
test2 <- information(e.lmm, p = coef(e.lmm), transform = FALSE, type = "expected") ## using sigma
GS <- -hessian(func = function(p){logLik(e.lmm, p = p, transform = FALSE)}, x = coef(e.lmm))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## difference: in one case (information)  : n / sigma^2
##             in the other case (hessian): n / sigma^2 - 3 * sum(residuals^2)/sigma^4
## which are equal when sum(residuals^2)/n = sigma^2


crossprod(X)/coef(e.lmm)["sigma"]^2

(nobs(e.lmm)["obs"]-length(coef(e.lmm, effects = "mean"))) / (2*coef(e.lmm)["sigma"]^4)

test <- vcov(e.lmm, p = coef(e.lmm), transform = FALSE) ## using sigma
test["sigma","sigma"] - coef(e.lmm)["sigma"]^6/(nobs(e.lmm)["obs"]-length(coef(e.lmm, effects = "mean")))

## [1] 145.568
## [1] 133.923

## [1] 145.568 -133.923

n.obs <- nobs(e.lmm)["obs"]
eps <- residuals(e.gls)
oo <- function(pp){    
    - 0.5 * n.obs * log(2*pi) - 0.5 * n.obs * log(pp^2) - 0.5 * sum(eps^2)/pp^2 -0.5 * log(det(t(X) %*% X * 1/pp^2)) + 4 * log(2*pi)/2
}
ss <- function(pp){
    - n.obs / pp + sum(eps^2)/pp^3 + 0.5 * tr(solve(t(X) %*% X * 1/pp^2) %*% (t(X) %*% X * 1/pp^2 * 2 * pp * 1/pp^2))
}
hh <- function(pp){
    term0 <- n.obs / pp^2 - 3 * sum(eps^2)/pp^4
    term1 <- (solve(t(X) %*% X * 1/pp^2) %*% (t(X) %*% X * 1/pp^2 * 2 * pp * 1/pp^2))^2
    term2 <- - 2 * solve(t(X) %*% X * 1/pp^2) %*% (t(X) %*% X * 1/pp^2 * 2 * pp * 1/pp^2 * 2 * pp * 1/pp^2)
    term3 <- solve(t(X) %*% X * 1/pp^2) %*% (t(X) %*% X * 1/pp^2 * 2 * 1/pp^2)
    print(c(term0, tr(term1),tr(term2),tr(term3)))
    return(term0 + 0.5 * tr(term1 + term2 + term3)  )
}

## vs
2 * n.obs / pp^2
n.obs/pp^2 * (1 - 3 * sum(eps^2)/n.obs/pp^2)
n.obs/pp^2 * (1 - 3 * (sum(eps^2)/n.obs/pp^2)

logLik(e.lmm, p = coef(e.lmm))
oo(coef(e.lmm, effects = "var"))
jacobian(oo, coef(e.lmm, effects = "var"))
ss(coef(e.lmm, effects = "var"))
hessian(oo, coef(e.lmm, effects = "var"))
xx <- hh(coef(e.lmm, effects = "var"))

logLik(e.lmm, p = coef(e.lmm))
score(e.lmm, p = coef(e.lmm))
information(e.lmm, p = coef(e.lmm))

test <- vcov(e.lmm, p = coef(e.lmm), transform = FALSE) ## using sigma
GS <- solve(-hessian(func = function(p){logLik(e.lmm, p = p, transform = FALSE)}, x = coef(e.lmm)))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

##             (Intercept)        X1         X2         X3   sigma
## (Intercept)    72.78400  4.963690   1.749840 -15.573352   0.000
## X1              4.96369 68.365383  -3.441063  -3.426039   0.000
## X2              1.74984 -3.441063  67.444721 -23.786061   0.000
## X3            -15.57335 -3.426039 -23.786061  69.314956   0.000
## sigma           0.00000  0.000000   0.000000   0.000000 145.568
##         [,1]   [,2]    [,3]    [,4]    [,5]
## [1,]  72.784  4.964   1.750 -15.573   0.000
## [2,]   4.964 68.365  -3.441  -3.426   0.000
## [3,]   1.750 -3.441  67.445 -23.786   0.000
## [4,] -15.573 -3.426 -23.786  69.315   0.000
## [5,]   0.000  0.000   0.000   0.000 133.923

## no transformation 
newp <- coef(e.lmm)+1
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- information(e.lmm, p = newp, transform = FALSE)
expect_equal(as.double(test[1:4,1:4]), as.double(GS[1:4,1:4]), tol = 1e-6) ## does not match as some terms do not cancel

## ** multiple variance parameters

## * No missing values
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
