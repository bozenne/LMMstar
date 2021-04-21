### test-lmm-examples.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (10:13) 
## Version: 
## Last-Updated: Apr 21 2021 (19:56) 
##           By: Brice Ozenne
##     Update #: 43
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

context("Check lmm on a simple example")


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
expect_equal(unname(coef(e.lmm, transform = 2)), unname(coef(e.lava)), tol = 1e-6)

## *** logLikelihood
expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)

## no transformation
newp <- coef(e.lmm, transform = FALSE)+1
newp.lava <- coef(e.lava) + 1 ; newp.lava["Y~~Y"] <- (sqrt(newp.lava["Y~~Y"]-1)+1)^2
expect_equal(logLik(e.lmm, p = newp), as.double(logLik(e.lava, p = newp.lava)), tol = 1e-6)

## *** score
expect_true(all(abs(score(e.lmm)) < 1e-6))

## no transformation
newp <- coef(e.lmm, transform = FALSE)+1
GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- score(e.lmm, p = newp, transform = FALSE)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p)}, x = newp.log)
test <- score(e.lmm, p = newp, transform = TRUE)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2
GS <- jacobian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(e.lmm, p = p)}, x = newp.2)
GS0 <- score(e.lava, p = newp.2)
test <- score(e.lmm, p = newp, transform = 2)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)
expect_equal(as.double(GS0), as.double(GS), tol = 1e-6)

## *** variance-covariance
test <- information(e.lmm, p = coef(e.lmm, transform = FALSE), transform = FALSE) ## using sigma
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, transform = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(e.lmm, p = coef(e.lmm, transform = FALSE), transform = TRUE) ## using log(sigma)
GS <- -hessian(func = function(p){p["sigma"]<-exp(p["sigma"]);logLik(e.lmm, p = p)}, x = coef(e.lmm, transform = TRUE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(e.lmm, p = coef(e.lmm, transform = FALSE), transform = 2) ## using sigma^2
GS <- -hessian(func = function(p){p["sigma"]<-sqrt(p["sigma"]);logLik(e.lmm, p = p)}, x = coef(e.lmm, transform = 2))
GS0 <- -hessian(func = function(p){logLik(e.lava, p = p)}, x = coef(e.lava))
expect_equal(as.double(test["sigma","sigma"]), as.double(nobs(e.lmm)[1]/(2*coef(e.lmm, effect = "variance", transform = 2)^2)), tol = 1e-6)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)
expect_equal(as.double(GS0), as.double(GS), tol = 1e-6)

test <- vcov(e.lmm, p = coef(e.lmm, transform = FALSE), transform = FALSE) ## using sigma
GS <- solve(-hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, transform = FALSE)))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## no transformation 
newp <- coef(e.lmm,transform = FALSE)+1
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
GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- score(e.lmm, p = newp, transform = FALSE)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## transformation
newp.log <- newp; newp.log["sigma"] <- log(newp.log["sigma"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p)}, x = newp.log)
test <- score(e.lmm, p = newp, transform = TRUE)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## *** variance-covariance
test1 <- information(e.lmm, p = coef(e.lmm), transform = FALSE, type = "observed") ## using sigma
test2 <- information(e.lmm, p = coef(e.lmm), transform = FALSE, type = "expected") ## using sigma
## difference: in one case (information)  : n / sigma^2
##             in the other case (hessian): n / sigma^2 - 3 * sum(residuals^2)/sigma^4
## which are equal when sum(residuals^2)/n = sigma^2
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, transform = FALSE))
expect_equal(as.double(test1), as.double(GS), tol = 1e-6)

test1 <- information(e.lmm, p = coef(e.lmm), transform = TRUE, type = "observed") ## using log(sigma)
test2 <- information(e.lmm, p = coef(e.lmm), transform = TRUE, type = "expected") ## using log(sigma)
GS <- -hessian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p)}, x = coef(e.lmm, transform = TRUE))
expect_equal(as.double(test1), as.double(GS), tol = 1e-6)

test1 <- information(e.lmm, p = coef(e.lmm), transform = 2, type = "observed") ## using sigma^2
test2 <- information(e.lmm, p = coef(e.lmm), transform = 2, type = "expected") ## using sigma^2
GS <- -hessian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(e.lmm, p = p)}, x = coef(e.lmm, transform = 2))
expect_equal(as.double(test1), as.double(GS), tol = 1e-6)
expect_equal(as.double(test1["sigma","sigma"]), as.double((nobs(e.lmm)[1]-length(coef(e.lmm, effect = "mean")))/(2*coef(e.lmm, effect = "variance", transform = 2)^2)), tol = 1e-6)

## no transformation 
newp <- coef(e.lmm, transform = FALSE)+1
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- information(e.lmm, p = newp, transform = FALSE)
expect_equal(as.double(test[1:4,1:4]), as.double(GS[1:4,1:4]), tol = 1e-6) ## does not match as some terms do not cancel

## ** multiple variance parameters (ML)
e.lmm <- lmm(Y ~ X1 + X2 + X3, variance = ~Gender|id, structure = "UN", data = d, debug = 2,
             method = "ML")
e.gls <- gls(Y ~ X1 + X2 + X3, data = dd, weights = varIdent(form=~1|Gender), method = "ML")

## *** coef
expect_equal(coef(e.lmm, effects = "mean"), coef(e.gls), tol = 1e-6)
## coef(e.lmm, transform = 2)
## sigma(e.gls)^2

## *** logLikelihood
expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)

## *** score
expect_true(all(abs(score(e.lmm, tranform = FALSE)) < 1e-4))

## no transformation
newp <- coef(e.lmm, transform = FALSE)+1
GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- score(e.lmm, p = newp, transform = FALSE)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"]); newp.log["k.0"] <- log(newp["k.0"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); p["k.0"] <- exp(p["k.0"]); logLik(e.lmm, p = p)}, x = newp.log)
test <- score(e.lmm, p = newp, transform = TRUE)
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

## *** variance-covariance
test <- information(e.lmm, p = coef(e.lmm, transform = FALSE), transform = FALSE) ## using sigma
GS <- -hessian(func = function(p){logLik(e.lmm, p = p, transform = FALSE)}, x = coef(e.lmm, transform = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(e.lmm, p = coef(e.lmm, transform = FALSE), transform = TRUE) ## using log(sigma)
GS <- -hessian(func = function(p){p["sigma"]<-exp(p["sigma"]);logLik(e.lmm, p = p, transform = FALSE)}, x = coef(e.lmm, transform = TRUE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

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
