### test-lmm-examples.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (10:13) 
## Version: 
## Last-Updated: May 10 2021 (19:35) 
##           By: Brice Ozenne
##     Update #: 79
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
    library(multcomp)
    library(emmeans)

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
distribution(m,~Gender.num) <- binomial.lvm()
transform(m,Gender~Gender.num) <- function(x){factor(x,levels=0:1,labels=c("M","F"))}
transform(m,id~Gender.num) <- function(x){paste0("id",1:NROW(x))}
latent(m) <- ~Gender.num
set.seed(10)
d <- lava::sim(m,n, latent = FALSE)
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
newp <- coef(e.lmm, transform.sigma = "none")+1
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- information(e.lmm, p = newp, transform.sigma = "none")
expect_equal(as.double(test[1:4,1:4]), as.double(GS[1:4,1:4]), tol = 1e-6) ## does not match as some terms do not cancel

## *** variance-covariance
test <- attr(vcov(e.lmm, df = TRUE, transform.sigma = "square"),"df")
GS <- c(rep(nobs(e.lmm)[1]-sum(e.lmm$param$type=="mu"),sum(e.lmm$param$type=="mu")), (nobs(e.lmm)[1]-sum(e.lmm$param$type=="mu"))/4)
expect_equal(unname(test), unname(GS))

## ** multiple variance parameters (ML)
e.lmm <- lmm(Y ~ -1 + Gender + (X1 + X2 + X3):Gender, variance = ~Gender|id, structure = "UN", data = d, debug = 2,
             method = "ML")

e.gls <- gls(Y ~ -1 + Gender + (X1 + X2 + X3):Gender, data = d, weights = varIdent(form=~1|Gender), method = "ML")
e.lmm2 <- lmm(Y ~(X1 + X2 + X3)*Gender, variance = Gender~time|id, structure = "CS", data = d, debug = 2,
             method = "ML")

## *** coef
expect_equal(coef(e.lmm, effects = "mean"), coef(e.gls), tol = 1e-6)
## coef(e.lmm, transform = 2)
## sigma(e.gls)^2

## *** logLikelihood
expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)
## coef(e.lmm, transform.sigma = "log", transform.k = "log")
## coef(e.lmm, transform.sigma = "square", transform.k = "square")
## coef(e.lmm, transform.sigma = "logsquare", transform.k = "logsquare")
expect_equal(logLik(e.lmm2),logLik(e.lmm))

## *** score
expect_true(all(abs(score(e.lmm, transform.sigma = "none", transform.k = "none")) < 1e-4))

## no transformation
newp <- coef(e.lmm, transform.sigma = "none", transform.k = "none")+1
GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- score(e.lmm, p = newp, transform.sigma = "none", transform.k = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"]); newp.log["k.M"] <- log(newp["k.M"])
GS <- jacobian(func = function(p){p[c("sigma","k.M")] <- exp(p[c("sigma","k.M")]); logLik(e.lmm, p = p)}, x = newp.log)
test <- score(e.lmm, p = newp, transform.sigma = "log", transform.k = "log")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2; newp.2["k.M"] <- newp["k.M"]^2*newp["sigma"]^2
GS <- jacobian(func = function(p){p[c("sigma","k.M")] <- c(sqrt(p["sigma"]),sqrt(p["k.M"]/p["sigma"])); logLik(e.lmm, p = p, transform.sigma = "none")}, x = newp.2)
test <- score(e.lmm, p = newp, transform.sigma = "square", transform.k = "var")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## TO KEEP: DEBUG TRANSFORMATION
## .transformDeriv(transform = 2, sigma = newp["sigma"], k = newp["k.M"], rho = NULL, pattern.param = attr(e.lmm$design$X.var,"Upattern.param")[[2]])
## FCT_TRANS <- function(p){
##     c(sqrt(p["sigma"]), sqrt(p["k.M"]/p["sigma"]))
## }
## jacobian(FCT_TRANS, newp.2[c("sigma","k.M")])
## c(1/(2*sqrt(newp.2["sigma"])),-sqrt(newp.2["k.M"])/(2*newp.2["sigma"]^(3/2)), 1/(2*sqrt(newp.2["sigma"])*sqrt(newp.2["k.M"])))

## *** information
test0 <- information(e.lmm, p = coef(e.lmm, transform.sigma = "none", transform.k = "none"), transform.sigma = "none", transform.k = "none", type.information = "expected") ## using sigma
test <- information(e.lmm, p = coef(e.lmm, transform.sigma = "none", transform.k = "none"), transform.sigma = "none", transform.k = "none", type.information = "observed") ## using sigma
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, transform.sigma = "none", transform.k = "none"))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test0 <- information(e.lmm, p = coef(e.lmm, transform.sigma = "none", transform.k = "none"), transform.sigma = "log", transform.k = "log", type.information = "expected") ## using log(sigma) and log(k)
test <- information(e.lmm, p = coef(e.lmm, transform.sigma = "none", transform.k = "none"), transform.sigma = "log", transform.k = "log", type.information = "observed") ## using log(sigma) and log(k)
GS <- -hessian(func = function(p){p[c("sigma","k.M")] <- exp(p[c("sigma","k.M")]);logLik(e.lmm, p = p, transform.sigma = "none", transform.k = "none")}, x = coef(e.lmm, transform.sigma = "log", transform.k = "log", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test0 <- information(e.lmm, p = coef(e.lmm, transform.sigma = "none", transform.k = "none"), transform.sigma = "square", transform.k = "var", type.information = "expected") ## using sigma^2 and sigma^2 k^2
test <- information(e.lmm, p = coef(e.lmm, transform.sigma = "none", transform.k = "none"), transform.sigma = "square", transform.k = "var", type.information = "observed") ## using sigma^2 and sigma^2 k^2
GS <- -hessian(func = function(p){p[c("sigma","k.M")] <- c(sqrt(p["sigma"]),sqrt(p["k.M"]/p["sigma"]));logLik(e.lmm, p = p, transform.sigma = "none", transform.k = "none")}, x = coef(e.lmm, transform.sigma = "square", transform.k = "var", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## *** variance-covariance
test <- attr(vcov(e.lmm, df = TRUE, transform.sigma = "square", transform.k = "var"),"df")
expect_equal(as.double(test[grep(names(test),pattern="M")]), c(sum(d$Gender=="M")/c(1,1,1,1,4)))
expect_equal(as.double(test[grep(names(test),pattern="F")]), c(sum(d$Gender=="F")/c(1,1,1,1,4)))

## *** variance-covariance
getVarCov(e.lmm)
getVarCov(e.lmm2)

## *** residuals
residuals(e.lmm, format = "long")
residuals(e.lmm, format = "wide")
residuals(e.lmm2)

residuals(e.lmm, format = "wide", type.residual = "normalized")

## *** confidence interval
confint(e.lmm, effects = "variance")
confint(e.lmm, effects = "variance", transform.sigma = "log", backtransform.sigma = "exp")


## *** print and summary
print(e.lmm)
summary(e.lmm)

## *** predictions
predict(e.lmm, newdata = d)

## *** interface to other packages
glht(e.lmm)
summary(e.lmm)
emmeans(e.lmm, specs = ~Gender)


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
