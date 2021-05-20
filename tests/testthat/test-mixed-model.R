### test-mixed-model.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 14 2021 (16:46) 
## Version: 
## Last-Updated: May 20 2021 (13:22) 
##           By: Brice Ozenne
##     Update #: 20
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
    library(lme4)
    library(lmerTest)

    library(LMMstar)
}

context("Check lmm on examples of mixed model")
LMMstar.options(method.numDeriv = "simple")

## * Simulate data
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
 

## * Random intercept model / Compound symmetry structure
## ** fit
eCS.lmm <- lmm(Y ~ visit + age + gender, variance = ~visit|id, structure = "CS", data = dL, debug = 2, method = "REML")
eCS.gls <- gls(Y ~ visit + age + gender, correlation = corCompSymm(form=~1|id), data = dL, method = "REML")

## ** coef
expect_equal(coef(eCS.lmm, effects = "mean"), coef(eCS.gls), tol = 1e-6)
coef(eCS.lmm, transform.rho = "cov")
coef(eCS.lmm, transform.sigma = "square")

## ** logLikelihood
expect_equal(logLik(eCS.lmm), as.double(logLik(eCS.gls)), tol = 1e-6)

## ** score
expect_true(all(abs(score(eCS.lmm)) < 1e-5))

## no transformation
newp <- coef(eCS.lmm, transform.sigma = "none")+0.1
GS <- jacobian(func = function(p){logLik(eCS.lmm, p = p)}, x = newp)
test <- score(eCS.lmm, p = newp, transform.sigma = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(eCS.lmm, p = p)}, x = newp.log)
test <- score(eCS.lmm, p = newp, transform.sigma = "log")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2
GS <- jacobian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(eCS.lmm, p = p)}, x = newp.2)
test <- score(eCS.lmm, p = newp, transform.sigma = "square")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## cov transformation
newp.cov <- newp; newp.cov[c("sigma","Rho")] <- c(newp["sigma"]^2,newp["Rho"]*newp["sigma"]^2)
GS <- jacobian(func = function(p){p[c("sigma","Rho")] <- c(sqrt(p["sigma"]),p["Rho"]/p["sigma"]); logLik(eCS.lmm, p = p)}, x = newp.cov)
test <- score(eCS.lmm, p = newp, transform.rho = "cov")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** information
test <- information(eCS.lmm, p = coef(eCS.lmm, transform.sigma = "none"), transform.sigma = "none")
GS <- -hessian(func = function(p){logLik(eCS.lmm, p = p)}, x = coef(eCS.lmm, transform.sigma = "none"))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(eCS.lmm, p = coef(eCS.lmm, transform.sigma = "none"), transform.sigma = "log")
GS <- -hessian(func = function(p){p["sigma"]<-exp(p["sigma"]);logLik(eCS.lmm, p = p)}, x = coef(eCS.lmm, transform.sigma = "log", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(eCS.lmm, p = coef(eCS.lmm, transform.sigma = "none"), transform.sigma = "square") 
GS <- -hessian(func = function(p){p["sigma"]<-sqrt(p["sigma"]);logLik(eCS.lmm, p = p)}, x = coef(eCS.lmm, transform.sigma = "square", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(eCS.lmm, p = coef(eCS.lmm, transform.sigma = "none"), transform.rho = "cov") 
GS <- -hessian(func = function(p){p[c("sigma","Rho")] <- c(sqrt(p["sigma"]),p["Rho"]/p["sigma"]); logLik(eCS.lmm, p = p)}, x = coef(eCS.lmm, transform.rho = "cov", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## no transformation 
newp <- coef(eCS.lmm, transform.sigma = "none")+0.1
GS <- -hessian(func = function(p){logLik(eCS.lmm, p = p)}, x = newp)
test <- information(eCS.lmm, p = newp, transform.sigma = "none")
expect_equal(as.double(test[1:6,1:6]), as.double(GS[1:6,1:6]), tol = 1e-6) ## does not match as some terms do not cancel

## ** variance-covariance
test <- vcov(eCS.lmm, p = coef(eCS.lmm, transform.sigma = "none"), transform.sigma = "none")
GS <- solve(information(eCS.lmm, p = coef(eCS.lmm, transform.sigma = "none"), transform.sigma = "none"))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** degree of freedom
test <- confint(eCS.lmm, transform.sigma = "square")$df
## anova(eCS.lmm)
## anova(eCS.gls)

## ** getVarCov
getVarCov(eCS.lmm)

## * Unstructed covariance matrix
## ** fit
eUNexp.lmm <- lmm(Y ~ visit + age + gender, variance = ~visit|id, structure = "UN", data = dL, debug = 2, method = "REML", type.information = "expected")
eUN.lmm <- lmm(Y ~ visit + age + gender, variance = ~visit|id, structure = "UN", data = dL, debug = 2, method = "REML", type.information = "observed")
eUN.gls <- gls(Y ~ visit + age + gender, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|visit), data = dL, method = "REML")

## ** coef
expect_equal(coef(eUN.lmm, effects = "mean"), coef(eUN.gls), tol = 1e-6)
expect_equal(coef(eUNexp.lmm, effects = "mean"), coef(eUN.gls), tol = 1e-6)

## ** logLikelihood
expect_equal(logLik(eUN.lmm), as.double(logLik(eUN.gls)), tol = 1e-6)
expect_equal(logLik(eUNexp.lmm), as.double(logLik(eUN.gls)), tol = 1e-6)

## ** score
## expect_true(all(abs(score(eUN.lmm)) < 1e-5))
## expect_true(all(abs(score(eUN.lmm, transform.sigma = "log", transform.k = "log", transform.rho = "atanh")) < 1e-5))

## no transformation
newp <- coef(eUN.lmm, transform.sigma = "none")+0.1
GS <- jacobian(func = function(p){logLik(eUN.lmm, p = p)}, x = newp)
test <- score(eUN.lmm, p = newp, transform.sigma = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(eUN.lmm, p = p)}, x = newp.log)
test <- score(eUN.lmm, p = newp, transform.sigma = "log")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2
GS <- jacobian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(eUN.lmm, p = p)}, x = newp.2)
test <- score(eUN.lmm, p = newp, transform.sigma = "square")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## cov transformation
## newp.cov <- newp; newp.cov[c("sigma","Rho")] <- c(newp["sigma"]^2,newp["Rho"]*newp["sigma"]^2)
## GS <- jacobian(func = function(p){p[c("sigma","Rho")] <- c(sqrt(p["sigma"]),p["Rho"]/p["sigma"]); logLik(eUN.lmm, p = p)}, x = newp.cov)
## test <- score(eUN.lmm, p = newp, transform.rho = "cov")
## expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** information
expect_equal(vcov(eUNexp.lmm, effects = "mean"), vcov(eUN.gls), tol = 1e-6)
test <- information(eUN.lmm, p = coef(eUN.lmm, transform.sigma = "none"), transform.sigma = "none")
GS <- -hessian(func = function(p){logLik(eUN.lmm, p = p)}, x = coef(eUN.lmm, transform.sigma = "none"))
expect_equal(unname(test),GS, tol = 1e-5)

range(test[8:16,8:16]-GS2[8:16,8:16])
round(unname(test),2)
round(GS,2)
round(GS2,2)

range(unname(round(100*(test-GS)/abs(GS),2))[8:16,8:16])
range(unname(round(test-GS,4)))

GS2 <- -jacobian(func = function(p){score(eUN.lmm, p = p)}, x = coef(eUN.lmm, transform.sigma = "none"))
range(GS-GS2)
range(test-GS)
round(GS2,2)

## *** debug
coef(eUN.lmm, transform.rho = "cov")
ls.grad <- lapply(1:16, function(i){
    pp <- coef(eUN.lmm, effects = c("variance","correlation"))
    iGrad <- jacobian(func = function(p){getVarCov(eUN.lmm, p = p)[i]},  x = pp)
    colnames(iGrad) <- names(pp)
    return(iGrad)
})

ls.grad
ls.GS <- lapply(1:10,function(iVar){
    matrix(sapply(ls.grad,function(i){i[,iVar]}),4,4)
})



ls.hess <- lapply(1:16, function(i){
    pp <- coef(eUN.lmm, effects = c("variance","correlation"))
    iHess <- hessian(func = function(p){getVarCov(eUN.lmm, p = p)[i]},  x = pp)
    dimnames(iHess) <- list(names(pp),names(pp))
    return(iHess)
})
ls.GS <- lapply(1:n.pair,function(iVar){
    matrix(sapply(ls.hess,function(i){i[pair[1,iVar],pair[2,iVar]]}),4,4)
})

range(do.call(rbind,lapply(1:length(out[[1]]),function(iIndex){range(out[[1]][[iIndex]] - ls.GS[[iIndex]])})))

test <- information(eUN.lmm, p = coef(eUN.lmm, transform.sigma = "none"), transform.sigma = "none")
expand.grid(names(coef(eUN.lmm, effects = c("variance","correlation"))),
            names(coef(eUN.lmm, effects = c("variance","correlation"))))
ls.derivative <- lapply(1:16, function(i){matrix(sapply(ls.hess,function(i){i[1,1]}),4,4)})

##         [,1]    [,2]    [,3]    [,4]
## [1,]  2.0000 -0.1053  0.0167  0.0975
## [2,] -0.1053  2.1336 -0.2790  0.0835
## [3,]  0.0167 -0.2790  2.1234 -0.0140
## [4,]  0.0975  0.0835 -0.0140  2.4544
round(matrix(sapply(ls.hess,function(i){i[1,2]}),4,4),4)
##         [,1]    [,2]    [,3]   [,4]
## [1,]  0.0000 -0.0962  0.0000 0.0000
## [2,] -0.0962  3.8961 -0.2547 0.0762
## [3,]  0.0000 -0.2547  0.0000 0.0000
## [4,]  0.0000  0.0762  0.0000 0.0000

2*diag(getVarCov(eUN.lmm))/coef(eUN.lmm)["k.Y4"]
getVarCov(eUN.lmm)/coef(eUN.lmm)["cor(Y4,Y3)"]
getVarCov(eUN.lmm)


fields::image.plot(test-GS)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(eUN.lmm, p = coef(eUN.lmm, transform.sigma = "none"), transform.sigma = "log")
GS <- -hessian(func = function(p){p["sigma"]<-exp(p["sigma"]);logLik(eUN.lmm, p = p)}, x = coef(eUN.lmm, transform.sigma = "log", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(eUN.lmm, p = coef(eUN.lmm, transform.sigma = "none"), transform.sigma = "square") 
GS <- -hessian(func = function(p){p["sigma"]<-sqrt(p["sigma"]);logLik(eUN.lmm, p = p)}, x = coef(eUN.lmm, transform.sigma = "square", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(eUN.lmm, p = coef(eUN.lmm, transform.sigma = "none"), transform.rho = "cov") 
GS <- -hessian(func = function(p){p[c("sigma","Rho")] <- c(sqrt(p["sigma"]),p["Rho"]/p["sigma"]); logLik(eUN.lmm, p = p)}, x = coef(eUN.lmm, transform.rho = "cov", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## no transformation 
newp <- coef(eUN.lmm, transform.sigma = "none")+0.1
GS <- -hessian(func = function(p){logLik(eUN.lmm, p = p)}, x = newp)
test <- information(eUN.lmm, p = newp, transform.sigma = "none")
expect_equal(as.double(test[1:6,1:6]), as.double(GS[1:6,1:6]), tol = 1e-6) ## does not match as some terms do not cancel

## ** variance-covariance
test <- vcov(eUN.lmm, p = coef(eUN.lmm, transform.sigma = "none"), transform.sigma = "none")
GS <- solve(information(eUN.lmm, p = coef(eUN.lmm, transform.sigma = "none"), transform.sigma = "none"))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** degree of freedom
test <- confint(eUN.lmm, transform.sigma = "square")$df
## anova(eUN.lmm)
## anova(eUN.gls)

## ** anova
anova(eUN.lmm, print.null = TRUE)
anova(eUN.lmm, transform.sigma = "none", transform.k = "none", print.null = TRUE)


## ** getVarCov
getVarCov(eUN.lmm)






##----------------------------------------------------------------------
### test-mixed-model.R ends here
