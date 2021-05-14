### test-mixed-model.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 14 2021 (16:46) 
## Version: 
## Last-Updated: May 14 2021 (17:43) 
##           By: Brice Ozenne
##     Update #: 6
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


## profvis::profvis(lmm(Y ~ visit + age + gender, variance = ~visit|id, structure = "CS", data = dL, debug = 2, method = "ML"))
eCS.lmm$param
coef(eCS.lmm)
summary(eCS.lmm)
confint(eCS.lmm)
anova(eCS.lmm)

## * Unstructed covariance matrix
## ** fit
eUN.lmm <- lmm(Y ~ visit + age + gender, variance = ~visit|id, structure = "UN", data = dL, debug = 2, method = "REML", df = FALSE, type.information = "expected")
eUN.gls <- gls(Y ~ visit + age + gender, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|visit), data = dL, method = "REML")

X <- model.matrix(eUN.lmm, effect = "mean")
X.sum <- Reduce("+",tapply(1:NROW(dL),dL$id,function(i){X[i,]}))
Omega <- getVarCov(eUN.lmm)
solve(t(X.sum) %*% solve(Omega) %*% X.sum)
information(eUN.lmm)

logLik(eUN.lmm)
logLik(eUN.gls)
as.double(getVarCov(eUN.lmm))-as.double(getVarCov(eUN.gls))
summary(eUN.lmm)
eUN.lmm$df

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
### test-mixed-model.R ends here
