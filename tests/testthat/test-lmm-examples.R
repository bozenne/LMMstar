### test-lmm-examples.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (10:13) 
## Version: 
## Last-Updated: mar 27 2021 (00:00) 
##           By: Brice Ozenne
##     Update #: 10
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


## ** single variance parameter
e.lmm <- lmm(Y ~ X1 + X2 + X3, variance = ~time|id, structure = "CS", data = d, debug = 2,
             method = "REML")

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
