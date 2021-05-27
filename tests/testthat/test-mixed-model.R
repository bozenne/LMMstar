### test-mixed-model.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 14 2021 (16:46) 
## Version: 
## Last-Updated: May 27 2021 (12:27) 
##           By: Brice Ozenne
##     Update #: 34
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
## LMMstar.options(method.numDeriv = "Richardson")
level.test <- 1

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
## eCS.lmm <- lmm(Y ~ visit + age + gender, repetition = ~visit|id, structure = "CS", data = dL, debug = 2, method = "ML")
## eCS.gls <- gls(Y ~ visit + age + gender, correlation = corCompSymm(form=~1|id), data = dL, method = "ML")
eCS.lmm <- lmm(Y ~ visit + age + gender, repetition = ~visit|id, structure = "CS", data = dL, debug = 2, method = "REML")
eCS.gls <- gls(Y ~ visit + age + gender, correlation = corCompSymm(form=~1|id), data = dL, method = "REML")

## ** coef
expect_equal(coef(eCS.lmm, effects = "mean"), coef(eCS.gls), tol = 1e-6)
coef(eCS.lmm, transform.rho = "cov")
coef(eCS.lmm, transform.sigma = "square")

## ** logLikelihood
expect_equal(logLik(eCS.lmm), as.double(logLik(eCS.gls)), tol = 1e-6)

## ** score
GS <- jacobian(func = function(p){logLik(eCS.lmm, p = p)}, x = coef(eCS.lmm))
test <- score(eCS.lmm)
expect_true(all(abs(test) < 1e-5))
expect_equal(as.double(GS), as.double(test), tol = 1e-5)

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
expect_equal(as.double(test), as.double(GS), tol = 1e-6) 

## ** degree of freedom
test <- confint(eCS.lmm, transform.sigma = "none")[,"df",drop=FALSE]
expect_equal(test[c("visitY2","visitY3","visitY4"),"df"], rep(297,3), tol = 1e-3)
expect_equal(test[c("age","genderfemale"),"df"], rep(97,2), tol = 1e-3)

## ** anova
anova(eCS.lmm)

## ** getVarCov
getVarCov(eCS.lmm)

## * Unstructed covariance matrix
## ** fit
eUNexp.lmm <- lmm(Y ~ visit + age + gender, repetition = ~visit|id, structure = "UN", data = dL, debug = 2, method = "REML", type.information = "expected", df = level.test)
eUN.lmm <- lmm(Y ~ visit + age + gender, repetition = ~visit|id, structure = "UN", data = dL, debug = 2, method = "REML", type.information = "observed")
eUN.gls <- gls(Y ~ visit + age + gender, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|visit), data = dL, method = "REML")

## ** coef
expect_equal(coef(eUN.lmm, effects = "mean"), coef(eUN.gls), tol = 1e-6)
expect_equal(coef(eUNexp.lmm, effects = "mean"), coef(eUN.gls), tol = 1e-6)

## ** logLikelihood
expect_equal(logLik(eUN.lmm), as.double(logLik(eUN.gls)), tol = 1e-6)
expect_equal(logLik(eUNexp.lmm), as.double(logLik(eUN.gls)), tol = 1e-6)

## ** score
test <- score(eUN.lmm)
GS <- jacobian(func = function(p){logLik(eUN.lmm, p = p)}, x = coef(eUN.lmm))
## expect_true(all(abs(test) < 1e-3))
expect_equal(as.double(GS), as.double(test), tol = 1e-5)

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

test <- information(eUN.lmm, p = coef(eUN.lmm, transform.sigma = "none"), transform.sigma = "log")
GS <- -hessian(func = function(p){p["sigma"]<-exp(p["sigma"]);logLik(eUN.lmm, p = p)}, x = coef(eUN.lmm, transform.sigma = "log", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(eUN.lmm, p = coef(eUN.lmm, transform.sigma = "none"), transform.sigma = "square") 
GS <- -hessian(func = function(p){p["sigma"]<-sqrt(p["sigma"]);logLik(eUN.lmm, p = p)}, x = coef(eUN.lmm, transform.sigma = "square", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## test <- information(eUN.lmm, p = coef(eUN.lmm, transform.sigma = "none"), transform.rho = "cov") 
## GS <- -hessian(func = function(p){p[c("sigma","Rho")] <- c(sqrt(p["sigma"]),p["Rho"]/p["sigma"]); logLik(eUN.lmm, p = p)}, x = coef(eUN.lmm, transform.rho = "cov", transform.names = FALSE))
## expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## no transformation 
newp <- coef(eUN.lmm, transform.sigma = "none")+0.1
GS <- -hessian(func = function(p){logLik(eUN.lmm, p = p)}, x = newp)
test <- information(eUN.lmm, p = newp, transform.sigma = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6) 

## ** degree of freedom
test <- confint(eUN.lmm)[,"df", drop=FALSE]
## anova(eUN.lmm)
## anova(eUN.gls)

## ** anova
anova(eUN.lmm, print.null = TRUE)


## ** getVarCov
getVarCov(eUN.lmm)

## * Stratified random intercept model / Compound symmetry structure
## ** fit
eCS.lmm <- lmm(Y ~ (visit + age)*gender, repetition = gender~visit|id, structure = "CS", data = dL, debug = 2, method = "REML")
eCS.gls <- list(male=gls(Y ~ visit + age, correlation = corCompSymm(form=~1|id), data = dL[dL$gender=="male",], method = "REML"),
                female=gls(Y ~ visit + age, correlation = corCompSymm(form=~1|id), data = dL[dL$gender=="female",], method = "REML"))

## ** coef
expect_equal(unname(coef(eCS.lmm, effects = "mean", strata = "male")), unname(coef(eCS.gls$male)), tol = 1e-6)
expect_equal(unname(coef(eCS.lmm, effects = "mean", strata = "female")), unname(coef(eCS.gls$female)), tol = 1e-6)
coef(eCS.lmm, transform.rho = "cov")
coef(eCS.lmm, transform.sigma = "square")

## ** logLikelihood
expect_equal(logLik(eCS.lmm), as.double(logLik(eCS.gls$male)+logLik(eCS.gls$female)), tol = 1e-6)

## ** score
test <- score(eCS.lmm)
GS <- jacobian(func = function(p){logLik(eCS.lmm, p = p)}, x = coef(eCS.lmm))
expect_true(all(abs(test) < 1e-3))
expect_equal(as.double(GS), as.double(test), tol = 1e-5)

## no transformation
newp <- coef(eCS.lmm, transform.sigma = "none")+0.1
GS <- jacobian(func = function(p){logLik(eCS.lmm, p = p)}, x = newp)
test <- score(eCS.lmm, p = newp, transform.sigma = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log[c("sigma:male","sigma:female")] <- log(newp[c("sigma:male","sigma:female")])
GS <- jacobian(func = function(p){p[c("sigma:male","sigma:female")] <- exp(p[c("sigma:male","sigma:female")]); logLik(eCS.lmm, p = p)}, x = newp.log)
test <- score(eCS.lmm, p = newp, transform.sigma = "log")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2[c("sigma:male","sigma:female")] <- newp[c("sigma:male","sigma:female")]^2
GS <- jacobian(func = function(p){p[c("sigma:male","sigma:female")] <- sqrt(p[c("sigma:male","sigma:female")]); logLik(eCS.lmm, p = p)}, x = newp.2)
test <- score(eCS.lmm, p = newp, transform.sigma = "square")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ## cov transformation
## newp.cov <- newp.2; newp.cov[c("Rho:male","Rho:female")] <- c(newp["Rho:male"]*newp["sigma:male"]^2, newp["Rho:female"]*newp["sigma:female"]^2)
## GS <- jacobian(func = function(p){p[c("sigma:male","Rho:male","sigma:female","Rho:female")] <- c(sqrt(p["sigma:male"]),p["Rho:male"]/p["sigma:male"],p["sigma:female"],p["Rho:female"]/p["sigma:female"]); logLik(eCS.lmm, p = p)}, x = newp.cov)
## test <- score(eCS.lmm, p = newp, transform.rho = "cov")
## expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** information
test <- information(eCS.lmm, p = coef(eCS.lmm, transform.sigma = "none"), transform.sigma = "none")
GS <- -hessian(func = function(p){logLik(eCS.lmm, p = p)}, x = coef(eCS.lmm, transform.sigma = "none"))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(eCS.lmm, p = coef(eCS.lmm, transform.sigma = "none"), transform.sigma = "log")
GS <- -hessian(func = function(p){p[c("sigma:male","sigma:female")]<-exp(p[c("sigma:male","sigma:female")]);logLik(eCS.lmm, p = p)}, x = coef(eCS.lmm, transform.sigma = "log", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(eCS.lmm, p = coef(eCS.lmm, transform.sigma = "none"), transform.sigma = "square") 
GS <- -hessian(func = function(p){p[c("sigma:male","sigma:female")]<-sqrt(p[c("sigma:male","sigma:female")]);logLik(eCS.lmm, p = p)}, x = coef(eCS.lmm, transform.sigma = "square", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## test <- information(eCS.lmm, p = coef(eCS.lmm, transform.sigma = "none"), transform.rho = "cov") 
## GS <- -hessian(func = function(p){p[c("sigma","Rho")] <- c(sqrt(p["sigma"]),p["Rho"]/p["sigma"]); logLik(eCS.lmm, p = p)}, x = coef(eCS.lmm, transform.rho = "cov", transform.names = FALSE))
## expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## no transformation 
newp <- coef(eCS.lmm, transform.sigma = "none")+0.1
GS <- -hessian(func = function(p){logLik(eCS.lmm, p = p)}, x = newp)
test <- information(eCS.lmm, p = newp, transform.sigma = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6) 

## ** degree of freedom
test <- confint(eCS.lmm, transform.sigma = "none")[,"df",drop=FALSE]

## ** anova
anova(eCS.lmm, print.null = TRUE)

## ** getVarCov
getVarCov(eCS.lmm)

## * Stratified unstructed covariance matrix
## ** fit
eUN.lmm <- lmm(Y ~ (visit + age)*gender, repetition = gender~visit|id, structure = "UN", data = dL, debug = 2, method = "REML")
eUN.gls <- list(male=gls(Y ~ visit + age, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|visit), data = dL[dL$gender=="male",], method = "REML"),
                female=gls(Y ~ visit + age, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|visit), data = dL[dL$gender=="female",], method = "REML"))
summary(eUN.lmm)

## ** coef
expect_equal(unname(coef(eUN.lmm, effects = "mean", strata = "male")), unname(coef(eUN.gls$male)), tol = 1e-6)
expect_equal(unname(coef(eUN.lmm, effects = "mean", strata = "female")), unname(coef(eUN.gls$female)), tol = 1e-6)
coef(eUN.lmm, transform.rho = "cov")
coef(eUN.lmm, transform.k = "var")
coef(eUN.lmm, transform.k = "log")

## lava trans
fct.trans <- function(p, inverse = FALSE){
    newp <- p
    if(inverse){
        newp[c("sigma:male","k.Y2:male","k.Y3:male","k.Y4:male")] <- c(sqrt(newp["sigma:male"]), sqrt(newp[c("k.Y2:male","k.Y3:male","k.Y4:male")]/newp["sigma:male"]));
        newp[c("sigma:female","k.Y2:female","k.Y3:female","k.Y4:female")] <- c(sqrt(newp["sigma:female"]), sqrt(newp[c("k.Y2:female","k.Y3:female","k.Y4:female")]/newp["sigma:female"]));
    }else{
        newp[c("sigma:male","k.Y2:male","k.Y3:male","k.Y4:male")] <- p["sigma:male"]^2 * c(1,p[c("k.Y2:male","k.Y3:male","k.Y4:male")]^2);
        newp[c("sigma:female","k.Y2:female","k.Y3:female","k.Y4:female")] <- p["sigma:female"]^2 * c(1,p[c("k.Y2:female","k.Y3:female","k.Y4:female")]^2);
    }
    return(newp)
}
expect_equal(fct.trans(coef(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none")),
             coef(eUN.lmm, transform.k = "var", transform.names=FALSE))

## ** logLikelihood
expect_equal(logLik(eUN.lmm), as.double(logLik(eUN.gls$male)+logLik(eUN.gls$female)), tol = 1e-6)

## ** score
test <- score(eUN.lmm)
GS <- jacobian(func = function(p){logLik(eUN.lmm, p = p)}, x = coef(eUN.lmm))
expect_true(all(abs(test) < 1e-3))
expect_equal(as.double(GS), as.double(test), tol = 1e-5)

## no transformation
newp <- coef(eUN.lmm, transform.sigma = "none")+0.1
GS <- jacobian(func = function(p){logLik(eUN.lmm, p = p)}, x = newp)
test <- score(eUN.lmm, p = newp, transform.sigma = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log[eUN.lmm$param$type %in% "sigma"] <- log(newp[eUN.lmm$param$type %in% "sigma"])
GS <- jacobian(func = function(p){p[eUN.lmm$param$type %in% "sigma"] <- exp(p[eUN.lmm$param$type %in% "sigma"]); logLik(eUN.lmm, p = p)}, x = newp.log)
test <- score(eUN.lmm, p = newp, transform.sigma = "log")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- fct.trans(newp); 
GS <- jacobian(func = function(p){logLik(eUN.lmm, p = fct.trans(p,inverse = TRUE))}, x = newp.2)
test <- score(eUN.lmm, p = newp, transform.k = "var")
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

test <- information(eUN.lmm, p = coef(eUN.lmm, transform.sigma = "none"), transform.sigma = "log")
GS <- -hessian(func = function(p){p[eUN.lmm$param$type %in% "sigma"]<-exp(p[eUN.lmm$param$type %in% "sigma"]);logLik(eUN.lmm, p = p)}, x = coef(eUN.lmm, transform.sigma = "log", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(eUN.lmm, p = coef(eUN.lmm, transform.k = "none"), transform.k = "var") 
GS <- -hessian(func = function(p){logLik(eUN.lmm, p = fct.trans(p,inverse = TRUE))}, x = coef(eUN.lmm, transform.k = "var", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## test <- information(eUN.lmm, p = coef(eUN.lmm, transform.sigma = "none"), transform.rho = "cov") 
## GS <- -hessian(func = function(p){p[c("sigma","Rho")] <- c(sqrt(p["sigma"]),p["Rho"]/p["sigma"]); logLik(eUN.lmm, p = p)}, x = coef(eUN.lmm, transform.rho = "cov", transform.names = FALSE))
## expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## no transformation 
newp <- coef(eUN.lmm, transform.sigma = "none")+0.1
GS <- -hessian(func = function(p){logLik(eUN.lmm, p = p)}, x = newp)
test <- information(eUN.lmm, p = newp, transform.sigma = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6) 

## ** degree of freedom
test <- confint(eUN.lmm)[,"df", drop=FALSE]
## anova(eUN.lmm)
## anova(eUN.gls)

## ** anova
anova(eUN.lmm, print.null = TRUE)

## ** getVarCov
getVarCov(eUN.lmm)



##----------------------------------------------------------------------
### test-mixed-model.R ends here
