### test-mixed-model.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 14 2021 (16:46) 
## Version: 
## Last-Updated: okt 15 2021 (17:35) 
##           By: Brice Ozenne
##     Update #: 98
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

context("Check lmm on examples of mixed model")
LMMstar.options(optimizer = "gls", method.numDeriv = "simple", precompute.moments = TRUE, # "Richardson"
                columns.confint = c("estimate","se","df","lower","upper","p.value"))

## * Simulate data
m <- lvm(c(Y1,Y2,Y3,Y4) ~ 0.05*age + gender)
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

test_that("Compound symmetry structure (REML)",{

## ** fit
eCS.lmm <- lmm(Y ~ visit + age + gender, repetition = ~visit|id, structure = "CS", data = dL, trace = 0, method = "REML")
eCS.gls <- gls(Y ~ visit + age + gender, correlation = corCompSymm(form=~1|id), data = dL, method = "REML")

## ** coef
expect_equal(coef(eCS.lmm, effects = "mean"), coef(eCS.gls), tol = 1e-6)
coef(eCS.lmm, transform.rho = "cov")
coef(eCS.lmm, transform.sigma = "square")

## ** logLikelihood
expect_equal(logLik(eCS.lmm), as.double(logLik(eCS.gls)), tol = 1e-6)

## ** score
GS <- jacobian(func = function(p){logLik(eCS.lmm, p = p)}, x = coef(eCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"))
test <- score(eCS.lmm, effects = "all")
expect_true(all(abs(test) < 1e-5))
expect_equal(as.double(GS), as.double(test), tol = 1e-5)

## no transformation
newp <- coef(eCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")+0.1
GS <- jacobian(func = function(p){logLik(eCS.lmm, p = p)}, x = newp)
test <- score(eCS.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(eCS.lmm, p = p)}, x = newp.log)
test <- score(eCS.lmm, effects = "all", p = newp , transform.sigma = "log", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2
GS <- jacobian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(eCS.lmm, p = p)}, x = newp.2)
test <- score(eCS.lmm, effects = "all", p = newp, transform.sigma = "square", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## cov transformation
newp.cov <- newp; newp.cov[c("sigma","rho")] <- c(newp["sigma"]^2,newp["rho"]*newp["sigma"]^2)
GS <- jacobian(func = function(p){p[c("sigma","rho")] <- c(sqrt(p["sigma"]),p["rho"]/p["sigma"]); logLik(eCS.lmm, p = p, transform.sigma = "none", transform.k = "none", transform.rho = "cov")}, x = newp.cov)
test <- score(eCS.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none", transform.rho = "cov")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** information
## no transformation
test <- information(eCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none",
                    p = coef(eCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"))
GS <- -hessian(func = function(p){logLik(eCS.lmm, p = p)},
               x = coef(eCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"])
GS <- -hessian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(eCS.lmm, p = p)}, x = newp.log)
test <- information(eCS.lmm, p = newp, effects = "all", transform.sigma = "log", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2
GS <- -hessian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(eCS.lmm, p = p)}, x = newp.2)
test <- information(eCS.lmm, p = newp, effects = "all", transform.sigma = "square", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## cov transformation
newp.cov <- newp; newp.cov[c("sigma","rho")] <- c(newp["sigma"]^2,newp["rho"]*newp["sigma"]^2)
GS <- -hessian(func = function(p){p[c("sigma","rho")] <- c(sqrt(p["sigma"]),p["rho"]/p["sigma"]); logLik(eCS.lmm, p = p)}, x = newp.cov)
test <- information(eCS.lmm, p = newp, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "cov")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## no transformation 
newp <- coef(eCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")+0.1
GS <- -hessian(func = function(p){logLik(eCS.lmm, p = p)}, x = newp)
test <- information(eCS.lmm, p = newp, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6) 

## ** degree of freedom
test <- confint(eCS.lmm, effects = "all", transform.sigma = "none")[,"df",drop=FALSE]
expect_equal(test[c("visitY2","visitY3","visitY4"),"df"], rep(297,3), tol = 1e-3)
expect_equal(test[c("age","genderfemale"),"df"], rep(97,2), tol = 1e-3)

## ** anova
eCS.lmm_anova <- anova(eCS.lmm, effects = "all")
expect_equal(eCS.lmm_anova$mean$df.denom, c(297.03106,  97.05084,  97.05084), tol = 1e-1)
expect_equal(eCS.lmm_anova$correlation$df.denom, c(14.7493), tol = 1e-1)

## ** getVarCov
getVarCov(eCS.lmm)

## ** prediction
test <- predict(eCS.lmm, newdata = dL)
index <- sample.int(NROW(dL))
test2 <- predict(eCS.lmm, newdata = dL[index,,drop=FALSE])
if(require(AICcmodavg)){
    GS <- AICcmodavg::predictSE(eCS.gls, newdata = dL)
    expect_equivalent(test$estimate,GS$fit, tol = 1e-7)
    expect_equivalent(test$se,GS$se.fit, tol = 1e-7)
    expect_equivalent(test[index,,drop=FALSE],test2, tol = 1e-7)
}

## ** ICC
## eCS0.lmm <- lmm(Y ~ 1, repetition = ~visit|id, structure = "CS", data = dL[dL$visit %in% c("Y1","Y2"),], trace = 0, method = "REML")
## psych::ICC(dW[,c("Y1","Y2")])$lme
## confint(eCS0.lmm, effect = "correlation")
})

## * Unstructed covariance matrix
test_that("Unstructured covariance matrix (REML)",{

## ** fit
eUNexp.lmm <- lmm(Y ~ visit + age + gender, repetition = ~visit|id, structure = "UN", data = dL, trace = 0, method = "REML", type.information = "expected")
eUN.lmm <- lmm(Y ~ visit + age + gender, repetition = ~visit|id, structure = "UN", data = dL, trace = 0, method = "REML", type.information = "observed")
eUN.gls <- gls(Y ~ visit + age + gender, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|visit), data = dL, method = "REML")

## ** coef
expect_equal(coef(eUN.lmm, effects = "mean"), coef(eUN.gls), tol = 1e-6)
expect_equal(coef(eUNexp.lmm, effects = "mean"), coef(eUN.gls), tol = 1e-6)

## ** logLikelihood
expect_equal(logLik(eUN.lmm), as.double(logLik(eUN.gls)), tol = 1e-6)
expect_equal(logLik(eUNexp.lmm), as.double(logLik(eUN.gls)), tol = 1e-6)

## ** score
test <- score(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")
GS <- jacobian(func = function(p){logLik(eUN.lmm, p = p)},
               x = coef(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"))
## expect_true(all(abs(test) < 1e-3))
expect_equal(as.double(GS), as.double(test), tol = 1e-5)

## no transformation
newp <- coef(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")+0.1
GS <- jacobian(func = function(p){logLik(eUN.lmm, p = p)}, x = newp)
test <- score(eUN.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(eUN.lmm, p = p)}, x = newp.log)
test <- score(eUN.lmm, effects = "all", p = newp, transform.sigma = "log", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2
GS <- jacobian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(eUN.lmm, p = p)}, x = newp.2)
test <- score(eUN.lmm, effects = "all", p = newp, transform.sigma = "square", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## cov transformation
## newp.cov <- newp; newp.cov[c("sigma","rho")] <- c(newp["sigma"]^2,newp["rho"]*newp["sigma"]^2)
## GS <- jacobian(func = function(p){p[c("sigma","rho")] <- c(sqrt(p["sigma"]),p["rho"]/p["sigma"]); logLik(eUN.lmm, p = p)}, x = newp.cov)
## test <- score(eUN.lmm, p = newp, transform.rho = "cov")
## expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** information
expect_equal(vcov(eUNexp.lmm, effects = "mean"), vcov(eUN.gls), tol = 1e-6)

## no transformation
newp <- coef(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")+0.1
GS <- -hessian(func = function(p){logLik(eUN.lmm, p = p)}, x = newp)
test <- information(eUN.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"])
GS <- -hessian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(eUN.lmm, p = p)}, x = newp.log)
test <- information(eUN.lmm, effects = "all", p = newp, transform.sigma = "log", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2
GS <- -hessian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(eUN.lmm, p = p)}, x = newp.2)
test <- information(eUN.lmm, effects = "all", p = newp, transform.sigma = "square", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## cov transformation

## no transformation 
newp <- coef(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")+0.1
GS <- -hessian(func = function(p){logLik(eUN.lmm, p = p)}, x = newp)
test <- information(eUN.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6) 

## ** degree of freedom and anova checked in test-ttest.R (section multiple t-test)
## anova(eUN.lmm, ci = TRUE)

eUN.lmm_anova <- anova(eUN.lmm, effects = "all", ci = TRUE)
expect_equal(eUN.lmm_anova$mean$df.denom, c(99.00144, 96.18348, 94.17755), tol = 1e-1)
expect_equal(eUN.lmm_anova$variance$df.denom, c(189.236), tol = 1e-1)
expect_equal(eUN.lmm_anova$correlation$df.denom, c(20.66968), tol = 1e-1)

## ** getVarCov
getVarCov(eUN.lmm)
})

## * Stratified random intercept model / Compound symmetry structure
test_that("Stratified compound symmetry structure (REML)",{

## ** fit
eCS0.lmm <- lmm(Y ~ gender, repetition = gender~visit|id, structure = "CS", data = dL, trace = 0, method = "REML") ## just to check that it runs
eCS.lmm <- lmm(Y ~ (visit + age)*gender, repetition = gender~visit|id, structure = "CS", data = dL, trace = 0, method = "REML")
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
test <- score(eCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")
GS <- jacobian(func = function(p){logLik(eCS.lmm, p = p)},
               x = coef(eCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"))
expect_true(all(abs(test) < 1e-3))
expect_equal(as.double(GS), as.double(test), tol = 1e-5)

## no transformation
newp <- coef(eCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")+0.1
GS <- jacobian(func = function(p){logLik(eCS.lmm, p = p)}, x = newp)
test <- score(eCS.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)


## log transformation
newp.log <- newp; newp.log[c("sigma:male","sigma:female")] <- log(newp[c("sigma:male","sigma:female")])
GS <- jacobian(func = function(p){p[c("sigma:male","sigma:female")] <- exp(p[c("sigma:male","sigma:female")]); logLik(eCS.lmm, p = p)}, x = newp.log)
test <- score(eCS.lmm, effects = "all", p = newp, transform.sigma = "log", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2[c("sigma:male","sigma:female")] <- newp[c("sigma:male","sigma:female")]^2
GS <- jacobian(func = function(p){p[c("sigma:male","sigma:female")] <- sqrt(p[c("sigma:male","sigma:female")]); logLik(eCS.lmm, p = p)}, x = newp.2)
test <- score(eCS.lmm, effects = "all", p = newp, transform.sigma = "square", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ## cov transformation
## newp.cov <- newp.2; newp.cov[c("rho:male","rho:female")] <- c(newp["rho:male"]*newp["sigma:male"]^2, newp["rho:female"]*newp["sigma:female"]^2)
## GS <- jacobian(func = function(p){p[c("sigma:male","rho:male","sigma:female","rho:female")] <- c(sqrt(p["sigma:male"]),p["rho:male"]/p["sigma:male"],p["sigma:female"],p["rho:female"]/p["sigma:female"]); logLik(eCS.lmm, p = p)}, x = newp.cov)
## test <- score(eCS.lmm, p = newp, transform.rho = "cov")
## expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** information
test <- information(eCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")
GS <- -hessian(func = function(p){logLik(eCS.lmm, p = p)},
               x = coef(eCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"))
expect_equal(as.double(GS), as.double(test), tol = 1e-5)

## no transformation
newp <- coef(eCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")+0.1
GS <- -hessian(func = function(p){logLik(eCS.lmm, p = p)}, x = newp)
test <- information(eCS.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log[c("sigma:male","sigma:female")] <- log(newp[c("sigma:male","sigma:female")])
GS <- -hessian(func = function(p){p[c("sigma:male","sigma:female")] <- exp(p[c("sigma:male","sigma:female")]); logLik(eCS.lmm, p = p)}, x = newp.log)
test <- information(eCS.lmm, effects = "all", p = newp, transform.sigma = "log", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2[c("sigma:male","sigma:female")] <- newp[c("sigma:male","sigma:female")]^2
GS <- -hessian(func = function(p){p[c("sigma:male","sigma:female")] <- sqrt(p[c("sigma:male","sigma:female")]); logLik(eCS.lmm, p = p)}, x = newp.2)
test <- information(eCS.lmm, effects = "all", p = newp, transform.sigma = "square", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** degree of freedom
test <- confint(eCS.lmm, effects = "all", transform.sigma = "none")[,"df",drop=FALSE]

## ** anova
eCS.lmm_anova <- anova(eCS.lmm, effects = "all")
expect_equal(eCS.lmm_anova$mean$df.denom, c(143.04097,  46.63249), tol = 1e-1)
expect_equal(eCS.lmm_anova$correlation$df.denom, c(7.450136), tol = 1e-1)

## ** getVarCov
getVarCov(eCS.lmm)
})

## * Stratified unstructed covariance matrix
test_that("Stratified unstructured (REML)",{

## ** fit
eUN.lmm <- lmm(Y ~ (visit + age)*gender, repetition = gender~visit|id, structure = "UN", data = dL, trace = 0, method = "REML")
eUN.gls <- list(male=gls(Y ~ visit + age, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|visit), data = dL[dL$gender=="male",], method = "REML"),
                female=gls(Y ~ visit + age, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|visit), data = dL[dL$gender=="female",], method = "REML"))
capture.output(summary(eUN.lmm))
## anova(eUN.lmm, effects = "all")

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
expect_equal(fct.trans(coef(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")),
             coef(eUN.lmm, effects = "all", transform.k = "var", transform.rho = "none", transform.names=FALSE))

## ** logLikelihood
expect_equal(logLik(eUN.lmm), as.double(logLik(eUN.gls$male)+logLik(eUN.gls$female)), tol = 1e-6)

## ** score
test <- score(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")
GS <- jacobian(func = function(p){logLik(eUN.lmm, p = p)}, x = coef(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"))
expect_true(all(abs(test) < 1e-3))
expect_equal(as.double(GS), as.double(test), tol = 1e-5)

## no transformation
newp <- coef(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")+0.1
GS <- jacobian(func = function(p){logLik(eUN.lmm, p = p)}, x = newp)
test <- score(eUN.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log[eUN.lmm$param$type %in% "sigma"] <- log(newp[eUN.lmm$param$type %in% "sigma"])
GS <- jacobian(func = function(p){p[eUN.lmm$param$type %in% "sigma"] <- exp(p[eUN.lmm$param$type %in% "sigma"]); logLik(eUN.lmm, p = p)}, x = newp.log)
test <- score(eUN.lmm, effects = "all", p = newp, transform.sigma = "log", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- fct.trans(newp); 
GS <- jacobian(func = function(p){logLik(eUN.lmm, p = fct.trans(p,inverse = TRUE))}, x = newp.2)
test <- score(eUN.lmm, effects = "all", p = newp, transform.k = "var", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## cov transformation
## newp.cov <- newp; newp.cov[c("sigma","rho")] <- c(newp["sigma"]^2,newp["rho"]*newp["sigma"]^2)
## GS <- jacobian(func = function(p){p[c("sigma","rho")] <- c(sqrt(p["sigma"]),p["rho"]/p["sigma"]); logLik(eUN.lmm, p = p)}, x = newp.cov)
## test <- score(eUN.lmm, p = newp, transform.rho = "cov")
## expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** information
## expect_equal(as.double(vcov(eUN.lmm, effects = "mean")), as.double(Matrix::bdiag(lapply(eUN.gls,vcov))), tol = 1e-6)

test <- information(eUN.lmm, effects = "all",
                    p = coef(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"),
                    transform.sigma = "none", transform.k = "none", transform.rho = "none")
GS <- -hessian(func = function(p){logLik(eUN.lmm, p = p)},
               x = coef(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"))
expect_equal(as.double(test),as.double(GS), tol = 1e-5)

test <- information(eUN.lmm, effects = "all",
                    p = coef(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"),
                    transform.sigma = "log", transform.k = "none", transform.rho = "none")
GS <- -hessian(func = function(p){p[eUN.lmm$param$type %in% "sigma"]<-exp(p[eUN.lmm$param$type %in% "sigma"]);logLik(eUN.lmm, p = p)},
               x = coef(eUN.lmm, effects = "all", transform.sigma = "log", transform.k = "none", transform.rho = "none", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(eUN.lmm, effects = "all",
                    p = coef(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"),
                    transform.k = "var", transform.rho = "none") 
GS <- -hessian(func = function(p){logLik(eUN.lmm, p = fct.trans(p,inverse = TRUE))},
               x = coef(eUN.lmm, effects = "all", transform.k = "var", transform.rho = "none", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## test <- information(eUN.lmm, p = coef(eUN.lmm, transform.sigma = "none"), transform.rho = "cov") 
## GS <- -hessian(func = function(p){p[c("sigma","rho")] <- c(sqrt(p["sigma"]),p["rho"]/p["sigma"]); logLik(eUN.lmm, p = p)}, x = coef(eUN.lmm, transform.rho = "cov", transform.names = FALSE))
## expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## no transformation 
newp <- coef(eUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")+0.1
GS <- -hessian(func = function(p){logLik(eUN.lmm, p = p)}, x = newp)
test <- information(eUN.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6) 

## ** degree of freedom
test <- confint(eUN.lmm, effects = "all")[,"df", drop=FALSE]
## anova(eUN.lmm)
## anova(eUN.gls)

## ** anova
eUN.lmm_anova <- anova(eUN.lmm, effects = "all", ci = TRUE)
expect_equal(eUN.lmm_anova$mean$df.denom, c(47.63744, 42.53266), tol = 1e-1)
expect_equal(eUN.lmm_anova$variance$df.denom, c(86.07826), tol = 1e-1)
expect_equal(eUN.lmm_anova$correlation$df.denom, c(9.100919), tol = 1e-1)


## ** getVarCov
getVarCov(eUN.lmm)
})

## * Missing data
test_that("missing values",{

## ** full cluster missing
    dL$Ymiss <- dL$Y
    dL$Ymiss[dL$id==1] <- NA

    eCS.lmm <- lmm(Ymiss ~ visit + age + gender, repetition = ~visit|id, structure = "CS", data = dL, trace = 0, method.fit = "ML")
    eCS2.lmm <- lmm(Ymiss ~ visit + age + gender, repetition = ~visit|id, structure = "CS", data = dL[dL$id!=1,], trace = 0, method.fit = "ML")
    expect_equal(logLik(eCS2.lmm), logLik(eCS.lmm))
    
    ## logLik(eCS.lmm, indiv = TRUE)
    ## score(eCS.lmm, indiv = TRUE)
    ## information(eCS.lmm, indiv = TRUE)
    
## ** only part of the cluster is missing
    set.seed(11)
    dL$Ymiss <- dL$Y
    dL$Ymiss[which(rbinom(NROW(dL), size = 1, prob = 0.1)==1)] <- NA

    eCS.lmm <- lmm(Ymiss ~ visit + age + gender, repetition = ~visit|id, structure = "CS", data = dL, trace = 0, method.fit = "ML")
    eCS.gls <- gls(Ymiss ~ visit + age + gender, correlation = corCompSymm(form=~1|id), data = dL, method = "ML", na.action = na.omit)
    expect_equal(as.double(logLik(eCS.lmm)), as.double(logLik(eCS.gls)))

    dL[dL$id==dL$id[1],,drop=FALSE]
})
##----------------------------------------------------------------------
### test-mixed-model.R ends here
