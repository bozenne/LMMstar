### test-auto-mixed-model.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 14 2021 (16:46) 
## Version: 
## Last-Updated: jul 11 2023 (15:33) 
##           By: Brice Ozenne
##     Update #: 182
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
LMMstar.options(optimizer = "FS", method.numDeriv = "simple", precompute.moments = TRUE, # "Richardson"
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
 

## * Compound symmetry structure

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
expect_true(all(abs(test) < LMMstar.options()$tol.score))
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
test <- model.tables(eCS.lmm, effects = "all", transform.sigma = "none")[,"df",drop=FALSE]
expect_equal(test[c("visitY2","visitY3","visitY4"),"df"], rep(297,3), tol = 1e-3)
expect_equal(test[c("age","genderfemale"),"df"], rep(97,2), tol = 1e-3)

## ** anova
eCS.lmm_anova <- anova(eCS.lmm, effects = "all")
expect_equal(eCS.lmm_anova$multivariate[eCS.lmm_anova$multivariate$type=="mu","df.denom"], c(297.03106,  97.05084,  97.05084), tol = 1e-1)
expect_equal(eCS.lmm_anova$multivariate[eCS.lmm_anova$multivariate$type=="rho","df.denom"], c(14.7493), tol = 1e-1)

## ** getVarCov
sigma(eCS.lmm)

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

})

## * Unstructed covariance matrix
test_that("Unstructured covariance matrix (REML)",{

## ** fit
eUNexp.lmm <- lmm(Y ~ visit + age + gender, repetition = ~visit|id, structure = "UN", data = dL, trace = 0, method = "REML", type.information = "expected")
eUN.lmm <- lmm(Y ~ visit + age + gender, repetition = ~visit|id, structure = "UN", data = dL, trace = 0, method = "REML", type.information = "observed")
eUN.gls <- gls(Y ~ visit + age + gender, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|visit), data = dL, method = "REML")

#p# ** coef
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
## name.sigma <- eUN.lmm$design$param[eUN.lmm$design$param$type=="rho","name"]
## name.rho <- eUN.lmm$design$param[eUN.lmm$design$param$type=="rho","name"]
## name.k1 <- eUN.lmm$design$param[eUN.lmm$design$param$type=="rho","k.x"]
## name.k2 <- eUN.lmm$design$param[eUN.lmm$design$param$type=="rho","k.y"]
## newp.cov <- newp
## newp.cov[name.rho] <- newp[name.rho]*ifelse(is.na(newp[name.k1]),1,newp[name.k1])*newp[name.k2]*newp["sigma"]^2
## newp.cov[name.rho] <- newp[name.rho]*ifelse(is.na(newp[name.k1]),1,newp[name.k1])*newp[name.k2]*newp["sigma"]^2
## GS <- jacobian(func = function(p){p[name.rho] <- c(p[name.rho]*); logLik(eUN.lmm, p = p)}, x = newp.cov)
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

eUN.lmm_anova <- anova(eUN.lmm, effects = "all", ci = TRUE)$multivariate
expect_equal(eUN.lmm_anova[eUN.lmm_anova$type=="mu","df.denom"], c(99.00144, 96.18348, 94.17755), tol = 1e-1)
expect_equal(eUN.lmm_anova[eUN.lmm_anova$type=="k","df.denom"], c(189.236), tol = 1e-1)
expect_equal(eUN.lmm_anova[eUN.lmm_anova$type=="rho","df.denom"], c(20.66968), tol = 1e-1)

## ** getVarCov
sigma(eUN.lmm)
})

## * Stratified compound symmetry structure
test_that("Stratified compound symmetry structure (REML)",{

## ** fit
eSCS0.lmm <- lmm(Y ~ gender, repetition = gender~visit|id, structure = "CS", data = dL, trace = 0, method = "REML") ## just to check that it runs
eSCS.lmm <- lmm(Y ~ (visit + age)*gender, repetition = gender~visit|id, structure = "CS", data = dL, trace = 0, method = "REML")
eSCS.gls <- list(male=gls(Y ~ visit + age, correlation = corCompSymm(form=~1|id), data = dL[dL$gender=="male",], method = "REML"),
                female=gls(Y ~ visit + age, correlation = corCompSymm(form=~1|id), data = dL[dL$gender=="female",], method = "REML"))

## ** coef
coef(eSCS.lmm, transform.rho = "cov")
coef(eSCS.lmm, transform.sigma = "square")

## ** logLikelihood
expect_equal(logLik(eSCS.lmm), as.double(logLik(eSCS.gls$male)+logLik(eSCS.gls$female)), tol = 1e-6)

## ** score
test <- score(eSCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")
GS <- jacobian(func = function(p){logLik(eSCS.lmm, p = p)},
               x = coef(eSCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"))
expect_true(all(abs(test) < 1e-3))
expect_equal(as.double(GS), as.double(test), tol = 1e-5)

## no transformation
newp <- coef(eSCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")+0.1
GS <- jacobian(func = function(p){logLik(eSCS.lmm, p = p)}, x = newp)
test <- score(eSCS.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log[c("sigma:male","sigma:female")] <- log(newp[c("sigma:male","sigma:female")])
GS <- jacobian(func = function(p){p[c("sigma:male","sigma:female")] <- exp(p[c("sigma:male","sigma:female")]); logLik(eSCS.lmm, p = p)}, x = newp.log)
test <- score(eSCS.lmm, effects = "all", p = newp, transform.sigma = "log", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2[c("sigma:male","sigma:female")] <- newp[c("sigma:male","sigma:female")]^2
GS <- jacobian(func = function(p){p[c("sigma:male","sigma:female")] <- sqrt(p[c("sigma:male","sigma:female")]); logLik(eSCS.lmm, p = p)}, x = newp.2)
test <- score(eSCS.lmm, effects = "all", p = newp, transform.sigma = "square", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## cov transformation
newp.cov <- newp.2; newp.cov[c("rho:male","rho:female")] <- c(newp["rho:male"]*newp["sigma:male"]^2, newp["rho:female"]*newp["sigma:female"]^2)
GS <- jacobian(func = function(p){
    p[c("sigma:male","rho:male","sigma:female","rho:female")] <- c(sqrt(p["sigma:male"]),p["rho:male"]/p["sigma:male"],sqrt(p["sigma:female"]),p["rho:female"]/p["sigma:female"])
    logLik(eSCS.lmm, p = p)
}, x = newp.cov)
test <- score(eSCS.lmm, p = newp, effects = "all", transform.rho = "cov")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** information
test <- information(eSCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")
GS <- -hessian(func = function(p){logLik(eSCS.lmm, p = p)},
               x = coef(eSCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"))
expect_equal(as.double(GS), as.double(test), tol = 1e-5)

## no transformation
newp <- coef(eSCS.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")+0.1
GS <- -hessian(func = function(p){logLik(eSCS.lmm, p = p)}, x = newp)
test <- information(eSCS.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log[c("sigma:male","sigma:female")] <- log(newp[c("sigma:male","sigma:female")])
GS <- -hessian(func = function(p){p[c("sigma:male","sigma:female")] <- exp(p[c("sigma:male","sigma:female")]); logLik(eSCS.lmm, p = p)}, x = newp.log)
test <- information(eSCS.lmm, effects = "all", p = newp, transform.sigma = "log", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2[c("sigma:male","sigma:female")] <- newp[c("sigma:male","sigma:female")]^2
GS <- -hessian(func = function(p){p[c("sigma:male","sigma:female")] <- sqrt(p[c("sigma:male","sigma:female")]); logLik(eSCS.lmm, p = p)}, x = newp.2)
test <- information(eSCS.lmm, effects = "all", p = newp, transform.sigma = "square", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## cov transformation
newp.cov <- newp.2; newp.cov[c("rho:male","rho:female")] <- c(newp["rho:male"]*newp["sigma:male"]^2, newp["rho:female"]*newp["sigma:female"]^2)
GS <- -hessian(func = function(p){
    p[c("sigma:male","rho:male","sigma:female","rho:female")] <- c(sqrt(p["sigma:male"]),p["rho:male"]/p["sigma:male"],sqrt(p["sigma:female"]),p["rho:female"]/p["sigma:female"])
    logLik(eSCS.lmm, p = p)
}, x = newp.cov)
test <- information(eSCS.lmm, p = newp, effects = "all", transform.rho = "cov")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** degree of freedom
test <- model.tables(eSCS.lmm, effects = "all", transform.sigma = "none")[,"df",drop=FALSE]

## ** anova
eSCS.lmm_anova <- anova(eSCS.lmm, effects = "all")$multivariate
if(eSCS.lmm$opt$name=="gls"){
    expect_equal(eSCS.lmm_anova[eSCS.lmm_anova$type=="mu","df.denom"], c(143.04096527, 46.63249177), tol = 1e-1)
}else{
    expect_equal(eSCS.lmm_anova[eSCS.lmm_anova$type=="mu","df.denom"], c(171.0218177, 56.03889211, 81.14993152, 285.63027224, 63.97167502), tol = 1e-1)
}
expect_equal(eSCS.lmm_anova[eSCS.lmm_anova$type=="rho","df.denom"], c(7.450136), tol = 1e-1)

## ** getVarCov
sigma(eSCS.lmm)
})

## * Stratified unstructed covariance matrix
test_that("Stratified unstructured (REML)",{

## ** fit
eSUN.lmm <- lmm(Y ~ (visit + age)*gender, repetition = gender~visit|id, structure = "UN", data = dL, trace = 0, method = "REML")
eSUN.gls <- list(male=gls(Y ~ visit + age, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|visit), data = dL[dL$gender=="male",], method = "REML"),
                female=gls(Y ~ visit + age, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|visit), data = dL[dL$gender=="female",], method = "REML"))
coef(eSUN.lmm, transform.sigma = "none", transform.k = "sd", effects = "variance")
coef(eSUN.lmm, transform.sigma = "none", effects = "variance")
capture.output(summary(eSUN.lmm))
## anova(eSUN.lmm, effects = "all")

## ** coef
coef(eSUN.lmm, transform.rho = "cov")
coef(eSUN.lmm, transform.k = "var")
coef(eSUN.lmm, transform.k = "log")

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
expect_equal(fct.trans(coef(eSUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")),
             coef(eSUN.lmm, effects = "all", transform.k = "var", transform.rho = "none", transform.names=FALSE))

## ** logLikelihood
expect_equal(logLik(eSUN.lmm), as.double(logLik(eSUN.gls$male)+logLik(eSUN.gls$female)), tol = 1e-6)

## ** score
test <- score(eSUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")
GS <- jacobian(func = function(p){logLik(eSUN.lmm, p = p)}, x = coef(eSUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"))
expect_true(all(abs(test) < 1e-3))
expect_equal(as.double(GS), as.double(test), tol = 1e-5)

## no transformation
newp <- coef(eSUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")+0.1
GS <- jacobian(func = function(p){logLik(eSUN.lmm, p = p)}, x = newp)
test <- score(eSUN.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log[eSUN.lmm$design$param$type %in% "sigma"] <- log(newp[eSUN.lmm$design$param$type %in% "sigma"])
GS <- jacobian(func = function(p){p[eSUN.lmm$design$param$type %in% "sigma"] <- exp(p[eSUN.lmm$design$param$type %in% "sigma"]); logLik(eSUN.lmm, p = p)}, x = newp.log)
test <- score(eSUN.lmm, effects = "all", p = newp, transform.sigma = "log", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- fct.trans(newp); 
GS <- jacobian(func = function(p){logLik(eSUN.lmm, p = fct.trans(p,inverse = TRUE))}, x = newp.2)
test <- score(eSUN.lmm, effects = "all", p = newp, transform.k = "var", transform.rho = "none", transform.names = FALSE)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## cov transformation
## newp.cov <- newp; newp.cov[c("sigma","rho")] <- c(newp["sigma"]^2,newp["rho"]*newp["sigma"]^2)
## GS <- jacobian(func = function(p){p[c("sigma","rho")] <- c(sqrt(p["sigma"]),p["rho"]/p["sigma"]); logLik(eSUN.lmm, p = p)}, x = newp.cov)
## test <- score(eSUN.lmm, p = newp, transform.rho = "cov")
## expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** information
## expect_equal(as.double(vcov(eSUN.lmm, effects = "mean")), as.double(Matrix::bdiag(lapply(eSUN.gls,vcov))), tol = 1e-6)

test <- information(eSUN.lmm, effects = "all",
                    p = coef(eSUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"),
                    transform.sigma = "none", transform.k = "none", transform.rho = "none")
GS <- -hessian(func = function(p){logLik(eSUN.lmm, p = p)},
               x = coef(eSUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"))
expect_equal(as.double(test),as.double(GS), tol = 1e-5)

test <- information(eSUN.lmm, effects = "all",
                    p = coef(eSUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"),
                    transform.sigma = "log", transform.k = "none", transform.rho = "none")
GS <- -hessian(func = function(p){p[eSUN.lmm$design$param$type %in% "sigma"]<-exp(p[eSUN.lmm$design$param$type %in% "sigma"]);logLik(eSUN.lmm, p = p)},
               x = coef(eSUN.lmm, effects = "all", transform.sigma = "log", transform.k = "none", transform.rho = "none", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(eSUN.lmm, effects = "all",
                    p = coef(eSUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none"),
                    transform.k = "var", transform.rho = "none", transform.names = FALSE) 
GS <- -hessian(func = function(p){logLik(eSUN.lmm, p = fct.trans(p,inverse = TRUE))},
               x = coef(eSUN.lmm, effects = "all", transform.k = "var", transform.rho = "none", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## test <- information(eSUN.lmm, p = coef(eSUN.lmm, transform.sigma = "none"), transform.rho = "cov") 
## GS <- -hessian(func = function(p){p[c("sigma","rho")] <- c(sqrt(p["sigma"]),p["rho"]/p["sigma"]); logLik(eSUN.lmm, p = p)}, x = coef(eSUN.lmm, transform.rho = "cov", transform.names = FALSE))
## expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## no transformation 
newp <- coef(eSUN.lmm, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")+0.1
GS <- -hessian(func = function(p){logLik(eSUN.lmm, p = p)}, x = newp)
test <- information(eSUN.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none", transform.rho = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6) 

## ** degree of freedom
test <- confint(eSUN.lmm, effects = "all")[,"df", drop=FALSE]
## anova(eSUN.lmm)
## anova(eSUN.gls)

## ** anova
eSUN.lmm_anova <- anova(eSUN.lmm, effects = "all", ci = TRUE)$multivariate
if(eSUN.lmm$opt$name=="gls"){
    expect_equal(eSUN.lmm_anova[eSUN.lmm_anova$type=="mu","df.denom"], c(47.63744, 42.53266), tol = 1e-1)
}else{
    expect_equal(eSUN.lmm_anova[eSUN.lmm_anova$type=="mu","df.denom"], c(57.00206139, 55.02679117, 67.96985137, 93.51660476, 54.92384757), tol = 1e-1)
}
expect_equal(eSUN.lmm_anova[eSUN.lmm_anova$type=="k","df.denom"], c(86.07826), tol = 1e-1)
expect_equal(eSUN.lmm_anova[eSUN.lmm_anova$type=="rho","df.denom"], c(9.100919), tol = 1e-1)


## ** getVarCov
sigma(eSUN.lmm)
})

## * Random intercept model
data(Orthodont,package="nlme")

test_that("Random intercept model",{
    eRI.lmer <- lmer(distance ~ age + (1|Subject), data=Orthodont)
    eRI.lmm <- lmm(distance ~ age + (1|Subject), data=Orthodont,
                   control = list(init = "lmer"))
    summary(eRI.lmm)
    ## likelihood
    expect_equal(as.double(logLik(eRI.lmer)), as.double(logLik(eRI.lmm)), tol = 1e-6)

    ## random effects (conditional mean)
    u.GS <- as.data.frame(ranef(eRI.lmer))    
    u.test <- ranef(eRI.lmm, effects = "mean", format = "long")
    expect_equal(as.double(u.GS$condval[match(u.test$level,u.GS$grp)]), as.double(u.test$estimate), tol = 1e-6)

    ## random effects (conditional variance)
    tau.GS <- as.data.frame(VarCorr(eRI.lmer))
    tau.test <- ranef(eRI.lmm, effects = "variance", format = "long")
    expect_equal(as.double(tau.GS[1,"vcov"]), as.double(tau.test[tau.test$type=="variance","estimate"]), tol = 1e-6)
})

## * Stratified random intercept model
Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")

test_that("Random intercept model",{
    eRI.mlmm <- mlmm(distance ~ 1 + (1|Subject), data = Orthodont, by = "nsex", trace = FALSE)
    eSRI.lmm <- lmm(distance ~ nsex + (1|Subject), repetition = nsex~1|Subject, data = Orthodont)
    
    ## likelihood
    expect_equal(as.double(sum(unlist(logLik(eRI.mlmm)))), as.double(logLik(eSRI.lmm)), tol = 1e-6)

    ## random effects (conditional mean)
    u.GS <- unlist(unname(unlist(unname(ranef(eRI.mlmm)), recursive = FALSE)), recursive = FALSE)
    u.test <- ranef(eSRI.lmm, effects = "mean", format = "long")
    expect_equal(as.double(u.GS[names(test)]), as.double(test), tol = 1e-6)

    ## random effects (conditional variance)
    tua.GS <- ranef(eRI.mlmm, effects = "variance")
    tau.test <- ranef(eSRI.lmm, effects = "variance", format = "long")
    expect_equal(as.double(GS), as.double(test), tol = 1e-5)
})

## * Crossed random intercept model (2 terms)
data(Penicillin, package = "lme4")
Penicillin$id <- 1


test_that("Crossed random intercept model (2 terms)",{

    eCRI2.lmer <- lmer(diameter ~ (1|plate) + (1|sample), data = Penicillin)
    ## eCRI2.lmm0 <- lmm(diameter ~ 1, repetition = ~1|id,
    ##                   structure = CS(list(~1,~plate+sample), type = "ho", group = 1:2),
    ##                   data = Penicillin, df = FALSE)
    ## eCRI2.lmm <- lmm(diameter ~ (1|plate) + (1|sample), data = Penicillin, df = FALSE,
    ##                  control = list(init = "lmer"), trace = 5)
    eCRI2.lmm <- lmm(diameter ~ (1|plate) + (1|sample), data = Penicillin, df = FALSE)

    ## likelihood
    expect_equal(as.double(logLik(eCRI2.lmer)), as.double(logLik(eCRI2.lmm)), tol = 1e-6)

    ## random effects (conditional mean)
    GS <- do.call(rbind,ranef(eCRI2.lmer))
    test <- do.call(c,ranef(eCRI2.lmm))
    expect_equal(as.double(GS[,1]), as.double(test), tol = 1e-6)

    ## random effects (conditional variance)
    GS <- as.data.frame(VarCorr(eCRI2.lmer))
    test <- ranef(eCRI2.lmm, effects = "variance")
    expect_equal(as.double(test), as.double(GS[1:2,"vcov"]), tol = 1e-3)
})

## * Crossed random intercept model (3 terms)
Sigma.CRI3 <- matrix(0, nrow = 6, ncol = 6)
Sigma.CRI3[1:3,1:3] <- Sigma.CRI3[1:3,1:3] <- 0.6^2
diag(Sigma.CRI3[1:3,4:6]) <- diag(Sigma.CRI3[4:6,1:3]) <- 0.4^2
Sigma.CRI3[1,6] <- Sigma.CRI3[2,4] <- Sigma.CRI3[3,5] <- Sigma.CRI3[6,1] <- Sigma.CRI3[4,2] <- Sigma.CRI3[5,3] <- 0.2^2
diag(Sigma.CRI3) <- 1

n <- 100
dfL.CRI3 <- reshape2::melt(data.frame(sample = 1:n,
                                      mvtnorm::rmvnorm(n, mean = 1:6, sigma = Sigma.CRI3)
                                      ),
                           id.vars = "sample")
dfL.CRI3$patient <- paste(dfL.CRI3$sample,dfL.CRI3$variable %in% paste0("X",4:6)+1, sep = ".")
dfL.CRI3$day <- paste(dfL.CRI3$sample, sapply(as.character(dfL.CRI3$variable), switch, "X1" = 1, "X2" = 2, "X3" = 3, "X4" = 1, "X5" = 2, "X6" = 3),sep=".")
dfL.CRI3$batch <- paste(dfL.CRI3$sample, sapply(as.character(dfL.CRI3$variable), switch, "X1" = 1, "X2" = 2, "X3" = 3, "X4" = 2, "X5" = 3, "X6" = 1),sep=".")
dfL.CRI3 <- dfL.CRI3[order(dfL.CRI3$sample),]

## head(dfL.CRI3,10)

test_that("Crossed random intercept model (3 terms)",{

    eCRI3.lmer <- lmer(value ~ 0 + variable + (1|patient) + (1|day) + (1|batch), data = dfL.CRI3)
    eCRI3.lmm <- lmm(value ~ 0 + variable + (1|patient) + (1|day) + (1|batch), data = dfL.CRI3, df = FALSE,
                     trace = 5)

    ## likelihood
    expect_equal(as.double(logLik(eCRI3.lmer)), as.double(logLik(eCRI3.lmm)), tol = 1e-6)

    ## random effects (conditional mean)
    GS <- do.call(rbind,ranef(eCRI3.lmer))
    test <- do.call(c,ranef(eCRI3.lmm))
    expect_equal(as.double(GS[,1]), as.double(test), tol = 1e-6)

    ## random effects (conditional variance)
    GS <- as.data.frame(VarCorr(eCRI3.lmer))
    test <- ranef(eCRI3.lmm, effects = "variance")
    expect_equal(as.double(test), as.double(GS[1:2,"vcov"]), tol = 1e-3)
})

## as.data.frame(ranef(eCRI3.lmer))


set.seed(10)
df.grow <- expand.grid(plant = as.factor(letters[1:5]),
                       location = as.factor(LETTERS[1:5]),
                       climate = 1:5)
u1 <- rnorm(length(unique(df.grow$plant)), sd = 0.5)
u2 <- rnorm(length(unique(df.grow$location)), sd = 1)
u3 <- rnorm(length(unique(df.grow$climate)), sd = 2)
df.grow$y <- rnorm(NROW(df.grow)) + u1[as.numeric(df.grow$plant)] + u2[as.numeric(df.grow$location)] + u3[df.grow$climate]
df.grow$id <- 1

eCRI3.lmm.bis <- lmm(y ~ 1, repetition = ~1|id, structure = CS(list(~1,~plant+location+climate), heterogeneous = -1), data = df.grow,
                     trace = 5)
coef(eCRI3.lmm.bis, effects = "all")

eCRI3.lmer <- lmer(y ~ (1|plant) + (1|location) + (1|climate), data = df.grow)
eCRI3.lmm <- lmm(y ~ (1|plant) + (1|location) + (1|climate), data = df.grow,
                 trace = 5)


summary(eCRI3.lmer)

## * Nested random intercept model (2 levels)
## https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified
## df.school <- read.table("http://bayes.acs.unt.edu:8083/BayesContent/class/Jon/R_SC/Module9/lmm.data.txt",
##                  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
## eNRI2.lmer <- lmer(extro ~ open + agree + social + (1 | school/class), data = dt)
## dt.red <- as.data.table(dt)[,.SD[1:10,.(id,extro=round(extro,2))],by = c("class","school")]
## setkeyv(dt.red, c("school","class"))
## dt.red$id <- 1:NROW(dt.red)

dt.red <- data.table("class" = rep(unlist(lapply(c("a","b","c","d"),rep,10)), 6), 
           "school" = unlist(lapply(c("I","II","III","IV","V","VI"),rep,40)), 
           "id" = 1:240, 
           "extro" = c(41.10, 43.24, 42.15, 43.72, 40.59, 41.38, 35.43, 39.46, 41.06, 38.96, 47.24, 45.77, 44.32, 47.34, 45.06, 47.01, 46.07, 44.41, 45.30, 45.20, 47.65, 47.76, 48.59, 47.46, 48.63, 48.31, 49.28, 47.81, 49.01, 48.00, 50.97, 49.28, 49.66, 50.32, 49.33, 49.87, 50.87, 51.43, 50.95, 50.87, 52.02, 51.67, 52.62, 51.96, 51.88, 52.56, 51.54, 52.08, 52.76, 52.69, 53.69, 54.10, 54.07, 53.70, 53.03, 53.88, 53.80, 53.50, 53.42, 53.17, 55.17, 54.58, 55.06, 54.40, 54.86, 54.45, 55.27, 54.60, 55.33, 54.81, 55.63, 55.65, 55.96, 56.01, 55.58, 56.01, 56.08, 55.47, 56.07, 55.56, 56.68, 57.30, 56.84, 57.24, 56.57, 56.94, 56.81, 57.11, 57.04, 56.35, 57.75, 57.84, 57.66, 58.03, 58.22, 57.50, 57.54, 57.42, 58.11, 58.08, 59.03, 58.41, 59.04, 59.09, 58.94, 58.39, 58.30, 58.80, 58.95, 58.35, 60.15, 59.45, 59.43, 59.98, 59.40, 60.05, 60.00, 59.90, 59.93, 59.80, 60.65, 60.84, 60.66, 60.31, 60.98, 60.75, 60.60, 60.42, 60.74, 60.17, 62.29, 61.84, 61.57, 61.17, 62.12, 61.59, 61.42, 61.36, 61.67, 62.12, 62.97, 63.15, 62.82, 62.84, 63.29, 62.68, 62.45, 63.21, 62.80, 62.75, 63.69, 64.25, 64.18, 64.14, 63.72, 63.42, 64.18, 63.68, 63.48, 64.04, 64.62, 65.13, 65.40, 65.23, 64.56, 64.63, 64.81, 64.31, 64.76, 64.44, 65.92, 66.49, 66.28, 65.68, 66.03, 66.46, 65.62, 66.01, 65.61, 66.10, 66.57, 66.98, 66.55, 66.96, 66.64, 67.24, 67.23, 67.34, 66.94, 67.30, 67.94, 68.01, 68.73, 68.61, 68.89, 68.91, 68.84, 68.06, 68.13, 68.93, 69.48, 70.12, 70.17, 70.67, 70.60, 69.62, 70.16, 70.65, 70.34, 69.87, 71.79, 71.36, 72.45, 71.07, 70.87, 70.94, 71.25, 72.48, 72.04, 71.67, 74.35, 74.70, 73.14, 74.44, 73.98, 75.03, 74.21, 76.51, 75.94, 73.42, 79.74, 78.15, 78.19, 83.34, 80.11, 79.73, 80.64, 77.07, 78.02, 82.25))


eNRI2.lmer <- lmer(extro ~ (1|school/class), data = dt.red)
ranef(eNRI2.lmer)



test_that("Nested random intercept model (2 levels)",{

    eNRI2.lmm <- lmm(angle ~ recipe * temperature + (1|recipe:replicate), data = cake, df = FALSE)

    ## likelihood
    expect_equal(as.double(logLik(eNRI2.lmer)), as.double(logLik(eNRI2.lmm)), tol = 1e-6)

    ## random effects
    GS <- as.data.frame(ranef(eNRI2.lmer))    
})


## * Nested random intercept model (3 levels)
Sigma.NRI3 <- matrix(0, nrow = 8, ncol = 8)
Sigma.NRI3[5:8,1:4] <- Sigma.NRI3[1:4,5:8] <- 0.25^2
Sigma.NRI3[1:2,3:4] <- Sigma.NRI3[3:4,1:2] <- Sigma.NRI3[5:6,7:8] <- Sigma.NRI3[7:8,5:6] <- 0.5^2
Sigma.NRI3[1:2,1:2] <- Sigma.NRI3[3:4,3:4] <- Sigma.NRI3[5:6,5:6] <- Sigma.NRI3[7:8,7:8] <- 0.8^2
diag(Sigma.NRI3) <- 1

## c(sqrt(0.25^2), sqrt(0.5^2-0.25^2), sqrt(0.8^2-0.5^2-0.25^2))
n <- 1000
dfL.NRI3 <- reshape2::melt(data.frame(patient = 1:n,
                                      mvtnorm::rmvnorm(n, mean = 1:8, sigma = Sigma.NRI3)
                                      ),
                           id.vars = "patient")
dfL.NRI3$day <- dfL.NRI3$variable %in% paste0("X",5:8)+1
dfL.NRI3$session <- paste(dfL.NRI3$day,dfL.NRI3$variable %in% c(paste0("X",c(3:4,7:8)))+1,sep=".")
dfL.NRI3 <- dfL.NRI3[order(dfL.NRI3$patient),]

## head(dfL.NRI3,10)

eNRI3.lmer <- lmer(value ~ session + (1|patient/day/session), data = dfL.NRI3)
## as.data.frame(ranef(eNRI3.lmer))

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
})

## * Baseline constraint
test_that("Baseline constraint",{

    dL$group <- as.factor(dL$id %% 2)
    dL$treat <- (dL$group==1)*(dL$visit!="Y1")
    table(dL$treat, baselineAdjustment(dL, variable = "group", repetition = ~visit|id, constrain = "Y1", new.level = "0"))
    dL$treat.visit <- baselineAdjustment(dL, variable = "group", repetition = ~visit|id, constrain = "Y1", collapse.time = ".")

    ## eUN.lmm <- lmm(Y ~ group*visit, repetition = ~group*visit|id, structure = "UN", data = dL, trace = 0, method = "REML")
    ## logLik(eUN.lmm)
    eCUN.lmm <- suppressMessages(lmm(Y ~ treat*visit, repetition = ~treat*visit|id, structure = "UN", data = dL, trace = 0, method = "REML", df = FALSE, control = list(optimizer = "FS")))
    eCUN2.lmm <- lmm(Y ~ treat.visit, repetition = ~treat.visit|id, structure = "UN", data = dL, trace = 0, method = "REML", df = FALSE, control = list(optimizer = "FS"))
    
    expect_equal(logLik(eCUN2.lmm), logLik(eCUN.lmm), tol = 1e-5)
    expect_equal(logLik(eCUN2.lmm), -618.14359397, tol = 1e-5)

    capture.output(summary(eCUN2.lmm))
    capture.output(summary(anova(eCUN2.lmm), method = "none"))
    ## plot(eCUN2.lmm, color = "group", time.var = "visit")

    ## baseline constrain for order 3 interaction
    eCUN.I2.lmm <- suppressMessages(lmm(Y ~ gender*treat*visit, repetition = ~treat*visit|id, structure = "UN", data = dL, trace = 0, method = "REML", df = FALSE, control = list(optimizer = "FS")))
    eCUN2.I2.lmm <- suppressMessages(lmm(Y ~ gender:treat.visit, repetition = ~treat.visit|id, structure = "UN", data = dL, trace = 0, method = "REML", df = FALSE, control = list(optimizer = "FS")))
    expect_equal(logLik(eCUN.I2.lmm), logLik(eCUN2.I2.lmm), tol = 1e-5)
    expect_equal(logLik(eCUN.I2.lmm), -598.96051963, tol = 1e-5)

})
##----------------------------------------------------------------------
### test-auto-mixed-model.R ends here
