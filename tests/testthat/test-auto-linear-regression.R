### test-auto-linear-regression.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (10:13) 
## Version: 
## Last-Updated: nov 13 2021 (18:18) 
##           By: Brice Ozenne
##     Update #: 188
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
    library(nlme)
    library(AICcmodavg)

    library(LMMstar)
}

context("Check lmm on examples of linear regression")
LMMstar.options(optimizer = "gls", method.numDeriv = "Richardson", precompute.moments = TRUE,
                columns.confint = c("estimate","se","df","lower","upper","p.value"))

## * simulate data
n <- 5e1
p <- 3
X.name <- paste0("X",1:p)
link.lvm <- paste0("Y~",X.name)
formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))

m <- lvm(formula.lvm)
distribution(m,~Id) <- Sequence.lvm(a = 1, b = n)
distribution(m,~Gender.num) <- binomial.lvm()
distribution(m,~X3) <- binomial.lvm(size=2)
transform(m,Gene~X3) <- function(x){factor(x,levels=0:2,labels=c("LL","LA","AA"))}
transform(m,Gender~Gender.num) <- function(x){factor(x,levels=0:1,labels=c("M","F"))}
transform(m,id~Gender.num) <- function(x){paste0("id",1:NROW(x))}
latent(m) <- ~Gender.num+X3
set.seed(10)
d <- lava::sim(m,n, latent = FALSE)
d$time <- "t1"


## * single variance parameter (ML)

test_that("single variance parameter (ML)",{

## ** fit
e.lmm <- lmm(Y ~ X1 + X2 + Gene, repetition = ~time|id, data = d, trace = 0, method = "ML")
e.gls <- gls(Y ~ X1 + X2 + Gene, data = d, method = "ML")
e.lava <- estimate(lvm(Y~X1+X2+Gene),data = d)

n.obs <- unname(nobs(e.lmm)[1])
n.mu <- length(coef(e.lmm, effects = "mean"))
n.sigma <- length(coef(e.lmm, effects = "variance"))
n.param <- length(coef(e.lmm, effects = "all"))

## ** coef
expect_equal(coef(e.lmm, effects = "mean"), coef(e.gls), tol = 1e-6)
expect_equal(unname(coef(e.lmm, transform.sigma = "square")), unname(coef(e.lava)), tol = 1e-6)

## ** logLikelihood
expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)
expect_equal(logLik(e.lmm), as.double(logLik(e.lava)), tol = 1e-6)
## coef(e.lmm, transform.sigma = "log")
## coef(e.lmm, transform.sigma = "square")
## coef(e.lmm, transform.sigma = "logsquare")

## no transformation
newp <- coef(e.lmm, effects = "all", transform.sigma = "none")+1
newp.lava <- coef(e.lava) + 1 ; newp.lava["Y~~Y"] <- (sqrt(newp.lava["Y~~Y"]-1)+1)^2
expect_equal(logLik(e.lmm, p = newp), as.double(logLik(e.lava, p = newp.lava)), tol = 1e-6)

## ** score
expect_true(all(abs(lava::score(e.lmm)) < 1e-6))
## score(e.lmm, transform.sigma = "log")
## score(e.lmm, transform.sigma = "square")
## score(e.lmm, transform.sigma = "logsquare")

## no transformation
newp <- coef(e.lmm, effects = "all", transform.sigma = "none")+1
GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p)}, x = newp.log)
test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "log")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)


## lava transformation
newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2
GS <- jacobian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(e.lmm, p = p)}, x = newp.2)
GS0 <- score(e.lava, p = newp.2)
test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "square")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)
expect_equal(as.double(GS0), as.double(GS), tol = 1e-6)

## ** information
test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "none", type.information = "observed")
testE <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "none", type.information = "expected")
testE2 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "square", type.information = "expected")
expect_equal(as.double(test),as.double(testE), tol = 1e-6)
expect_equal(unname(testE["sigma","sigma"]),unname(2*n.obs/coef(e.lmm, effects = "all", transform.sigma = "none")["sigma"]^2), tol = 1e-6)
expect_equal(unname(testE2["sigma^2","sigma^2"]),unname(n.obs/(2*coef(e.lmm, effects = "all", transform.sigma = "none")["sigma"]^4)), tol = 1e-6)

GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none"))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "log")
GS <- -hessian(func = function(p){p["sigma"]<-exp(p["sigma"]);logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "log", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "square") 
GS <- -hessian(func = function(p){p["sigma"]<-sqrt(p["sigma"]);logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "square", transform.names = FALSE))
GS0 <- -hessian(func = function(p){logLik(e.lava, p = p)}, x = coef(e.lava))
expect_equal(as.double(test["sigma^2","sigma^2"]), as.double(nobs(e.lmm)[1]/(2*coef(e.lmm, effect = "variance", transform.sigma = "square")^2)), tol = 1e-6)
expect_equal(as.double(test), as.double(GS), tol = 1e-6)
expect_equal(as.double(GS0), as.double(GS), tol = 1e-6)

## no transformation 
newp <- coef(e.lmm, effects = "all", transform.sigma = "none")+1
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- information(e.lmm, effects = "all", p = newp, transform.sigma = "none")
expect_equal(as.double(test[1:4,1:4]), as.double(GS[1:4,1:4]), tol = 1e-6) ## does not match as some terms do not cancel

## ** variance-covariance
test <- vcov(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "none")
GS <- solve(-hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none")))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test <- vcov(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "square", df = TRUE)
GS <- solve(-hessian(func = function(p){p["sigma"] <- sqrt(p["sigma"]);logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "square", transform.names = FALSE)))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

expect_equal(unname(sqrt(diag(vcov(e.lmm, effects = "all", transform.sigma="square",robust=TRUE)))), unname(estimate(e.lava, robust = TRUE)$coefmat[,"Std.Err"]),
             tol = 1e-4)

## ** degree of freedom
test <- confint(e.lmm, effects = "all", transform.sigma = "log", type.information = "observed")$df
expect_equal(test, rep(n.obs,n.param), tol = 1e-6)

test <- suppressWarnings(confint(e.lmm, effects = "all", transform.sigma = "square", type.information = "expected")$df)
expect_equal(test, c(rep(n.obs,n.mu),n.obs/4), tol = 1e-6)

## ## numerical derivative with appropriate transformation
## FF <- function(p){p["sigma"] <- sqrt(p["sigma"]);diag(vcov(e.lmm, p = p, transform.sigma = "square"))}
## GG <- jacobian(func = FF, x = coef(e.lmm, transform.sigma = "square", transform.names = FALSE))
## VV <- vcov(e.lmm, p = coef(e.lmm), transform.sigma = "square")
## test2 <- sapply(1:NROW(GG),function(gg){2*VV[gg,gg]^2/ (GG[gg,,drop=FALSE] %*% VV %*% t(GG[gg,,drop=FALSE]))})
## expect_equal(unname(test2), unname(GS))

## ## numerical derivative on another scale
## FF.bis <- function(p){diag(vcov(e.lmm, p = p, transform.sigma = "none"))}
## GG.bis <- jacobian(func = FF.bis, x = coef(e.lmm, transform.sigma = "none"))
## VV.bis <- vcov(e.lmm, p = coef(e.lmm), transform.sigma = "none")
## test2 <- sapply(1:NROW(GG.bis),function(gg){2*VV.bis[gg,gg]^2/ (GG.bis[gg,,drop=FALSE] %*% VV.bis %*% t(GG.bis[gg,,drop=FALSE]))})

## ** anova
test <- anova(e.lmm)
expect_equal(test$mean["Gene","statistic"]*test$mean["Gene","df.num"],
             unname(lava::compare(e.lava, par = c("Y~GeneLA","Y~GeneAA"))$statistic), tol = 1e-6)
expect_equal(test$mean["Gene","df.denom"], NROW(d))

test <- anova(e.lmm, effect = c("GeneLA=0","GeneAA=0"))
test2 <-  anova(e.lmm, effect = c("GeneLA=0","GeneAA=0"))

})

## * single variance parameter (REML)
test_that("single variance parameter (REML)",{

## ** fit
e.lmm <- lmm(Y ~ X1 + X2 + Gene, repetition = ~time|id, data = d, trace = 0, method = "REML", df = TRUE)
e.gls <- gls(Y ~ X1 + X2 + Gene, data = d, method = "REML")
e.lm <- lm(Y ~ X1 + X2 + Gene, data = d)

n.obs <- unname(nobs(e.lmm)[1])
n.mu <- length(coef(e.lmm, effects = "mean"))
n.sigma <- length(coef(e.lmm, effects = "variance"))
n.param <- length(coef(e.lmm, effects = "all"))

## ** coef
expect_equal(coef(e.lmm, effects = "mean"), coef(e.gls), tol = 1e-6)

## ** logLikelihood
expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)

## ** score
expect_true(all(abs(score(e.lmm, effects = "all")) < 1e-6))

## no transformation
newp <- coef(e.lmm, effects = "all", transform.sigma = "none")+1
GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## transformation
newp.log <- newp; newp.log["sigma"] <- log(newp.log["sigma"])
GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p)}, x = newp.log)
test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "log")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** information
test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "none", type.information = "observed")
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none"))

expect_equal(as.double(test), as.double(GS), tol = 1e-6)
expect_equal(unname(test["sigma","sigma"]),unname(2*(n.obs-n.mu)/coef(e.lmm, effects = "all", transform.sigma = "none")["sigma"]^2), tol = 1e-6)

testE <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "none", type.information = "expected")
## difference: in one case (information)  : n / sigma^2
##             in the other case (hessian): n / sigma^2 - 3 * sum(residuals^2)/sigma^4
## which are equal when sum(residuals^2)/n = sigma^2

test2 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "square", type.information = "observed")
expect_equal(unname(test2["sigma^2","sigma^2"]),unname((n.obs-n.mu)/(2*coef(e.lmm, effects = "all", transform.sigma = "none")["sigma"]^4)), tol = 1e-6)

test1 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "log", type = "observed")
test2 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "log", type = "expected")
GS <- -hessian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p)}, x = coef(e.lmm, transform.sigma = "log", transform.names = FALSE))
expect_equal(as.double(test1), as.double(GS), tol = 1e-6)

test1 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "square", type = "observed")
test2 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "square", type = "expected")
GS <- -hessian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(e.lmm, p = p)}, x = coef(e.lmm, transform.sigma = "square", transform.names = FALSE))
expect_equal(as.double(test1), as.double(GS), tol = 1e-6)
expect_equal(as.double(test1["sigma^2","sigma^2"]), as.double((nobs(e.lmm)[1]-length(coef(e.lmm, effect = "mean")))/(2*coef(e.lmm, effect = "variance", transform.sigma = "square")^2)), tol = 1e-6)

## no transformation 
newp <- coef(e.lmm, effects = "all", transform.sigma = "none")+1
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- information(e.lmm, effects = "all", p = newp, transform.sigma = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** degree of freedom
test <- confint(e.lmm, effects = "all", transform.sigma = "log", type.information = "observed")$df
expect_equal(test, rep(n.obs-n.mu,n.param), tol = 1e-6)

## numerical derivative
test <- confint(e.lmm, effects = "all", transform.sigma = "log", type.information = "observed")$df
FF.bis <- function(p){p["sigma"] <- exp(p["sigma"])   ; diag(vcov(e.lmm, effects = "all", p = p, transform.sigma = "log", type.information = "observed"))}
GG.bis <- jacobian(func = FF.bis, x = coef(e.lmm, effects = "all", transform.sigma = "log", transform.names = FALSE), method = LMMstar.options()$method.numDeriv)
VV.bis <- vcov(e.lmm, effects = "all", transform.sigma = "log", type.information = "observed")
GS <- sapply(1:NROW(GG.bis),function(gg){2*VV.bis[gg,gg]^2/ (GG.bis[gg,,drop=FALSE] %*% VV.bis %*% t(GG.bis[gg,,drop=FALSE]))})
expect_equal(unname(test), unname(GS), tol = 1e-6)

## test <- vcov(e.lmm, df = 2, transform.sigma = "log", type.information = "observed")
## all(abs(test-VV.bis)<1e-10)
## apply(attr(test,"dVcov"),3,diag)
## GG.bis

## numerical derivative on another scale
test <- confint(e.lmm, effects = "all", transform.sigma = "square", type.information = "observed")$df
FF.bis <- function(p){p["sigma"] <- sqrt(p["sigma"]); diag(vcov(e.lmm, effects = "all", p = p, transform.sigma = "square", type.information = "observed"))}
GG.bis <- jacobian(func = FF.bis, x = coef(e.lmm, effects = "all", transform.sigma = "square", transform.names = FALSE), method = LMMstar.options()$method.numDeriv)
VV.bis <- vcov(e.lmm, effects = "all", transform.sigma = "square", type.information = "observed")
GS <- sapply(1:NROW(GG.bis),function(gg){2*VV.bis[gg,gg]^2/ (GG.bis[gg,,drop=FALSE] %*% VV.bis %*% t(GG.bis[gg,,drop=FALSE]))})
expect_equal(unname(test), unname(GS), tol = 1e-6)

## ** anova
test <- anova(e.lmm)
GS <- anova(e.gls, type = "marginal")
expect_equal(test$mean$statistic,GS[["F-value"]][-1], tol = 1e-6)
expect_equal(test$mean$p.value,GS[["p-value"]][-1], tol = 1e-6)

## ** predictions
test <- predict(e.lmm, newdata = d)
GS <- predict(e.lm, newdata = d, se = TRUE)
expect_equal(test$estimate,as.double(GS$fit), tol = 1e-7)
expect_equal(test$se,as.double(GS$se.fit), tol = 1e-7)
expect_true(all(abs(test$df-GS$df) < 1e-7))

})


## * multiple variance parameters (ML)
test_that("multiple variance parameter (ML)",{

## ** fit
e.lmm <- lmm(Y ~ -1 + Gender + (X1 + X2 + Gene):Gender, repetition = ~Gender|id, structure = "UN", data = d, trace = 0,
             method = "ML")
e.gls <- gls(Y ~ -1 + Gender + (X1 + X2 + Gene):Gender, data = d, weights = varIdent(form=~1|Gender), method = "ML")
e.lmm2 <- lmm(Y ~(X1 + X2 + Gene)*Gender, repetition = Gender~time|id, structure = "CS", data = d, trace = 0,
             method = "ML")

n.obs <- unname(nobs(e.lmm)[1])
n.mu <- length(coef(e.lmm, effects = "mean"))
n.sigma <- length(coef(e.lmm, effects = "variance"))
n.param <- length(coef(e.lmm, effects = "all"))

## ** coef
expect_equal(coef(e.lmm, effects = "mean")[names(coef(e.gls))], coef(e.gls), tol = 1e-6)
## coef(e.lmm, transform = 2)
## sigma(e.gls)^2

## ** logLikelihood
expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)
## coef(e.lmm, transform.sigma = "log", transform.k = "log")
## coef(e.lmm, transform.sigma = "square", transform.k = "square")
## coef(e.lmm, transform.sigma = "logsquare", transform.k = "logsquare")
expect_equal(logLik(e.lmm2),logLik(e.lmm))

## ** score
test <- score(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none")
GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)
expect_true(all(abs(test) < 1e-4))

## no transformation
newp <- coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none")+1
GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"]); newp.log["k.F"] <- log(newp["k.F"])
GS <- jacobian(func = function(p){p[c("sigma","k.F")] <- exp(p[c("sigma","k.F")]); logLik(e.lmm, p = p)}, x = newp.log)
test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "log", transform.k = "log")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2; newp.2["k.F"] <- newp["k.F"]^2*newp["sigma"]^2
GS <- jacobian(func = function(p){p[c("sigma","k.F")] <- c(sqrt(p["sigma"]),sqrt(p["k.F"]/p["sigma"])); logLik(e.lmm, p = p, transform.sigma = "none")}, x = newp.2)
test <- score(e.lmm, effects = "all", p = newp, transform.k = "var")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## TO KEEP: DEBUG TRANSFORMATION
## .transformDeriv(transform = 2, sigma = newp["sigma"], k = newp["k.M"], rho = NULL, pattern.param = attr(e.lmm$design$X.var,"Upattern.param")[[2]])
## FCT_TRANS <- function(p){
##     c(sqrt(p["sigma"]), sqrt(p["k.M"]/p["sigma"]))
## }
## jacobian(FCT_TRANS, newp.2[c("sigma","k.M")])
## c(1/(2*sqrt(newp.2["sigma"])),-sqrt(newp.2["k.M"])/(2*newp.2["sigma"]^(3/2)), 1/(2*sqrt(newp.2["sigma"])*sqrt(newp.2["k.M"])))

## ** information
test0 <- information(e.lmm, effects = "all",
                     p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"),
                     transform.sigma = "none", transform.k = "none", type.information = "expected") ## using sigma
                     
                      test <- information(e.lmm, effects = "all",
                                         p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"),
                                         transform.sigma = "none", transform.k = "none", type.information = "observed") ## using sigma
                                         GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"))
                                         expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test0 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.sigma = "log", transform.k = "log", type.information = "expected") ## using log(sigma) and log(k)
test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.sigma = "log", transform.k = "log", type.information = "observed") ## using log(sigma) and log(k)
GS <- -hessian(func = function(p){p[c("sigma","k.F")] <- exp(p[c("sigma","k.F")]);logLik(e.lmm, p = p, transform.sigma = "none", transform.k = "none")}, x = coef(e.lmm, transform.sigma = "log", transform.k = "log", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test0 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.k = "var", type.information = "expected") ## using sigma^2 and sigma^2 k^2
test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.k = "var", type.information = "observed") ## using sigma^2 and sigma^2 k^2
GS <- -hessian(func = function(p){p[c("sigma","k.F")] <- c(sqrt(p["sigma"]),sqrt(p["k.F"]/p["sigma"]));logLik(e.lmm, p = p, transform.sigma = "none", transform.k = "none")}, x = coef(e.lmm, transform.k = "var", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** degree of freedom

test <- confint(e.lmm, effects = "all", transform.k = "logsd", type.information = "observed", backtransform = TRUE)[,"df",drop=FALSE]
name.coefM <- grep(rownames(test),pattern="M",value=TRUE)
name.coefF <- grep(rownames(test),pattern="F",value=TRUE)
expect_equal(test[name.coefM,], rep(sum(d$Gender=="M"),length(name.coefM)), tol = 1e-3)
expect_equal(test[name.coefF,], rep(sum(d$Gender=="F"),length(name.coefM)), tol = 1e-3)

test <- suppressWarnings(confint(e.lmm, effects = "all", transform.k = "var", type.information = "expected")[,"df",drop=FALSE])
name.coefM <- grep(rownames(test),pattern="M",value=TRUE)
name.coefF <- grep(rownames(test),pattern="F",value=TRUE)
expect_equal(test[name.coefM,], c(rep(sum(d$Gender=="M"),length(name.coefM)-1),sum(d$Gender=="M")/4), tol = 1e-6)
expect_equal(test[name.coefF,], c(rep(sum(d$Gender=="F"),length(name.coefM)-1),sum(d$Gender=="F")/4), tol = 1e-6)

## ** variance-covariance
getVarCov(e.lmm)
getVarCov(e.lmm2)

## ** residuals
residuals(e.lmm, format = "long")
residuals(e.lmm, format = "wide")
residuals(e.lmm2)

residuals(e.lmm, format = "wide", type = "normalized")

## ** confidence interval
expect_equal(confint(e.lmm, transform.k = "sd", effects = "variance")$estimate,
             confint(e.lmm, effects = "variance", transform.k = "logsd", backtransform = TRUE)$estimate,
             tol = 1e-6)
vcov(e.lmm, effects = "variance", df = TRUE, p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"))

## anova(e.lmm)
## anova(e.gls)

## ** print and summary
capture.output(print(e.lmm))
capture.output(summary(e.lmm))

## ** interface to other packages
glht(e.lmm)
## capture.output(emmeans(e.lmm, specs = ~Gender))

})

## * multiple variance parameters (REML)

test_that("multiple variance parameters (REML)",{

## ** fit
e.lmm <- lmm(Y ~ -1 + Gender + (X1 + X2 + Gene):Gender, repetition = ~Gender|id, structure = "UN", data = d, trace = 0,
             method = "REML")
e.lmm2 <- lmm(Y ~(X1 + X2 + Gene)*Gender, repetition = Gender~time|id, structure = "CS", data = d, trace = 0,
             method = "REML")
e.gls <- gls(Y ~ -1 + Gender + (X1 + X2 + Gene):Gender, data = d, weights = varIdent(form=~1|Gender), method = "REML")
e.gls2 <- gls(Y ~ (X1 + X2 + Gene)*Gender, data = d, weights = varIdent(form=~1|Gender), method = "REML")

n.obs <- unname(nobs(e.lmm)[1])
n.mu <- length(coef(e.lmm, effects = "mean"))
n.sigma <- length(coef(e.lmm, effects = "variance"))
n.param <- length(coef(e.lmm, effects = "all"))

## ** coef
expect_equal(coef(e.lmm, effects = "mean")[names(coef(e.gls))], coef(e.gls), tol = 1e-6)
## coef(e.lmm, transform = 2)
## sigma(e.gls)^2

## ** logLikelihood
expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)
expect_equal(logLik(e.lmm2), as.double(logLik(e.lmm)), tol = 1e-6)
## coef(e.lmm, transform.sigma = "log", transform.k = "log")
## coef(e.lmm, transform.sigma = "square", transform.k = "square")
## coef(e.lmm, transform.sigma = "logsquare", transform.k p= "logsquare")

## ** score
test <- score(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none")
GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"))
expect_true(all(abs(test) < 1e-6))
expect_equal(as.double(GS),as.double(test))

## no transformation
newp <- coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none")+1
GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## log transformation
newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"]); newp.log["k.F"] <- log(newp["k.F"])
GS <- jacobian(func = function(p){p[c("sigma","k.F")] <- exp(p[c("sigma","k.F")]); logLik(e.lmm, p = p)}, x = newp.log)
test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "log", transform.k = "log")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## lava transformation
newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2; newp.2["k.F"] <- newp["k.F"]^2*newp["sigma"]^2
GS <- jacobian(func = function(p){p[c("sigma","k.F")] <- c(sqrt(p["sigma"]),sqrt(p["k.F"]/p["sigma"])); logLik(e.lmm, p = p, transform.sigma = "none")}, x = newp.2)
test <- score(e.lmm, effects = "all", p = newp, transform.k = "var")
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** information
test0 <- information(e.lmm, effects = "all",
                     p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"),
                     transform.sigma = "none", transform.k = "none", type.information = "expected") ## using sigma
                     test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"),
                                         transform.sigma = "none", transform.k = "none", type.information = "observed") ## using sigma
GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test0 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.sigma = "log", transform.k = "log", type.information = "expected") ## using log(sigma) and log(k)
test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.sigma = "log", transform.k = "log", type.information = "observed") ## using log(sigma) and log(k)
GS <- -hessian(func = function(p){p[c("sigma","k.F")] <- exp(p[c("sigma","k.F")]);logLik(e.lmm, p = p, transform.sigma = "none", transform.k = "none")}, x = coef(e.lmm, effects = "all", transform.sigma = "log", transform.k = "log", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

test0 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.k = "var", type.information = "expected") ## using sigma^2 and sigma^2 k^2
test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.k = "var", type.information = "observed") ## using sigma^2 and sigma^2 k^2
GS <- -hessian(func = function(p){p[c("sigma","k.F")] <- c(sqrt(p["sigma"]),sqrt(p["k.F"]/p["sigma"]));logLik(e.lmm, p = p, transform.sigma = "none", transform.k = "none")}, x = coef(e.lmm, effects = "all", transform.k = "var", transform.names = FALSE))
expect_equal(as.double(test), as.double(GS), tol = 1e-6)

## ** degree of freedom
test <- confint(e.lmm, effects = "all", transform.k = "logsd", type.information = "observed", backtransform = TRUE)[,"df",drop=FALSE]
name.coefM <- grep(rownames(test),pattern="M",value=TRUE)
name.coefF <- grep(rownames(test),pattern="F",value=TRUE)
expect_equal(test[name.coefM,], rep(sum(d$Gender=="M")-length(name.coefM)+1,length(name.coefM)), tol = 1e-6)
expect_equal(test[name.coefF,], rep(sum(d$Gender=="F")-length(name.coefM)+1,length(name.coefM)), tol = 1e-6)

## ** confidence interval
capture.output(summary(e.lmm))
capture.output(anova(e.lmm))
## anova(e.gls,type="marginal")

## ** predictions
if(require(AICcmodavg)){
    GS <- AICcmodavg::predictSE(e.gls2, newdata = d)
    pp1 <- predict(e.lmm, newdata = d)
    pp2 <- predict(e.lmm2, newdata = d)
    set.seed(10)
    index <- sample.int(NROW(d))
    pp3 <- predict(e.lmm, newdata = d[index,,drop=FALSE])
    expect_equivalent(pp1$estimate,GS$fit, tol = 1e-5)
    expect_equivalent(pp1$se,GS$se.fit, tol = 1e-5)
    expect_equal(pp1,pp2, tol = 1e-5)
    expect_equivalent(pp1[index,,drop=FALSE],pp3, tol = 1e-5) ## different rownames
    ## .getUVarCov(e.lmm)
}

})

## * missing values
test_that("missing values",{
    set.seed(11)
    d$Ymiss <- d$Y
    d$Ymiss[which(rbinom(NROW(d), size = 1, prob = 0.1)==1)] <- NA
    
    ## ** fit
    e.lmm <- lmm(Ymiss ~ X1 + X2 + Gene, repetition = ~time|id, structure = "CS", data = d, trace = 0,
                 method = "ML", df = TRUE)
    e.gls <- gls(Ymiss ~ X1 + X2 + Gene, data = d, method = "ML", na.action = na.omit)
    expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)
})


##----------------------------------------------------------------------
### test-auto-linear-regression.R ends here
