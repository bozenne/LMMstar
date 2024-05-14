### test-auto-linear-regression.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (10:13) 
## Version: 
## Last-Updated: May 12 2024 (20:37) 
##           By: Brice Ozenne
##     Update #: 222
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
    library(nlme)

    library(LMMstar)
}

context("Check lmm on linear regressions")
LMMstar.options(optimizer = "FS", method.numDeriv = "Richardson", precompute.moments = TRUE,
                columns.confint = c("estimate","se","df","lower","upper","p.value"), df = TRUE)

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

    ## ** iteration
    expect_equal(e.lmm$opt$n.iter,0)

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
    ## GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
    GS <- rbind(c(-25.87627956, -15.79878208, -11.78725113, -12.99280729, -7.23996332, 18.51219198))
    test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "none")
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    ## log transformation
    newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"])
    ## GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p)}, x = newp.log)
    GS <- rbind(c(-25.87627956, -15.79878208, -11.78725113, -12.99280729, -7.23996332, 34.25729653))
    test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "log")
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    ## lava transformation
    newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2
    ## GS <- jacobian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(e.lmm, p = p)}, x = newp.2)
    GS <- rbind(c(-25.87627956, -15.79878208, -11.78725113, -12.99280729, -7.23996332, 5.0018724))
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

    ## GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none"))
    GS <- cbind(c(69.11854669, 4.7137157, 1.66171656, 30.41216054, 16.58845121, 0), 
                c(4.7137157, 64.92245702, -3.26776856, 2.49477216, 5.92598878, 0), 
                c(1.66171656, -3.26776856, 64.04816027, -1.81300977, -4.82994473, 0), 
                c(30.41216054, 2.49477216, -1.81300977, 30.41216054, 0, 0), 
                c(16.58845121, 5.92598878, -4.82994473, 0, 16.58845121, 0), 
                c(0, 0, 0, 0, 0, 138.23709338)
                )
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "log")
    ## GS <- -hessian(func = function(p){p["sigma"]<-exp(p["sigma"]);logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "log", transform.names = FALSE))
    GS <- cbind(c(69.11854669, 4.7137157, 1.66171656, 30.41216054, 16.58845121, 0), 
                c(4.7137157, 64.92245702, -3.26776856, 2.49477216, 5.92598878, 0), 
                c(1.66171656, -3.26776856, 64.04816027, -1.81300977, -4.82994473, 0), 
                c(30.41216054, 2.49477216, -1.81300977, 30.41216054, 0, 0), 
                c(16.58845121, 5.92598878, -4.82994473, 0, 16.58845121, 0), 
                c(0, 0, 0, 0, 0, 100)
                )
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "square") 
    ## GS <- -hessian(func = function(p){p["sigma"]<-sqrt(p["sigma"]);logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "square", transform.names = FALSE))
    GS <- cbind(c(69.11854669, 4.7137157, 1.66171656, 30.41216054, 16.58845121, 0), 
                c(4.7137157, 64.92245702, -3.26776856, 2.49477216, 5.92598878, 0), 
                c(1.66171656, -3.26776856, 64.04816027, -1.81300977, -4.82994473, 0), 
                c(30.41216054, 2.49477216, -1.81300977, 30.41216054, 0, 0), 
                c(16.58845121, 5.92598878, -4.82994473, 0, 16.58845121, 0), 
                c(0, 0, 0, 0, 0, 47.77373496)
                )
    ## GS0 <- -hessian(func = function(p){logLik(e.lava, p = p)}, x = coef(e.lava))
    GS0 <- cbind(c(69.11854669, 4.7137157, 1.66171656, 30.41216054, 16.58845121, 0), 
                 c(4.7137157, 64.92245702, -3.26776856, 2.49477216, 5.92598878, 0), 
                 c(1.66171656, -3.26776856, 64.04816027, -1.81300977, -4.82994473, 0), 
                 c(30.41216054, 2.49477216, -1.81300977, 30.41216054, 0, 0), 
                 c(16.58845121, 5.92598878, -4.82994473, 0, 16.58845121, 0), 
                 c(0, 0, 0, 0, 0, 47.77373496)
                 )
    expect_equal(as.double(test["sigma^2","sigma^2"]), as.double(nobs(e.lmm)[1]/(2*coef(e.lmm, effect = "variance", transform.sigma = "square")^2)), tol = 1e-6)
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)
    expect_equal(as.double(GS0), as.double(GS), tol = 1e-6)

    ## no transformation 
    newp <- coef(e.lmm, effects = "all", transform.sigma = "none")+1
    ## GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
    GS <- cbind(c(14.60089647, 0.9957454, 0.35102809, 6.42439445, 3.50421515, -27.96640152), 
                c(0.9957454, 13.7144966, -0.69029736, 0.52700631, 1.25183113, -17.07490763), 
                c(0.35102809, -0.69029736, 13.52980643, -0.38298791, -1.02029812, -12.73935061), 
                c(6.42439445, 0.52700631, -0.38298791, 6.42439445, 0, -14.04228396), 
                c(3.50421515, 1.25183113, -1.02029812, 0, 3.50421515, -7.8247617), 
                c(-27.96640152, -17.07490763, -12.73935061, -14.04228396, -7.8247617, 59.21302733)
                )
    test <- information(e.lmm, effects = "all", p = newp, transform.sigma = "none")
    expect_equal(as.double(test[1:4,1:4]), as.double(GS[1:4,1:4]), tol = 1e-6) ## does not match as some terms do not cancel

    ## ** variance-covariance
    test <- vcov(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "none")
    ## GS <- solve(-hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none")))
    GS <- cbind(c(0.04803812, 0.002704, -0.00631934, -0.04863666, -0.05084404, 0), 
                c(0.002704, 0.01613436, 4.7e-07, -0.00402751, -0.00846764, 0), 
                c(-0.00631934, 4.7e-07, 0.0168306, 0.00732265, 0.01121962, 0), 
                c(-0.04863666, -0.00402751, 0.00732265, 0.08228516, 0.05220751, 0), 
                c(-0.05084404, -0.00846764, 0.01121962, 0.05220751, 0.11741863, 0), 
                c(0, 0, 0, 0, 0, 0.00723395)
                )
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    test <- vcov(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "square", df = TRUE)
    ## GS <- solve(-hessian(func = function(p){p["sigma"] <- sqrt(p["sigma"]);logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "square", transform.names = FALSE)))
    GS <- cbind(c(0.04803812, 0.002704, -0.00631934, -0.04863666, -0.05084404, 0), 
                c(0.002704, 0.01613436, 4.7e-07, -0.00402751, -0.00846764, 0), 
                c(-0.00631934, 4.7e-07, 0.0168306, 0.00732265, 0.01121962, 0), 
                c(-0.04863666, -0.00402751, 0.00732265, 0.08228516, 0.05220751, 0), 
                c(-0.05084404, -0.00846764, 0.01121962, 0.05220751, 0.11741863, 0), 
                c(0, 0, 0, 0, 0, 0.020932)
                )
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    expect_equal(unname(sqrt(diag(vcov(e.lmm, effects = "all", transform.sigma="square",robust=TRUE)))), unname(estimate(e.lava, robust = TRUE)$coefmat[,"Std.Err"]),
                 tol = 1e-4)

    ## ** degree of freedom
    test <- model.tables(e.lmm, effects = "all", transform.sigma = "log", type.information = "observed")$df
    expect_equal(test, rep(n.obs,n.param), tol = 1e-6)

    test <- model.tables(e.lmm, effects = "all", transform.sigma = "square", type.information = "expected")$df
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
    ## summary(test)
    ## summary(anova(e.lmm, effects = "all"))
   
    expect_equal(prod(test$multivariate[test$multivariate$test=="Gene",c("statistic","df.num")]),
                 unname(lava::compare(e.lava, par = c("Y~GeneLA","Y~GeneAA"))$statistic), tol = 1e-6)
    expect_equal(test$multivariate[test$multivariate$test=="Gene","df.denom"], NROW(d), tol = 1e-6)

    test2 <- anova(e.lmm, effect = c("GeneLA=0","GeneAA=0"))
    ## summary(test2)
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

    ## ** iteration
    expect_equal(e.lmm$opt$n.iter,0)

    ## ** coef
    expect_equal(coef(e.lmm, effects = "mean"), coef(e.gls), tol = 1e-6)

    ## ** logLikelihood
    expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)

    ## ** score
    expect_true(all(abs(score(e.lmm, effects = "all")) < 1e-6))

    ## no transformation
    newp <- coef(e.lmm, effects = "all", transform.sigma = "none")+1
    ## GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
    GS <- rbind(c(-24.6360655, -15.04156844, -11.22230459, -12.37008012, -6.89296196, 18.57017787))
    test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "none")
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    ## transformation
    newp.log <- newp; newp.log["sigma"] <- log(newp.log["sigma"])
    ## GS <- jacobian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p)}, x = newp.log)
    GS <- rbind(c(-24.6360655, -15.04156844, -11.22230459, -12.37008012, -6.89296196, 35.2189616))
    test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "log")
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    ## ** information
    test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "none", type.information = "observed")
    ## GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none"))
    GS <- cbind(c(62.20669202, 4.24234413, 1.4955449, 27.37094449, 14.92960608, 0), 
                c(4.24234413, 58.43021132, -2.94099171, 2.24529494, 5.33338991, 0), 
                c(1.4955449, -2.94099171, 57.64334425, -1.6317088, -4.34695026, 0), 
                c(27.37094449, 2.24529494, -1.6317088, 27.37094449, 0, 0), 
                c(14.92960608, 5.33338991, -4.34695026, 0, 14.92960608, 0), 
                c(0, 0, 0, 0, 0, 111.97204564)
                )

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
    ## GS <- -hessian(func = function(p){p["sigma"] <- exp(p["sigma"]); logLik(e.lmm, p = p)}, x = coef(e.lmm, transform.sigma = "log", transform.names = FALSE))
    GS <- cbind(c(62.20669202, 4.24234413, 1.4955449, 27.37094449, 14.92960608, 0), 
                c(4.24234413, 58.43021132, -2.94099171, 2.24529494, 5.33338991, 0), 
                c(1.4955449, -2.94099171, 57.64334425, -1.6317088, -4.34695026, 0), 
                c(27.37094449, 2.24529494, -1.6317088, 27.37094449, 0, 0), 
                c(14.92960608, 5.33338991, -4.34695026, 0, 14.92960608, 0), 
                c(0, 0, 0, 0, 0, 89.99999999)
                )
    expect_equal(as.double(test1), as.double(GS), tol = 1e-6)

    test1 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "square", type = "observed")
    test2 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none"), transform.sigma = "square", type = "expected")
    ## GS <- -hessian(func = function(p){p["sigma"] <- sqrt(p["sigma"]); logLik(e.lmm, p = p)}, x = coef(e.lmm, transform.sigma = "square", transform.names = FALSE))
    GS <- cbind(c(62.20669202, 4.24234413, 1.4955449, 27.37094449, 14.92960608, 0), 
                c(4.24234413, 58.43021132, -2.94099171, 2.24529494, 5.33338991, 0), 
                c(1.4955449, -2.94099171, 57.64334425, -1.6317088, -4.34695026, 0), 
                c(27.37094449, 2.24529494, -1.6317088, 27.37094449, 0, 0), 
                c(14.92960608, 5.33338991, -4.34695026, 0, 14.92960608, 0), 
                c(0, 0, 0, 0, 0, 34.82705279)
                )
    expect_equal(as.double(test1), as.double(GS), tol = 1e-6)
    expect_equal(as.double(test1["sigma^2","sigma^2"]), as.double((nobs(e.lmm)[1]-length(coef(e.lmm, effect = "mean")))/(2*coef(e.lmm, effect = "variance", transform.sigma = "square")^2)), tol = 1e-6)

    ## no transformation 
    newp <- coef(e.lmm, effects = "all", transform.sigma = "none")+1
    ## GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
    GS <- cbind(c(13.90109583, 0.94802071, 0.3342038, 6.11648216, 3.336263, -25.98010262), 
                c(0.94802071, 13.05717987, -0.65721237, 0.50174763, 1.19183261, -15.86217132), 
                c(0.3342038, -0.65721237, 12.88134164, -0.36463184, -0.97139665, -11.83454496), 
                c(6.11648216, 0.50174763, -0.36463184, 6.11648216, 0, -13.04493816), 
                c(3.336263, 1.19183261, -0.97139665, 0, 3.336263, -7.26901213), 
                c(-25.98010262, -15.86217132, -11.83454496, -13.04493816, -7.26901213, 54.3969021)
                )
    test <- information(e.lmm, effects = "all", p = newp, transform.sigma = "none")
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    ## ** degree of freedom
    test <- model.tables(e.lmm, effects = "all", transform.sigma = "log", type.information = "observed")$df
    expect_equal(test, rep(n.obs-n.mu,n.param), tol = 1e-6)

    ## numerical derivative
    test <- model.tables(e.lmm, effects = "all", transform.sigma = "log", type.information = "observed")$df
    ## FF.bis <- function(p){p["sigma"] <- exp(p["sigma"])   ; diag(vcov(e.lmm, effects = "all", p = p, transform.sigma = "log", type.information = "observed"))}
    ## GG.bis <- jacobian(func = FF.bis, x = coef(e.lmm, effects = "all", transform.sigma = "log", transform.names = FALSE), method = LMMstar.options()$method.numDeriv)
    GG.bis <- cbind(c(0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0), 
                    c(0.10675137, 0.03585414, 0.03740132, 0.18285592, 0.26093029, 0.02222222)
                    )
    VV.bis <- vcov(e.lmm, effects = "all", transform.sigma = "log", type.information = "observed")
    GS <- sapply(1:NROW(GG.bis),function(gg){2*VV.bis[gg,gg]^2/ (GG.bis[gg,,drop=FALSE] %*% VV.bis %*% t(GG.bis[gg,,drop=FALSE]))})
    expect_equal(unname(test), unname(GS), tol = 1e-6)

    ## test <- vcov(e.lmm, df = 2, transform.sigma = "log", type.information = "observed")
    ## all(abs(test-VV.bis)<1e-10)
    ## apply(attr(test,"dVcov"),3,diag)
    ## GG.bis

    ## numerical derivative on another scale
    test <- model.tables(e.lmm, effects = "all", transform.sigma = "square", type.information = "observed")$df
    ## FF.bis <- function(p){p["sigma"] <- sqrt(p["sigma"]); diag(vcov(e.lmm, effects = "all", p = p, transform.sigma = "square", type.information = "observed"))}
    ## GG.bis <- jacobian(func = FF.bis, x = coef(e.lmm, effects = "all", transform.sigma = "square", transform.names = FALSE), method = LMMstar.options()$method.numDeriv)
    GG.bis <- cbind(c(0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0), 
                    c(0, 0, 0, 0, 0, 0), 
                    c(0.0664065, 0.02230367, 0.02326613, 0.11374862, 0.1623161, 0.14289281)
                    )
    VV.bis <- vcov(e.lmm, effects = "all", transform.sigma = "square", type.information = "observed")
    GS <- sapply(1:NROW(GG.bis),function(gg){2*VV.bis[gg,gg]^2/ (GG.bis[gg,,drop=FALSE] %*% VV.bis %*% t(GG.bis[gg,,drop=FALSE]))})
    expect_equal(unname(test), unname(GS), tol = 1e-6)

    ## ** anova
    test <- anova(e.lmm)
    GS <- anova(e.gls, type = "marginal")
    expect_equal(test$multivariate$statistic,GS[["F-value"]][-1], tol = 1e-6)
    expect_equal(test$multivariate$p.value,GS[["p-value"]][-1], tol = 1e-6)

    ## ** predictions
    test <- predict(e.lmm, newdata = d, se = TRUE)
    GS <- predict(e.lm, newdata = d, se = TRUE)
    expect_equal(test$estimate,as.double(GS$fit), tol = 1e-7)
    expect_equal(test$se,as.double(GS$se.fit), tol = 1e-7)
    expect_true(all(abs(test$df-GS$df) < 1e-7))

})


## * multiple variance parameters (ML)
test_that("multiple variance parameter (ML)",{

    ## ** fit
    e.lmm <- lmm(Y ~ -1 + Gender + (X1 + X2 + Gene):Gender, repetition = ~Gender|id, structure = "IND", data = d, trace = 0,
                 method = "ML")
    e.gls <- gls(Y ~ -1 + Gender + (X1 + X2 + Gene):Gender, data = d, weights = varIdent(form=~1|Gender), method = "ML")
    e.lmm2 <- lmm(Y ~(X1 + X2 + Gene)*Gender, repetition = Gender~time|id, structure = "ID", data = d, trace = 0,
                  method = "ML")

    n.obs <- unname(nobs(e.lmm)[1])
    n.mu <- length(coef(e.lmm, effects = "mean"))
    n.sigma <- length(coef(e.lmm, effects = "variance"))
    n.param <- length(coef(e.lmm, effects = "all"))

    ## ** iteration
    expect_equal(e.lmm$opt$n.iter,0)
    expect_equal(e.lmm2$opt$n.iter,0)

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
    ## GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"))
    GS <- rbind(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)
    expect_true(all(abs(test) < 1e-4))

    ## no transformation
    newp <- coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none")+1
    ## GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
    GS <- rbind(c(-12.59534519, -3.67518576, -4.45620371, -3.12477721, -5.02017891, -1.86951941, -7.25091419, -1.59294925, -4.04666275, -0.88593539, 0.42454894, -7.10232389))
    test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none")
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    ## log transformation
    newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"]); newp.log["k.F"] <- log(newp["k.F"])
    ## GS <- jacobian(func = function(p){p[c("sigma","k.F")] <- exp(p[c("sigma","k.F")]); logLik(e.lmm, p = p)}, x = newp.log)
    GS <- rbind(c(-12.59534519, -3.67518576, -4.45620371, -3.12477721, -5.02017891, -1.86951941, -7.25091419, -1.59294925, -4.04666275, -0.88593539, 0.78245096, -13.60825916))
    test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "log", transform.k = "log")
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    ## lava transformation
    newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2; newp.2["k.F"] <- newp["k.F"]^2*newp["sigma"]^2
    ## GS <- jacobian(func = function(p){p[c("sigma","k.F")] <- c(sqrt(p["sigma"]),sqrt(p["k.F"]/p["sigma"])); logLik(e.lmm, p = p, transform.sigma = "none")}, x = newp.2)
    GS <- rbind(c(-12.59534519, -3.67518576, -4.45620371, -3.12477721, -5.02017891, -1.86951941, -7.25091419, -1.59294925, -4.04666275, -0.88593539, 2.11832942, -0.54564433))
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
                         transform.sigma = "none", transform.k = "none", type.information = "expected")
                    
    test <- information(e.lmm, effects = "all",
                        p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"),
                        transform.sigma = "none", transform.k = "none", type.information = "observed")
    ## GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"))
    GS <- cbind(c(33.77057622, 0, 0.88985514, 0, -1.19553421, 0, 16.88528811, 0, 9.8497514, 0, 0, 0), 
                c(0, 43.59955426, 0, 4.65756814, 0, 3.44053616, 0, 16.76905933, 0, 8.38452966, 0, 0), 
                c(0.88985514, 0, 21.70660151, 0, -3.46195721, 0, 0.25820444, 0, 1.90589312, 0, 0, 0), 
                c(0, 4.65756814, 0, 52.88654407, 0, 0.16173647, 0, 2.7186085, 0, 4.91728133, 0, 0), 
                c(-1.19553421, 0, -3.46195721, 0, 30.28860275, 0, 0.62724252, 0, -2.2642157, 0, 0, 0), 
                c(0, 3.44053616, 0, 0.16173647, 0, 41.59845209, 0, -2.94680881, 0, -3.16068124, 0, 0), 
                c(16.88528811, 0, 0.25820444, 0, 0.62724252, 0, 16.88528811, 0, 0, 0, 0, 0), 
                c(0, 16.76905933, 0, 2.7186085, 0, -2.94680881, 0, 16.76905933, 0, 0, 0, 0), 
                c(9.8497514, 0, 1.90589312, 0, -2.2642157, 0, 0, 0, 9.8497514, 0, 0, 0), 
                c(0, 8.38452966, 0, 4.91728133, 0, -3.16068124, 0, 0, 0, 8.38452966, 0, 0), 
                c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 140.71073425, 67.33760941), 
                c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 67.33760941, 61.97047366)
                )
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    test0 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.sigma = "log", transform.k = "log", type.information = "expected") ## using log(sigma) and log(k)
    test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.sigma = "log", transform.k = "log", type.information = "observed") ## using log(sigma) and log(k)
    ## GS <- -hessian(func = function(p){p[c("sigma","k.F")] <- exp(p[c("sigma","k.F")]);logLik(e.lmm, p = p, transform.sigma = "none", transform.k = "none")}, x = coef(e.lmm, transform.sigma = "log", transform.k = "log", transform.names = FALSE))
    GS <- cbind(c(33.77057622, 0, 0.88985514, 0, -1.19553421, 0, 16.88528811, 0, 9.8497514, 0, 0, 0), 
                c(0, 43.59955426, 0, 4.65756814, 0, 3.44053616, 0, 16.76905933, 0, 8.38452966, 0, -1e-08), 
                c(0.88985514, 0, 21.70660151, 0, -3.46195721, 0, 0.25820444, 0, 1.90589312, 0, 0, 0), 
                c(0, 4.65756814, 0, 52.88654407, 0, 0.16173647, 0, 2.7186085, 0, 4.91728133, 0, 0), 
                c(-1.19553421, 0, -3.46195721, 0, 30.28860275, 0, 0.62724252, 0, -2.2642157, 0, 0, 0), 
                c(0, 3.44053616, 0, 0.16173647, 0, 41.59845209, 0, -2.94680881, 0, -3.16068124, 0, 0), 
                c(16.88528811, 0, 0.25820444, 0, 0.62724252, 0, 16.88528811, 0, 0, 0, 0, 0), 
                c(0, 16.76905933, 0, 2.7186085, 0, -2.94680881, 0, 16.76905933, 0, 0, 0, 0), 
                c(9.8497514, 0, 1.90589312, 0, -2.2642157, 0, 0, 0, 9.8497514, 0, 0, 0), 
                c(0, 8.38452966, 0, 4.91728133, 0, -3.16068124, 0, 0, 0, 8.38452966, 0, 0), 
                c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 99.99999999, 52), 
                c(0, -1e-08, 0, 0, 0, 0, 0, 0, 0, 0, 52, 52)
                )
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    test0 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.k = "var", type.information = "expected") ## using sigma^2 and sigma^2 k^2
    test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.k = "var", type.information = "observed") ## using sigma^2 and sigma^2 k^2
    ## GS <- -hessian(func = function(p){p[c("sigma","k.F")] <- c(sqrt(p["sigma"]),sqrt(p["k.F"]/p["sigma"]));logLik(e.lmm, p = p, transform.sigma = "none", transform.k = "none")}, x = coef(e.lmm, transform.k = "var", transform.names = FALSE))
    GS <- cbind(c(33.77057622, 0, 0.88985514, 0, -1.19553421, 0, 16.88528811, 0, 9.8497514, 0, 0, 0), 
                c(0, 43.59955426, 0, 4.65756814, 0, 3.44053616, 0, 16.76905933, 0, 8.38452966, 0, 0), 
                c(0.88985514, 0, 21.70660151, 0, -3.46195721, 0, 0.25820444, 0, 1.90589312, 0, 0, 0), 
                c(0, 4.65756814, 0, 52.88654407, 0, 0.16173647, 0, 2.7186085, 0, 4.91728133, 0, 0), 
                c(-1.19553421, 0, -3.46195721, 0, 30.28860275, 0, 0.62724252, 0, -2.2642157, 0, 0, 0), 
                c(0, 3.44053616, 0, 0.16173647, 0, 41.59845209, 0, -2.94680881, 0, -3.16068124, 0, 0), 
                c(16.88528811, 0, 0.25820444, 0, 0.62724252, 0, 16.88528811, 0, 0, 0, 0, 0), 
                c(0, 16.76905933, 0, 2.7186085, 0, -2.94680881, 0, 16.76905933, 0, 0, 0, 0), 
                c(9.8497514, 0, 1.90589312, 0, -2.2642157, 0, 0, 0, 9.8497514, 0, 0, 0), 
                c(0, 8.38452966, 0, 4.91728133, 0, -3.16068124, 0, 0, 0, 8.38452966, 0, 0), 
                c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 23.75941288, 0), 
                c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 36.5561756)
                )
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    ## ** degree of freedom
    test <- model.tables(e.lmm, effects = "all", transform.k = "logsd", type.information = "observed", backtransform = TRUE)
    name.coefM <- grep(rownames(test),pattern="M",value=TRUE)
    name.coefF <- grep(rownames(test),pattern="F",value=TRUE)
    expect_equal(test[name.coefM,"df"], rep(sum(d$Gender=="M"),length(name.coefM)), tol = 1e-3)
    expect_equal(test[name.coefF,"df"], rep(sum(d$Gender=="F"),length(name.coefM)), tol = 1e-3)

    test <- model.tables(e.lmm, effects = "all", transform.k = "var", type.information = "expected")
    name.coefM <- grep(rownames(test),pattern="M",value=TRUE)
    name.coefF <- grep(rownames(test),pattern="F",value=TRUE)
    expect_equal(test[name.coefM,"df"], c(rep(sum(d$Gender=="M"),length(name.coefM)-1),sum(d$Gender=="M")/4), tol = 1e-6)
    expect_equal(test[name.coefF,"df"], c(rep(sum(d$Gender=="F"),length(name.coefM)-1),sum(d$Gender=="F")/4), tol = 1e-6)

    ## ** variance-covariance
    sigma(e.lmm)
    sigma(e.lmm2)

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
    capture.output(anova(e.lmm))
})

## * multiple variance parameters (REML)

test_that("multiple variance parameters (REML)",{

    ## ** fit
    e.lmm <- lmm(Y ~ -1 + Gender + (X1 + X2 + Gene):Gender, repetition = ~Gender|id, structure = "IND", data = d, trace = 0,
                 method = "REML")
    e.lmm2 <- lmm(Y ~(X1 + X2 + Gene)*Gender, repetition = Gender~time|id, structure = "ID", data = d, trace = 0,
                  method = "REML")
    e.gls <- gls(Y ~ -1 + Gender + (X1 + X2 + Gene):Gender, data = d, weights = varIdent(form=~1|Gender), method = "REML")
    e.gls2 <- gls(Y ~ (X1 + X2 + Gene)*Gender, data = d, weights = varIdent(form=~1|Gender), method = "REML")

    n.obs <- unname(nobs(e.lmm)[1])
    n.mu <- length(coef(e.lmm, effects = "mean"))
    n.sigma <- length(coef(e.lmm, effects = "variance"))
    n.param <- length(coef(e.lmm, effects = "all"))

    ## ** iteration
    expect_equal(e.lmm$opt$n.iter,0)
    expect_equal(e.lmm2$opt$n.iter,0)

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
    ## GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"))
    GS <- rbind(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    expect_true(all(abs(test) < 1e-6))
    expect_equal(as.double(GS),as.double(test))

    ## no transformation
    newp <- coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none")+1
    ## GS <- jacobian(func = function(p){logLik(e.lmm, p = p)}, x = newp)
    GS <- rbind(c(-11.28047812, -3.32312688, -3.99100681, -2.82544388, -4.4961069, -1.69043162, -6.49396881, -1.44035508, -3.6242191, -0.80106854, 2.86919921, -5.13678117))
    test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "none", transform.k = "none")
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    ## log transformation
    newp.log <- newp; newp.log["sigma"] <- log(newp["sigma"]); newp.log["k.F"] <- log(newp["k.F"])
    ## GS <- jacobian(func = function(p){p[c("sigma","k.F")] <- exp(p[c("sigma","k.F")]); logLik(e.lmm, p = p)}, x = newp.log)
    GS <- rbind(c(-11.28047812, -3.32312688, -3.99100681, -2.82544388, -4.4961069, -1.69043162, -6.49396881, -1.44035508, -3.6242191, -0.80106854, 5.58767775, -9.79530729))
    test <- score(e.lmm, effects = "all", p = newp, transform.sigma = "log", transform.k = "log")
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    ## lava transformation
    newp.2 <- newp; newp.2["sigma"] <- newp["sigma"]^2; newp.2["k.F"] <- newp["k.F"]^2*newp["sigma"]^2
    ## GS <- jacobian(func = function(p){p[c("sigma","k.F")] <- c(sqrt(p["sigma"]),sqrt(p["k.F"]/p["sigma"])); logLik(e.lmm, p = p, transform.sigma = "none")}, x = newp.2)
    GS <- rbind(c(-11.28047812, -3.32312688, -3.99100681, -2.82544388, -4.4961069, -1.69043162, -6.49396881, -1.44035508, -3.6242191, -0.80106854, 2.02800642, -0.35513442))
    test <- score(e.lmm, effects = "all", p = newp, transform.k = "var")
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    ## ** information
    test0 <- information(e.lmm, effects = "all",
                         p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"),
                         transform.sigma = "none", transform.k = "none", type.information = "expected")
                         
    test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"),
                        transform.sigma = "none", transform.k = "none", type.information = "observed")
    ## GS <- -hessian(func = function(p){logLik(e.lmm, p = p)}, x = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"))
    GS <- cbind(c(26.73503951, 0, 0.70446865, 0, -0.94646459, 0, 13.36751975, 0, 7.79771986, 0, 0, 0), 
                c(0, 35.21502459, 0, 3.76188196, 0, 2.77889459, 0, 13.54424023, 0, 6.77212011, 0, 0), 
                c(0.70446865, 0, 17.18439286, 0, -2.74071613, 0, 0.20441185, 0, 1.50883205, 0, 0, 0), 
                c(0, 3.76188196, 0, 42.71605483, 0, 0.1306333, 0, 2.19579918, 0, 3.97165031, 0, 0), 
                c(-0.94646459, 0, -2.74071613, 0, 23.97847718, 0, 0.496567, 0, -1.79250409, 0, 0, 0), 
                c(0, 2.77889459, 0, 0.1306333, 0, 33.59874977, 0, -2.38011481, 0, -2.55285792, 0, 0), 
                c(13.36751975, 0, 0.20441185, 0, 0.496567, 0, 13.36751975, 0, 0, 0, 0, 0), 
                c(0, 13.54424023, 0, 2.19579918, 0, -2.38011481, 0, 13.54424023, 0, 0, 0, 0), 
                c(7.79771986, 0, 1.50883205, 0, -1.79250409, 0, 0, 0, 7.79771986, 0, 0, 0), 
                c(0, 6.77212011, 0, 3.97165031, 0, -2.55285792, 0, 0, 0, 6.77212011, 0, 0), 
                c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 89.11679836, 48.8794842), 
                c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 48.8794842, 51.06629502)
                )
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    test0 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.sigma = "log", transform.k = "log", type.information = "expected") ## using log(sigma) and log(k)
    test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.sigma = "log", transform.k = "log", type.information = "observed") ## using log(sigma) and log(k)
    ## GS <- -hessian(func = function(p){p[c("sigma","k.F")] <- exp(p[c("sigma","k.F")]);logLik(e.lmm, p = p, transform.sigma = "none", transform.k = "none")}, x = coef(e.lmm, effects = "all", transform.sigma = "log", transform.k = "log", transform.names = FALSE))
    GS <- cbind(c(26.73503951, 0, 0.70446865, 0, -0.94646459, 0, 13.36751975, 0, 7.79771986, 0, 0, 0), 
                c(0, 35.21502459, 0, 3.76188196, 0, 2.77889459, 0, 13.54424023, 0, 6.77212011, 1e-08, 0), 
                c(0.70446865, 0, 17.18439286, 0, -2.74071613, 0, 0.20441185, 0, 1.50883205, 0, 0, 0), 
                c(0, 3.76188196, 0, 42.71605483, 0, 0.1306333, 0, 2.19579918, 0, 3.97165031, 0, 0), 
                c(-0.94646459, 0, -2.74071613, 0, 23.97847718, 0, 0.496567, 0, -1.79250409, 0, 0, 0), 
                c(0, 2.77889459, 0, 0.1306333, 0, 33.59874977, 0, -2.38011481, 0, -2.55285792, 0, 0), 
                c(13.36751975, 0, 0.20441185, 0, 0.496567, 0, 13.36751975, 0, 0, 0, 0, 0), 
                c(0, 13.54424023, 0, 2.19579918, 0, -2.38011481, 0, 13.54424023, 0, 0, 0, 0), 
                c(7.79771986, 0, 1.50883205, 0, -1.79250409, 0, 0, 0, 7.79771986, 0, 0, 0), 
                c(0, 6.77212011, 0, 3.97165031, 0, -2.55285792, 0, 0, 0, 6.77212011, 0, 0), 
                c(0, 1e-08, 0, 0, 0, 0, 0, 0, 0, 0, 79.99999991, 42.00000002), 
                c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 42.00000002, 41.99999999)
                )
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    test0 <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.k = "var", type.information = "expected") ## using sigma^2 and sigma^2 k^2
    test <- information(e.lmm, effects = "all", p = coef(e.lmm, effects = "all", transform.sigma = "none", transform.k = "none"), transform.k = "var", type.information = "observed") ## using sigma^2 and sigma^2 k^2
    ## GS <- -hessian(func = function(p){p[c("sigma","k.F")] <- c(sqrt(p["sigma"]),sqrt(p["k.F"]/p["sigma"]));logLik(e.lmm, p = p, transform.sigma = "none", transform.k = "none")}, x = coef(e.lmm, effects = "all", transform.k = "var", transform.names = FALSE))
    GS <- cbind(c(26.73503951, 0, 0.70446865, 0, -0.94646459, 0, 13.36751975, 0, 7.79771986, 0, 0, 0), 
                c(0, 35.21502459, 0, 3.76188196, 0, 2.77889459, 0, 13.54424023, 0, 6.77212011, 0, 0), 
                c(0.70446865, 0, 17.18439286, 0, -2.74071613, 0, 0.20441185, 0, 1.50883205, 0, 0, 0), 
                c(0, 3.76188196, 0, 42.71605483, 0, 0.1306333, 0, 2.19579918, 0, 3.97165031, 0, 0), 
                c(-0.94646459, 0, -2.74071613, 0, 23.97847718, 0, 0.496567, 0, -1.79250409, 0, 0, 0), 
                c(0, 2.77889459, 0, 0.1306333, 0, 33.59874977, 0, -2.38011481, 0, -2.55285792, 0, 0), 
                c(13.36751975, 0, 0.20441185, 0, 0.496567, 0, 13.36751975, 0, 0, 0, 0, 0), 
                c(0, 13.54424023, 0, 2.19579918, 0, -2.38011481, 0, 13.54424023, 0, 0, 0, 0), 
                c(7.79771986, 0, 1.50883205, 0, -1.79250409, 0, 0, 0, 7.79771986, 0, 0, 0), 
                c(0, 6.77212011, 0, 3.97165031, 0, -2.55285792, 0, 0, 0, 6.77212011, 0, 0), 
                c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11.78861494, 0), 
                c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 19.26187655)
                )
    expect_equal(as.double(test), as.double(GS), tol = 1e-6)

    ## ** degree of freedom
    test <- model.tables(e.lmm, effects = "all", transform.k = "logsd", type.information = "observed", backtransform = TRUE)
    name.coefM <- grep(rownames(test),pattern="M",value=TRUE)
    name.coefF <- grep(rownames(test),pattern="F",value=TRUE)
    expect_equal(test[name.coefM,"df"], rep(sum(d$Gender=="M")-length(name.coefM)+1,length(name.coefM)), tol = 1e-6)
    expect_equal(test[name.coefF,"df"], rep(sum(d$Gender=="F")-length(name.coefM)+1,length(name.coefM)), tol = 1e-6)

    ## ** confidence interval
    capture.output(summary(e.lmm))
    capture.output(anova(e.lmm))
    ## anova(e.gls,type="marginal")

    ## ** predictions
    ## GS <- lapply(AICcmodavg::predictSE(e.gls2, newdata = d),head)
    GS <- list(fit = c("1" = -0.58075328, "2" = -0.10811341, "3" = 0.90634198, "4" = 2.98227219, "5" = 1.41661756, "6" = 2.27857674) ,
     se.fit = c("1" = 0.30371618, "2" = 0.26101415, "3" = 0.37908306, "4" = 0.43713441, "5" = 0.37370611, "6" = 0.35737563) )

    pp1 <- predict(e.lmm, newdata = d, se = TRUE)
    pp2 <- predict(e.lmm2, newdata = d, se = TRUE)
    set.seed(10)
    index <- sample.int(NROW(d))
    pp3 <- predict(e.lmm, newdata = d[index,,drop=FALSE], se = TRUE)
    expect_equivalent(head(pp1$estimate),GS$fit, tol = 1e-5)
    expect_equivalent(head(pp1$se),GS$se.fit, tol = 1e-5)
    expect_equal(pp1,pp2, tol = 1e-5)
    expect_equivalent(pp1[index,,drop=FALSE],pp3, tol = 1e-5) ## different rownames
    ## .getUVarCov(e.lmm)
})

## * missing values
test_that("missing values",{
    set.seed(11)
    d$Ymiss <- d$Y
    d$Ymiss[which(rbinom(NROW(d), size = 1, prob = 0.1)==1)] <- NA
    
    ## ** fit
    e.lmm <- suppressWarnings(lmm(Ymiss ~ X1 + X2 + Gene, repetition = ~time|id,
                                  structure = "ID", data = d, trace = 0,
                                  method = "ML", df = TRUE))
    e.gls <- gls(Ymiss ~ X1 + X2 + Gene, data = d, method = "ML", na.action = na.omit)
    expect_equal(logLik(e.lmm), as.double(logLik(e.gls)), tol = 1e-6)
})


##----------------------------------------------------------------------
### test-auto-linear-regression.R ends here
