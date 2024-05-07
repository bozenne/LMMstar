### test-auto-residuals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  4 2021 (11:49) 
## Version: 
## Last-Updated: maj  7 2024 (13:58) 
##           By: Brice Ozenne
##     Update #: 37
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
    library(ggplot2)

    library(LMMstar)
}


context("Check residuals and partial residuals calculation")
LMMstar.options(optimizer = "FS", method.numDeriv = "simple", precompute.moments = TRUE,
                columns.confint = c("estimate","se","df","lower","upper","p.value"))

## * Linear model
set.seed(10)
dL <- sampleRem(100, n.times = 3, format = "long")

test_that("linear model",{

    e.lm <- lm(Y~visit+X1+X2+X6, data = dL)
    GS <- residuals(e.lm, type = "partial")

    ## single variable
    e.lmm <- lmm(Y~visit+X1+X2+X6, data = dL)
    test1 <- residuals(e.lmm, type = "partial-center", variable = "X1")
    expect_equal(as.double(test1), as.double(GS[,"X1"]), tol = 1e-6)

    ## same but with another reference
    test2 <- residuals(e.lmm, type = "partial", variable = "X1")
    e.diff <- mean(e.lmm$design$mean[,"X1"] * coef(e.lmm)["X1"]) ## expected difference
    expect_equal(as.double(test2-test1), rep(e.diff,length(test1)), tol = 1e-6)

    ## note: only match with continuous covariate
    test1.bis <- residuals(e.lmm, type = "partial-center", variable = "visit")
    e.diff.bis <- mean(e.lmm$design$mean[,c("visit2","visit3"),drop=FALSE] %*% coef(e.lmm)[c("visit2","visit3")]) ## expected difference
    expect_equal(as.double(test1.bis - GS[,"visit"]), rep(e.diff.bis,length(test1)), tol = 1e-6)
    ## plot(e.lmm, type = "partial", variable = "X1")
    ## plot(e.lmm, type = "partial", variable = c("(Intercept)","X1"))
})

## ** Linear model with interaction
test_that("linear model with interaction",{

    e.lmm <- lmm(Y~visit*X6+X2+X5, data = dL)
    ## plot(e.lmm, type = "partial", variable = c("visit","X6"))
    ## plot(e.lmm, type = "partial", variable = c("visit","X6","(Intercept)"))
    
    test1 <- residuals(e.lmm, type = "partial", variable = c("visit","X6"))
    GS1 <- dL$Y - coef(e.lmm)["(Intercept)"] - coef(e.lmm)["X2"]*dL$X2 - coef(e.lmm)["X5"]*dL$X5
    expect_equal(as.double(test1), as.double(GS1), tol = 1e-6)

    test2 <- residuals(e.lmm, type = "partial", variable = c("visit","X6","(Intercept)"))
    GS2 <- dL$Y - coef(e.lmm)["X2"]*dL$X2 - coef(e.lmm)["X5"]*dL$X5
    expect_equal(as.double(test2), as.double(GS2), tol = 1e-6)

})

## ** Linear model with splines

test_that("linear model with splines",{

    ## 1- poly
    ePOLY.lm <- lm(Y~visit+X1+stats::poly(X6,4), data = dL)
    ePOLY.lmm <- lmm(Y~visit+X1+stats::poly(X6, 4), data = dL)
    ## plot(ePOLY.lmm, type = "partial", variable = "X6")

    ## compare predictions
    GS.POLY <- predict(ePOLY.lm, newdata = dL[1:2,])
    test.POLY <- predict(ePOLY.lmm, newdata = dL[1:2,])
    expect_equal(as.double(test.POLY), as.double(GS.POLY), tol = 1e-6)
    
    ## compare partial residuals
    test.POLY <- residuals(ePOLY.lmm, type = "partial", variable = "X6")
    GS.POLY <- dL$Y - coef(ePOLY.lm)["(Intercept)"] - dL$X1 * coef(ePOLY.lm)["X1"] - c(0,coef(ePOLY.lm)[c("visit2","visit3")])[as.numeric(dL$visit)]
    expect_equal(as.double(GS.POLY), as.double(test.POLY), tol = 1e-6)

    ## ## 2- ns
    eNS.lm <- lm(Y~visit+X1+splines::ns(X6,4), data = dL)
    eNS.lmm <- lmm(Y~visit+X1+splines::ns(X6, 4), data = dL)
    ## plot(eNS.lmm, type = "partial", variable = "X6")

    ## ## compare predictions
    GS.NS <- predict(eNS.lm, newdata = dL[1:2,])
    test.NS <- predict(eNS.lmm, newdata = dL[1:2,])
    expect_equal(as.double(test.NS), as.double(GS.NS), tol = 1e-6)
    
    ## ## compare partial residuals
    test.NS <- residuals(eNS.lmm, type = "partial", variable = "X6")
    GS.NS <- dL$Y - coef(eNS.lm)["(Intercept)"] - dL$X1 * coef(eNS.lm)["X1"] - c(0,coef(eNS.lm)[c("visit2","visit3")])[as.numeric(dL$visit)]
    expect_equal(as.double(GS.NS), as.double(test.NS), tol = 1e-6)
})


## * Linear mixed model
test_that("linear mixed model",{

    ## without missing values
    e.lmm <- lmm(Y~visit+X1+X2+X5, repetition = ~visit|id, data = dL)
    pRes <- residuals(e.lmm, type = "partial", variable = c("(Intercept)","visit","X1"), keep.data = TRUE)
    GS <- coef(e.lmm)
    test <- coef(lmm(r.partial ~ visit + X1, repetition = ~visit|id, data = pRes))
    expect_equal(GS[names(test)],test, tol = 1e-4)

    ## only match lm when full interaction with time due to unstructure covariance matrix
    expect_equal(coef(lm(r.partial ~ visit + X1, data = pRes))[c("visit2","visit3")],
                 test[c("visit2","visit3")], tol = 1e-3)

    ## with missing values (only approximate equality)
    dLNA <- dL
    dLNA$Y[c(79, 13, 191, 139, 75, 213, 78, 4, 246, 274)] <- NA

    e.lmmNA <- lmm(Y~visit:X1+X2+X5, repetition = ~visit|id, data = dLNA)
    pResNA <- residuals(e.lmmNA, type = "partial", variable = c("(Intercept)","visit","X1"), keep.data = TRUE)
    GS <- coef(e.lmmNA)
    test <- coef(lmm(r.partial ~ visit:X1, repetition = ~visit|id, data = pResNA))
    expect_equal(GS[names(test)],test, tol = 1e-2)

    ## pResNA2 <- fitted(e.lmmNA, type = "outcome", keep.data = TRUE)
    ## pResNA2$r.partial <- pResNA2$Y - coef(e.lmmNA)["X2"]*pResNA2$X2 - coef(e.lmmNA)["X5"]*pResNA2$X5
    ## test2 <- coef(lmm(r.partial ~ visit:X1, repetition = ~visit|id, data = pResNA2))
    ## expect_equal(GS[names(test2)],test2, tol = 1e-3)
})

## * Linear mixed model (gradient)
data(gastricbypassL, package = "LMMstar")
e.lmm <- lmm(glucagonAUC ~ weight + visit + (1|id), data = gastricbypassL)

test_that("gradient of the residuals in lmm",{

    test <- residuals(e.lmm, type = c("response","pearson","studentized","normalized","normastudentized"), simplify = FALSE, keep.data = TRUE, fitted.ci = TRUE) 
    GS <- estimate(e.lmm, transform.sigma = "log", transform.rho = "atanh", function(p){ ## p <- NULL
        iRes <- residuals(e.lmm, type = c("response","pearson","studentized","normalized","normastudentized"), p = p, keep.data = TRUE)
        c(iRes$fitted,iRes$r.response,iRes$r.pearson,iRes$r.studentized,iRes$r.normalized,iRes$r.normastudentized)
    }, method.numDeriv = "Richardson")
    expect_equivalent(attr(GS,"grad")[1:80,], attr(test,"grad")[,,"fitted"], tol = 1e-6)
    expect_equivalent(attr(GS,"grad")[81:160,], attr(test,"grad")[,,"r.response"], tol = 1e-6)
    expect_equivalent(attr(GS,"grad")[161:240,], attr(test,"grad")[,,"r.pearson"], tol = 1e-6)
    expect_equivalent(attr(GS,"grad")[241:320,], attr(test,"grad")[,,"r.studentized"], tol = 1e-6)
    expect_equivalent(attr(GS,"grad")[321:400,], attr(test,"grad")[,,"r.normalized"], tol = 1e-6)
    expect_equivalent(attr(GS,"grad")[401:480,], attr(test,"grad")[,,"r.normastudentized"], tol = 1e-6)

    test.se <- sqrt(diag(attr(test,"grad")[,,"r.normalized"] %*% vcov(e.lmm, effects = "all") %*% t(attr(test,"grad")[,,"r.normalized"])))
    expect_equivalent(GS$se[321:400], test.se, tol = 1e-6)

    test2 <- residuals(e.lmm, variable = "weight", type = "partial", simplify = FALSE, keep.data = TRUE, fitted.ci = TRUE)
    grad.NA <- attr(test2,"grad")[is.na(gastricbypassL$glucagonAUC),,"r.partial"]
    grad.NNA <- attr(test2,"grad")[!is.na(gastricbypassL$glucagonAUC),,"r.partial"]
    expect_true(all(is.na(grad.NA)))
    expect_true(all(grad.NNA[,"(Intercept)"]==-1))
    expect_true(all(grad.NNA[,"visit2"]==-(gastricbypassL$visit[!is.na(gastricbypassL$glucagonAUC)]=="2")))
    expect_true(all(grad.NNA[,"visit3"]==-(gastricbypassL$visit[!is.na(gastricbypassL$glucagonAUC)]=="3")))
    expect_true(all(grad.NNA[,"visit4"]==-(gastricbypassL$visit[!is.na(gastricbypassL$glucagonAUC)]=="4")))
    expect_true(all(grad.NNA[,c("weight","sigma","rho(id)")]==0))
    
})
##----------------------------------------------------------------------
### test-auto-residuals.R ends here
