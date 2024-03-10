### test-auto-partial-residuals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  4 2021 (11:49) 
## Version: 
## Last-Updated: Mar 10 2024 (16:20) 
##           By: Brice Ozenne
##     Update #: 27
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


context("Check partial residuals calculation")
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

##----------------------------------------------------------------------
### test-auto-partial-residuals.R ends here
