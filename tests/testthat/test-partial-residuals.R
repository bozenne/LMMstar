### test-partial-residuals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  4 2021 (11:49) 
## Version: 
## Last-Updated: nov  4 2021 (16:34) 
##           By: Brice Ozenne
##     Update #: 11
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
LMMstar.options(optimizer = "gls", method.numDeriv = "simple", precompute.moments = TRUE,
                columns.confint = c("estimate","se","df","lower","upper","p.value"))

## * Linear model
set.seed(10)
dL <- sampleRem(100, n.times = 3, format = "long")

test_that("linear model",{

    e.lm <- lm(Y~visit+X1+X2+X6, data = dL)
    GS <- residuals(e.lm, type = "partial")

    ## single variable
    e.lmm <- lmm(Y~visit+X1+X2+X6, data = dL)
    test1 <- residuals(e.lmm, type = "partial", var = "visit")
    test1.bis <- residuals(e.lmm, type = "partial", var = "visit", format = "wide")
    test2 <- residuals(e.lmm, type = "partial", var = "X1")
    expect_equal(as.double(test1), as.double(GS[,"visit"]), tol = 1e-6)
    expect_equal(as.double(test2), as.double(GS[,"X1"]), tol = 1e-6)

    ## same but with another reference
    test3 <- residuals(e.lmm, type = "partial-ref", var = "visit")
    e.diff <- mean(e.lmm$design$mean[,c("visit2","visit3")] %*% coef(e.lmm)[c("visit2","visit3")]) ## expected difference
    expect_equal(as.double(test3-test1), rep(e.diff,length(test1)), tol = 1e-6)

    ## plot(e.lmm, type = "partial", var = "X6")
    
    ## ## manual plot
    ## rr <- residuals(e.lmm, type = "partial-ref", var = "X6", keep.data = TRUE)
    ## dd.gg <- predict(e.lmm, newdata = rr, keep.newdata = TRUE, type = "static0")
    ## gg <- ggplot(dd.gg)
    ## gg <- gg + geom_point(aes(x = X6, y = r.partial)) + geom_line(aes(x = X6, y = estimate))
    ## gg <- gg + geom_ribbon(aes(x = X6, ymin = lower, ymax = upper), alpha = 0.4)
    ## gg
    ## gg + geom_smooth(aes(x = X6, y = r.partial), method = "lm")
    ## termplot(e.lm, terms = c("X6"), partial.resid = TRUE, se = TRUE)
})

test_that("linear model with interaction",{

    e.lm <- lm(Y~visit*X6+X2+X5, data = dL)
    GS <- residuals(e.lm, type = "partial")
    
    ## single variable
    e.lmm <- lmm(Y~visit*X6+X2+X5, data = dL)
    test1 <- residuals(e.lmm, type = "partial", var = "visit")
    test1.bis <- residuals(e.lmm, type = "partial", var = "visit", format = "wide")
    test2 <- residuals(e.lmm, type = "partial", var = "X6")
    expect_equal(as.double(test1), as.double(GS[,"visit"]), tol = 1e-6)
    expect_equal(as.double(test2), as.double(GS[,"X6"]), tol = 1e-6)

    ## multiple variables
    test3 <- residuals(e.lmm, type = "partial", var = c("visit","X6"))
    
    ## plot(e.lmm, type = "partial", var = c("visit","X6"))

    ## ## manual plot
    ## rr <- residuals(e.lmm, type = c("partial-ref"), var = c("visit","X6"), keep.data = TRUE)
    ## dd.gg <- predict(e.lmm, newdata = rr, keep.newdata = TRUE, type = "static0")
    ## gg <- ggplot(dd.gg)
    ## gg <- gg + geom_point(aes(x = X6, y = r.partial)) + geom_line(aes(x = X6, y = estimate, group = visit, color = visit))
    ## gg <- gg + geom_ribbon(aes(x = X6, ymin = lower, ymax = upper, group = visit), alpha = 0.2)
    ## gg

})

test_that("linear model with splines",{
    ## 1- poly
    e.lm <- lm(Y~visit+X1+stats::poly(X6,4), data = dL)
    e.lmm <- lmm(Y~visit+X1+stats::poly(X6, 4), data = dL)

    ## compare predictions
    GS <- predict(e.lm, newdata = dL[1:2,])
    test <- predict(e.lmm, newdata = dL[1:2,])
    expect_equal(as.double(test$estimate), as.double(GS), tol = 1e-6)
    
    ## compare partial residuals
    GS <- residuals(e.lm, type = "partial")
    test <- residuals(e.lmm, type = "partial", var = "X6")
    expect_equal(as.double(test), as.double(GS[,"stats::poly(X6, 4)"]), tol = 1e-6)

    ## plot(e.lmm, type = "partial", var = "X6")

    ## ## manual plot
    ## rr <- residuals(e.lmm, type = "partial-ref", var = "X6", keep.data = TRUE)
    ## dd.gg <- predict(e.lmm, newdata = rr, keep.newdata = TRUE, type = "static0")
    ## gg <- ggplot(dd.gg)
    ## gg <- gg + geom_point(aes(x = X6, y = r.partial)) + geom_line(aes(x = X6, y = estimate))
    ## gg <- gg + geom_ribbon(aes(x = X6, ymin = lower, ymax = upper), alpha = 0.4)
    ## gg
    ## gg + geom_smooth(aes(x = X6, y = r.partial), method = "lm")
    ## termplot(e.lm, terms = c("X6"), partial.resid = TRUE, se = TRUE)

    ## ## 2- ns
    ## e.lm <- lm(Y~visit+X1+splines::ns(X6,4), data = dL)
    ## e.lmm <- lmm(Y~visit+X1+splines::ns(X6, 4), data = dL)

    ## ## compare predictions
    ## GS <- predict(e.lm, newdata = dL[1:2,])
    ## test <- predict(e.lmm, newdata = dL[1:2,])
    ## expect_equal(as.double(test$estimate), as.double(GS), tol = 1e-6)
    
    ## ## compare partial residuals
    ## GS <- residuals(e.lm, type = "partial")
    ## test <- residuals(e.lmm, type = "partial", var = "X6")
    ## expect_equal(as.double(test), as.double(GS[,"splines::ns(X6, 4)"]), tol = 1e-6)
})

##----------------------------------------------------------------------
### test-partial-residuals.R ends here
