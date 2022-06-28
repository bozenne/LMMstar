### profile.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 20 2022 (14:45) 
## Version: 
## Last-Updated: jun 28 2022 (11:47) 
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
}

context("Profile likelihood")
LMMstar.options(df = FALSE, optimizer = "FS")

## * lm
data(gastricbypassW, package = "LMMstar")

test_that("Profile likelihood (lm)",{

    e.full <- lm(weight2 ~ weight1 + glucagonAUC1, data = gastricbypassW)
    e.full2 <- lmm(weight2 ~ weight1 + glucagonAUC1, data = gastricbypassW, method.fit = "ML", control = list(optimizer = "FS"))
    expect_equal(logLik(e.full2), as.double(logLik(e.full)), tol = 1e-5)
    MLE <- coef(e.full2, effects = "all")

    e0.profile <- lm(weight2 ~ 0 + weight1 + glucagonAUC1, data = gastricbypassW)
    e0.profile2 <- lmm(weight2 ~ 0 + weight1 + glucagonAUC1, data = gastricbypassW, method.fit = "ML", control = list(optimizer = "FS"))
    expect_equal(logLik(e0.profile2), as.double(logLik(e0.profile)), tol = 1e-5)
    cMLE <- coef(e0.profile2, effects = "all")

    e4.profile <- lm(weight2 ~ 0 + weight1 + glucagonAUC1, data = cbind(gastricbypassW, off = 4), offset = off)
    test <- profile(e.full2, effects = "(Intercept)", profile.likelihood = TRUE, plot = FALSE, maxpts = c(0,4))
    expect_equal(test$logLik[1:2], as.double(c(logLik(e0.profile),logLik(e4.profile))), tol = 1e-5)

})

##----------------------------------------------------------------------
### profile.R ends here
