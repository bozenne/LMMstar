### profile.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 20 2022 (14:45) 
## Version: 
## Last-Updated: jun 20 2022 (17:52) 
##           By: Brice Ozenne
##     Update #: 6
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

    e.profile <- lm(weight2 ~ 0 + weight1 + glucagonAUC1, data = gastricbypassW)
    e.profile2 <- lmm(weight2 ~ 0 + weight1 + glucagonAUC1, data = gastricbypassW, method.fit = "ML", control = list(optimizer = "FS"))
    expect_equal(logLik(e.profile2), as.double(logLik(e.profile)), tol = 1e-5)


    profile(e.full2, profile.likelihood = TRUE)

    e.full2$design$precompute.XY
    e.full2$design$precompute.XX

    cMLE <- coef(e.profile2, effects = "all")

    test5 <- .constrain.lmm(e.full2, effects = stats::setNames(5, "(Intercept)"), init = MLE, trace = 5)
    testBest <- .constrain.lmm(e.full2, effects = stats::setNames(4.6655966323, "(Intercept)"), init = MLE, trace = 5)
    test4 <- .constrain.lmm(e.full2, effects = stats::setNames(4, "(Intercept)"), init = MLE, trace = 5)

    e.off <- lm(weight2 ~ 0 + weight1 + glucagonAUC1, data = cbind(gastricbypassW, off = 4), offset = off)
    score(e.full2, p = test4$param, effects = "all")

    X <- cbind(gastricbypassW$weight1,gastricbypassW$glucagonAUC1)
    solve(t(X) %*% X) %*% t(X) %*% cbind(gastricbypassW$weight2-4)
    solve(t(X) %*% X) %*% t(X) %*% cbind(gastricbypassW$weight2-4)
    
    newcoef <- c("(Intercept)" = 4, coef(e.off), sigma = sigma(e.off))
    rbind(value = newcoef,
          score = score(e.full2, p = newcoef, effects = "all"))
    

    test3 <- .constrain.lmm(e.full2, effects = stats::setNames(3, "(Intercept)"), init = MLE, trace = 5)
    test2 <- .constrain.lmm(e.full2, effects = stats::setNames(2, "(Intercept)"), init = MLE, trace = 5)
    test1 <- .constrain.lmm(e.full2, effects = stats::setNames(1, "(Intercept)"), init = MLE, trace = 5)

    logLik(e.full2, p = c("(Intercept)" = 0, MLE[-1]))
    SSS <- score(e.full2, p = c("(Intercept)" = 0, MLE[-1]), effects = "all")
    III <- information(e.full2, p = c("(Intercept)" = 0, MLE[-1]), effects = "all", type.information = "expected")
    MLE[4]*exp(SSS[4]/III[4,4])
    (1-MLE[4])*exp(SSS[4]/III[4,4])

    XX <- design$mean[,2:3]
    solve(crossprod(XX)) %*% t(XX) %*% design$Y
    ## Initialization:
    ##   (Intercept)       weight1  glucagonAUC1         sigma 
    ##  0.0000000000  0.9236808578 -0.0003247549  2.3487164956 

    ## Loop:
    ## iteration 1: logLik=-84.91579652
    ##          (Intercept)       weight1  glucagonAUC1    sigma
    ## estimate           0  9.236809e-01 -3.247549e-04 16.89219
    ## diff               0 -3.330669e-16 -1.897354e-18 14.54347
    ## score             NA            NA            NA 78.91929


    logLik(e.full2, p = test$param)
    score(e.full2, p = test$param)
    score(e.full2, p = cMLE)
    ## expect_equal(logLik(e.full2, p = cMLE), as.double(logLik(e.profile)), tol = 1e-5)
    
    e.full2 <- lmm(weight2 ~ weight1 + glucagonAUC1, data = gastricbypassW, method.fit = "ML")
    expect_equal(logLik(e.full2), as.double(logLik(e.full)), tol = 1e-5)
})

##----------------------------------------------------------------------
### profile.R ends here
