### test-auto-anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul 13 2022 (13:55) 
## Version: 
## Last-Updated: Jul 14 2022 (12:38) 
##           By: Brice Ozenne
##     Update #: 4
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

if(FALSE){
    library(lava)
    library(reshape2)
    library(testthat)

    library(LMMstar)
}

context("Check anova method")

## * Likelihood ratio tests

## ** simulate data
set.seed(10)
dL <- sampleRem(100, n.times = 3, format = "long")

dL$X1 <- as.factor(dL$X1)
dL$X2 <- as.factor(dL$X2)

## ** LRT
test_that("LRT"{

    ## remove variance factor
    e0 <- suppressMessages(anova(lmm(Y ~ X1 + X2, data = dL),
                                 lmm(Y ~ X1 + X2, repetition = ~visit, structure = "IND", data = dL)))
    expect_equal(list(e0[,c("null")], e0[,c("df")], e0[,c("statistic")], e0[,c("p.value")]),
                 list("k.2==0, k.3==0", 2, 0.1421563, 0.9313891), tol = 1e-4)

    dL$id2 <- 1:NROW(dL)
    dL$time2 <- 1
    e00 <- suppressMessages(anova(lmm(Y ~ X1 + X2, data = dL),
                                  lmm(Y ~ X1 + X2, repetition = ~time2|id2, structure = IND(visit~1), data = dL, control = list(optimizer = "FS"))))
    expect_equal(list(e00[,c("null")], e00[,c("df")], e00[,c("statistic")], e00[,c("p.value")]),
                 list("sigma:1==sigma, sigma:2==sigma, sigma:3==sigma", 2, 0.1421563, 0.9313891), tol = 1e-4)

    ## remove mean factor
    e1 <- anova(lmm(Y ~ X1 + X2, repetition = ~visit|id, structure = "CS", data = dL, method.fit = "ML"),
                lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "CS", data = dL, method.fit = "ML"))

    expect_equal(list(e1[,c("null")], e1[,c("df")], e1[,c("statistic")], e1[,c("p.value")]),
                 list("X5==0", 1, 0.4686210, 0.4936222), tol = 1e-4)

    ## swap mean factor
    expect_error(anova(lmm(Y ~ X1 + X2, repetition = ~visit|id, structure = "CS", data = dL, method.fit = "ML"),
                       lmm(Y ~ X1 + X3, repetition = ~visit|id, structure = "UN", data = dL, method.fit = "ML")))

    ## remove correlation factor
    ## via structure
    e2 <- anova(lmm(Y ~ X1:X2 + X5, repetition = ~visit|id, structure = "CS", data = dL, method.fit = "ML"),
                lmm(Y ~ X1*X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, method.fit = "ML"))

    expect_equal(list(e2[,c("null")], e2[,c("df")], e2[,c("statistic")], e2[,c("p.value")]),
                 list("k.2==0, k.3==0, rho(1,2)==rho, rho(1,3)==rho, rho(2,3)==rho", 4, 29.66939, 5.714571e-06), tol = 1e-4)

    ## via strata
    e3 <- anova(lmm(Y ~ X1 + X5, repetition = ~visit|id, structure = "UN", data = dL, method.fit = "ML", control = list(optimizer = "FS")),
                lmm(Y ~ X1 + X5, repetition = X2~visit|id, structure = "UN", data = dL, method.fit = "ML", control = list(optimizer = "FS")))

    expect_equal(list(e3[,c("null")], e3[,c("df")], e3[,c("statistic")], e3[,c("p.value")]),
                 list("sigma:0==sigma, sigma:1==sigma, k.2:0==k.2, k.3:0==k.3, k.2:1==k.2, k.3:1==k.3, rho(1,2):0==rho(1,2), rho(1,3):0==rho(1,3), rho(2,3):0==rho(2,3), rho(1,2):1==rho(1,2), rho(1,3):1==rho(1,3), rho(2,3):1==rho(2,3)", 6, 0.5038597, 0.9977912), tol = 1e-4)

    ## both (does not work for now)
    ## anova(lmm(Y ~ X1 + X5, repetition = ~visit|id, structure = "CS", data = dL, method.fit = "ML", control = list(optimizer = "FS")),
    ##       lmm(Y ~ X1 + X5, repetition = X2~visit|id, structure = "UN", data = dL, method.fit = "ML", control = list(optimizer = "FS")))

})

##----------------------------------------------------------------------
### test-auto-anova.R ends here
