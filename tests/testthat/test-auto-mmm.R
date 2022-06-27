### test-auto-mmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:20) 
## Version: 
## Last-Updated: jun 27 2022 (11:55) 
##           By: Brice Ozenne
##     Update #: 52
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
    library(lava)
    library(multcomp)
    library(nlme)

    library(LMMstar)
}

context("Check lmm ")
LMMstar.options(optimizer = "gls", method.numDeriv = "Richardson", precompute.moments = TRUE,
                columns.confint = c("estimate","se","df","lower","upper","p.value"))


## * Adjusting for multiple comparisons over different mixed models
## simulated data
set.seed(10)
dL <- sampleRem(1e2, n.times = 3, format = "long")
dL$Xcat <- as.character(rbinom(NROW(dL),size=3,prob = .5))

test_that("rbind for anova",{
    e.lmm1 <- lmm(Y ~ X1+X2+X3, repetition = ~visit|id, data = dL)
    e.lmm2 <- lmm(Y ~ X1+X8+X9, repetition = ~visit|id, data = dL)
    e.lmm3 <- lmm(Y ~ X1+Xcat, repetition = ~visit|id, data = dL)

    AAA <- anova(e.lmm1, ci = TRUE, effect = c("X1|X2,X3"="X1=0","X2|X1,X3"="X2=0"))
    BBB <- anova(e.lmm2, ci = TRUE, effect = c("X1|X8,X9"="X1=0"))
    CCC <- anova(e.lmm3, ci = TRUE, effect = c("X1|Xcat"="X1=0"))
    ZZZ <- rbind(AAA,BBB,CCC)
    test1 <- rbind(confint(AAA, method = "none"),
                   confint(BBB, method = "none"),
                   confint(CCC, method = "none"))
    test2 <- confint(ZZZ, method = "none")
    test <- capture.output(summary(ZZZ))

    expect_equal(as.double(unlist(test1[1:2,c("estimate","se","df","lower","upper","p.value")])), as.double(unlist(model.tables(e.lmm1)[c("X1","X2"),])), tol = 1e-6)
    expect_equal(as.double(unlist(test1[,c("estimate","se","df","lower","upper","p.value")])), as.double(unlist(test2[,c("estimate","se","df","lower","upper","p.value")])), tol = 1e-6)


    ## outplot <- plot(e.lmm2, var = "X8", type = "partial")
    ## outplot$data
    ## xxx <- residuals(e.lmm2, var = "X8", type = "partial", format = "long", keep.data = TRUE)
})





##----------------------------------------------------------------------
### test-auto-mmm.R ends here
