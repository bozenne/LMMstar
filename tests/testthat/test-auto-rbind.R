### test-auto-rbind.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:20) 
## Version: 
## Last-Updated: jul 17 2025 (17:39) 
##           By: Brice Ozenne
##     Update #: 69
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

context("Check rbind(lmm)")
LMMstar.options(optimizer = "FS", method.numDeriv = "simple", precompute.moments = TRUE,
                columns.confint = c("estimate","se","df","lower","upper","p.value"))


## * Adjusting for multiple comparisons over different models
## simulated data
set.seed(10)
dL <- sampleRem(1e2, n.times = 3, format = "long")
dL$Xcat <- as.character(rbinom(NROW(dL),size=3,prob = .5))

test_that("rbind for anova",{
    e.lmm1 <- lmm(Y ~ X1+X2+X3, repetition = ~visit|id, data = dL)
    e.lmm2 <- lmm(Y ~ X1+X8+X9, repetition = ~visit|id, data = dL)
    e.lmm3 <- lmm(Y ~ X1+Xcat, repetition = ~visit|id, data = dL)

    AAA <- anova(e.lmm1, effect = c("X1|X2,X3"="X1=0","X2|X1,X3"="X2=0"), robust = TRUE, df = FALSE, simplify = FALSE)
    BBB <- anova(e.lmm2, effect = c("X1|X8,X9"="X1=0"), robust = TRUE, df = FALSE, simplify = FALSE)
    CCC <- anova(e.lmm3, effect = c("X1|Xcat"="X1=0"), robust = TRUE, df = FALSE, simplify = FALSE)
    ZZZ <- rbind(AAA,BBB,CCC)
    test1 <- rbind(confint(AAA, method = "none"),
                   confint(BBB, method = "none"),
                   confint(CCC, method = "none"))
    test2 <- confint(ZZZ, method = "none")
    test <- capture.output(summary(ZZZ))

    expect_equal(as.double(unlist(test1[,c("estimate","se","df","lower","upper","p.value")])),
                 as.double(unlist(test2[,c("estimate","se","df","lower","upper","p.value")])), tol = 1e-6)
})



##----------------------------------------------------------------------
### test-auto-rbind.R ends here
