### test-auto-mmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:20) 
## Version: 
## Last-Updated: mar 12 2024 (09:46) 
##           By: Brice Ozenne
##     Update #: 60
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
LMMstar.options(optimizer = "FS", method.numDeriv = "Richardson", precompute.moments = TRUE,
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

    AAA <- anova(e.lmm1, ci = TRUE, effect = c("X1|X2,X3"="X1=0","X2|X1,X3"="X2=0"), robust = TRUE, df = FALSE)
    BBB <- anova(e.lmm2, ci = TRUE, effect = c("X1|X8,X9"="X1=0"), robust = TRUE, df = FALSE)
    CCC <- anova(e.lmm3, ci = TRUE, effect = c("X1|Xcat"="X1=0"), robust = TRUE, df = FALSE)
    ZZZ <- rbind(AAA,BBB,CCC)
    test1 <- rbind(confint(AAA, method = "none"),
                   confint(BBB, method = "none"),
                   confint(CCC, method = "none"))
    test2 <- confint(ZZZ, method = "none")
    test <- capture.output(summary(ZZZ))

    expect_equal(as.double(unlist(test1[,c("estimate","se","df","lower","upper","p.value")])),
                 as.double(unlist(test2[,c("estimate","se","df","lower","upper","p.value")])), tol = 1e-6)


    ## outplot <- plot(e.lmm2, var = "X8", type = "partial")
    ## outplot$data
    ## xxx <- residuals(e.lmm2, var = "X8", type = "partial", format = "long", keep.data = TRUE)
})


## * Delta method over differents models
set.seed(10)
dL <- sampleRem(1e2, n.times = 1, format = "long")

test_that("estimate for rbind.anova",{

    m.lvm <- lvm(Y ~ X5 + X1,
                 X5 ~ X1)
    e.lvm <- estimate(m.lvm, data = dL)
    GS <- effects(e.lvm, Y~X1)

    e.lmm1 <- lmm(Y ~ X5 + X1, repetition = ~visit|id, data = dL)
    e.aov1 <- anova(e.lmm1, effects = c("X1=0","X5=0"))
    e.lmm2 <- lmm(X5 ~ X1, repetition = ~visit|id, data = dL)
    e.aov2 <- anova(e.lmm2, effects = "X1=0")
    e.maov <- rbind(e.aov1, e.aov2)

    test <- estimate(e.maov, f = function(p){
        DE <- as.double(p["Y: X1"])
        IE <- as.double(p["Y: X5"]*p["X5: X1"])
        out <- c(direct = DE,
                 indirect = IE,
                 total = DE + IE,
                 proportion = 1/(1 + DE/IE))
    })

    expect_equal(as.double(GS$coef)[1:3], as.double(test$estimate)[c(3,1,2)], tol = 1e-5)
})

##----------------------------------------------------------------------
### test-auto-mmm.R ends here
