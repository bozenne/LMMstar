### test-auto-anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul 13 2022 (13:55) 
## Version: 
## Last-Updated: jul 24 2025 (16:54) 
##           By: Brice Ozenne
##     Update #: 72
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
    library(testthat)
    library(mice)
    library(Matrix)

    library(LMMstar)
}

context("Check LRT and Wald tests")
LMMstar.options(optimizer = "FS", precompute.moments = TRUE,
                columns.confint = c("estimate","se","df","lower","upper","p.value"))

## * Simulate data
set.seed(10)
dL <- sampleRem(100, n.times = 3, format = "long")

dL$X1 <- as.factor(dL$X1)
dL$X2 <- as.factor(dL$X2)

## * Likelihood ratio tests
test_that("LRT", {

    ## remove variance factor
    e0 <- anova(lmm(Y ~ X1 + X2, data = dL),
                lmm(Y ~ X1 + X2, repetition = ~visit, structure = "IND", data = dL))
    expect_equal(list(gsub(" ","",e0$null), e0$df, e0$statistic, e0$p.value),
                 list("k.2==1\nk.3==1", 2, 0.1421563, 0.9313891), tol = 1e-4)

    dL$id2 <- 1:NROW(dL)
    dL$time2 <- 1
    e00 <- anova(lmm(Y ~ X1 + X2, data = dL),
                 lmm(Y ~ X1 + X2, repetition = ~time2|id2, structure = ID(visit~1), data = dL, control = list(optimizer = "FS")))
    expect_equal(list(gsub(" ","",e00$null), e00$df, e00$statistic, e00$p.value),
                 list("sigma:2==sigma:3==sigma:1", 2, 0.1421563, 0.9313891), tol = 1e-4)

    ## remove mean factor
    e1 <- anova(lmm(Y ~ X1 + X2, repetition = ~visit|id, structure = "CS", data = dL, method.fit = "ML"),
                lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "CS", data = dL, method.fit = "ML"))
    expect_equal(list(gsub(" ","",e1$null), e1$df, e1$statistic, e1$p.value),
                 list("X5==0", 1, 0.4686210, 0.4936222), tol = 1e-4)

    ## swap mean factor
    expect_error(anova(lmm(Y ~ X1 + X2, repetition = ~visit|id, structure = "CS", data = dL, method.fit = "ML"),
                       lmm(Y ~ X1 + X3, repetition = ~visit|id, structure = "UN", data = dL, method.fit = "ML")))

    ## remove correlation factor
    ## via structure
    e2 <- anova(lmm(Y ~ X1:X2 + X5, repetition = ~visit|id, structure = "CS", data = dL, method.fit = "ML"),
                lmm(Y ~ X1*X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, method.fit = "ML"))

    expect_equal(list(gsub(" ","",e2$null), e2$df, e2$statistic, e2$p.value),
                 list("k.2==1\nk.3==1\nrho(1,2)==rho(1,3)==rho(2,3)", 4, 29.66939, 5.714571e-06), tol = 1e-4)

    ## via strata
    e3 <- anova(lmm(Y ~ X1 + X5, repetition = ~visit|id, structure = "UN", data = dL, method.fit = "ML", control = list(optimizer = "FS")),
                lmm(Y ~ X1 + X5, repetition = X2~visit|id, structure = "UN", data = dL, method.fit = "ML", control = list(optimizer = "FS")))

    expect_equal(list(gsub(" ","",e3$null), e3$df, e3$statistic, e3$p.value),
                 list("sigma:1==sigma:0\nk.2:1==k.2:0\nk.3:1==k.3:0\nrho(1,2):1==rho(1,2):0\nrho(1,3):1==rho(1,3):0\nrho(2,3):1==rho(2,3):0", 6, 0.5038597, 0.9977912), tol = 1e-4)

    ## both (does not work for now)
    ## anova(lmm(Y ~ X1 + X5, repetition = ~visit|id, structure = "CS", data = dL, method.fit = "ML", control = list(optimizer = "FS")),
    ##       lmm(Y ~ X1 + X5, repetition = X2~visit|id, structure = "UN", data = dL, method.fit = "ML", control = list(optimizer = "FS")))

})

## * Wald test (single model)
e.lmm1 <- lmm(Y ~ X1+X2+X3, repetition = ~visit|id, data = dL,
              structure = "UN", df = FALSE)

dL2 <- dL
dL2$id <- paste0(dL$id,".bis")
dL2$Y2 <- dL$Y
e.lmm2 <- lmm(Y2 ~ X1+X2+X3, repetition = ~visit|id, data = dL2,
              structure = "UN", df = FALSE)

## ** back-transform
test_that("anova_lmm vs. confint", {

    ## anova(e.lmm1, effects = "X11=0")

    ## check back-transform is working like confint
    GS <- model.tables(e.lmm1, effects = "all")
    test1 <- model.tables(anova(e.lmm1, effects = "all"), method = "none")
    expect_equivalent(test1, GS[sapply(strsplit(rownames(test1), split = "="),"[[",1),], tol = 1e-5)

    test2 <- model.tables(anova(e.lmm1, effects = c("k.2=1","X11=0")), method = "none")
    expect_equivalent(test2, GS[c("k.2","X11"),], tol = 1e-5)

    ## default no back
    test1.bis <- model.tables(anova(e.lmm1, effects = "all"), backtransform = FALSE, method = "none")
    expect_equal(test1$p.value, test1.bis$p.value, tol = 1e-5)
    expect_equal(confint(e.lmm1, backtransform = FALSE, effects = "all", transform.names = FALSE)[sapply(strsplit(rownames(test1), split = "="),"[[",1),"estimate"],
                 unname(test1.bis$estimate), tol = 1e-5)

})

## ** pooling
e.lmm3 <- lmm(Y ~ X1*X2+X5, repetition = ~visit|id, data = dL,
              structure = "UN", df = TRUE)

test_that("pooling anova_lmm ", {

    e.anova <- anova(e.lmm3, effects = c("X11=0","X21=0","X5=0","X11:X21=0"), simplify = FALSE)

    ## average
    test.average <- model.tables(e.anova, method = "average")
    test.average2 <- estimate(e.lmm3, function(p){mean(p[2:5])})
    GS.average <- anova(e.lmm3, effects = c("average" = "0.25*X11 + 0.25*X21 + 0.25*X5 + 0.25*X11:X21=0"))
    expect_equal(as.double(model.tables(GS.average)), as.double(test.average), tol = 1e-6)
    expect_equal(as.double(model.tables(GS.average)[,c("estimate","se","df")]),
                 as.double(test.average2[,c("estimate","se","df")]), tol = 1e-4)

    ## pool.se
    test.se <- model.tables(e.anova, method = "pool.se")
    test.se2 <- estimate(e.lmm3, function(p){
        weighted.mean(p[2:5], 1/diag(vcov(e.lmm3, p = p))[2:5])
    })
    expect_equivalent(c(0.163520459587233, 0.0132631732100504, 0.816759191973988, 0.00645717522872792),
                      weights(e.anova, method = "pool.se")[1,], tol = 1e-6) ## diag(1/vcov(e.lmm3))[-1]/sum(1/diag(vcov(e.lmm3))[-1])
    ## pool.fixse (no uncertainty about the weights)
    GS.fixse <- anova(e.lmm3, effects = c("pool.se" = "0.163520470*X11 + 0.013263174*X21 + 0.816759181*X5 + 0.006457176*X11:X21=0"))

    ## no uncertainty about the weights
    expect_equal(test.se2$estimate, test.se$estimate, tol = 1e-6) ## same
    expect_equal(test.se2$se, test.se$se, tol = 1e-5) ## small difference (1st order approximation of the delta method)
    expect_equal(test.se2$df, test.se$df, tol = 1) ## not too big difference (1st order approximation of the delta method)
    
    expect_equal(model.tables(GS.fixse)$estimate, test.se$estimate, tol = 1e-6)
    expect_true(test.se$se > model.tables(GS.fixse)$se) ## too low variance when ignoring the variability of the weights (small difference in practice)
  
    ## pool.gls
    test.gls <- model.tables(e.anova, method = "pool.gls")    
    expect_equal(as.double(test.gls), c(0.24327505, 0.25866845, 90.40240892, -0.2705838, 0.7571339, 0.34947263), tol = 1e-3)
    ## estimate(e.anova, function(object){coef(object, method = "pool.gls")})

    test.rubin <- model.tables(e.anova, method = "pool.rubin")
    expect_equal(test.rubin$estimate,test.average$estimate, tol = 1e-6)
})


## * Wald test (multiple models via rbind)

test_that("rbind.anova_lmm", {

    ## same clusters - model-based
    A123.bis <- rbind(anova(e.lmm1, effect = c("X11=0"), simplify = FALSE),
                      anova(e.lmm1, effect = c("X21=0"), simplify = FALSE),
                      anova(e.lmm1, effect = c("X3=0"), simplify = FALSE))
    A123 <- anova(e.lmm1, effect = c("X11=0","X21=0","X3=0"))

    expect_equal(coef(A123), coef(A123.bis), tol = 1e-5)
    expect_equal(diag(vcov(A123)), diag(vcov(A123.bis)), tol = 1e-5)
    ## extra-diagonal elements differ as some are model-based and other are robust-based

    ## same clusters - robust
    A123.bis_robust <- rbind(anova(e.lmm1, effect = c("X11=0"), robust = TRUE, simplify = FALSE),
                             anova(e.lmm1, effect = c("X21=0"), robust = TRUE, simplify = FALSE),
                             anova(e.lmm1, effect = c("X3=0"), robust = TRUE, simplify = FALSE))
    A123_robust <- anova(e.lmm1, effect = c("X11=0","X21=0","X3=0"), robust = TRUE)

    expect_equal(coef(A123_robust), coef(A123.bis_robust), tol = 1e-5)
    expect_equal(vcov(A123_robust), vcov(A123.bis_robust), tol = 1e-5)

    ## different clusters
    A1 <- anova(e.lmm1, effect = c("X11=0","X21=0","X3=0"), simplify = FALSE)
    A2 <- anova(e.lmm2, effect = c("X11=0","X21=0","X3=0"), simplify = FALSE)
    A12.ter <- rbind(A1,A2)
    A123 <- anova(e.lmm1, effect = c("X11=0","X21=0","X3=0"))

    expect_equal(unname(rep(coef(A123),2)), unname(coef(A12.ter)), tol = 1e-5)
    expect_equal(unname(as.matrix(bdiag(vcov(A123),vcov(A123)))), unname(vcov(A12.ter)), tol = 1e-5)
    expect_equal(A12.ter$multivariate$statistic, A123$multivariate$statistic, tol = 1e-5)

})

## * Rubin's rule
set.seed(123)
df.NA <- mice(nhanes, printFlag = FALSE)

test_that("Rubin's rule", {

    ## mice
    ls.mice <- with(data=df.NA,exp=lm(chl~bmi))
    GS <- summary(pool(ls.mice))

    ## lmm
    df.NNA <- complete(df.NA, action = "long")
    e.mlmm <- mlmm(chl~bmi, by = ".imp", repetition = ~1|.id, data = df.NNA, effects = "bmi=0", trace = FALSE)
    ## coef(e.mlmm)
    ## vcov(e.mlmm)
    test <- model.tables(e.mlmm, method = "pool.rubin")

    expect_equal(as.double(GS[GS$term=="bmi",c("estimate","std.error")]),
                 as.double(test[,c("estimate","se")]), tol = 1e-6)
    expect_equal(as.double(GS[GS$term=="bmi",c("df","p.value")]),
                 as.double(test[,c("df","p.value")]), tol = 1e-2)

})

##----------------------------------------------------------------------
### test-auto-anova.R ends here
