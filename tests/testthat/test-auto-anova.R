### test-auto-anova.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul 13 2022 (13:55) 
## Version: 
## Last-Updated: maj  7 2024 (10:27) 
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
    library(lava)
    library(reshape2)
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
    expect_equal(list(e0[,c("null")], e0[,c("df")], e0[,c("statistic")], e0[,c("p.value")]),
                 list("k.2==0, k.3==0", 2, 0.1421563, 0.9313891), tol = 1e-4)

    dL$id2 <- 1:NROW(dL)
    dL$time2 <- 1
    e00 <- anova(lmm(Y ~ X1 + X2, data = dL),
                 lmm(Y ~ X1 + X2, repetition = ~time2|id2, structure = ID(visit~1), data = dL, control = list(optimizer = "FS")))
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
                 list("k.2==0, k.3==0, rho(1,2)==rho(id), rho(1,3)==rho(id), rho(2,3)==rho(id)", 4, 29.66939, 5.714571e-06), tol = 1e-4)

    ## via strata
    e3 <- anova(lmm(Y ~ X1 + X5, repetition = ~visit|id, structure = "UN", data = dL, method.fit = "ML", control = list(optimizer = "FS")),
                lmm(Y ~ X1 + X5, repetition = X2~visit|id, structure = "UN", data = dL, method.fit = "ML", control = list(optimizer = "FS")))

    expect_equal(list(e3[,c("null")], e3[,c("df")], e3[,c("statistic")], e3[,c("p.value")]),
                 list("sigma:0==sigma, sigma:1==sigma, k.2:0==k.2, k.3:0==k.3, k.2:1==k.2, k.3:1==k.3, rho(1,2):0==rho(1,2), rho(1,3):0==rho(1,3), rho(2,3):0==rho(2,3), rho(1,2):1==rho(1,2), rho(1,3):1==rho(1,3), rho(2,3):1==rho(2,3)", 6, 0.5038597, 0.9977912), tol = 1e-4)

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
    expect_equal(test1, GS[rownames(test1),], tol = 1e-5)

    test2 <- model.tables(anova(e.lmm1, effects = c("k.2=0","X11=0")), method = "none")
    expect_equal(test2, GS[rownames(test2),], tol = 1e-5)

    ## default no back
    test1.bis <- model.tables(anova(e.lmm1, effects = "all"), backtransform = FALSE, method = "none")
    expect_equal(test1$p.value, test1.bis$p.value, tol = 1e-5)
    expect_equal(confint(e.lmm1, backtransform = FALSE, effects = "all", transform.names = FALSE)[rownames(test1.bis),"estimate"],
                 unname(test1.bis$estimate), tol = 1e-5)

})

## ** pooling
e.lmm3 <- lmm(Y ~ X1*X2+X5, repetition = ~visit|id, data = dL,
              structure = "UN", df = TRUE)

test_that("pooling anova_lmm ", {

    e.anova <- anova(e.lmm3, effects = c("X11=0","X21=0","X5=0","X11:X21=0"))

    ## average
    test.average <- model.tables(e.anova, method = "average")
    test.average2 <- estimate(e.lmm3, function(p){mean(p[2:5])})
    GS.average <- anova(e.lmm3, effects = c("average" = "0.25*X11 + 0.25*X21 + 0.25*X5 + 0.25*X11:X21=0"))
    expect_equal(as.double(model.tables(GS.average)), as.double(test.average), tol = 1e-6)
    expect_equal(as.double(model.tables(GS.average)[,c("estimate","se","df")]),
                 as.double(test.average2[,c("estimate","se","df")]), tol = 1e-4)

    ## pool.se
    test.se <- model.tables(e.anova, method = "pool.fixse")
    test.se2 <- estimate(e.lmm3, function(p){
        weighted.mean(p[2:5], 1/diag(vcov(e.lmm3))[2:5])
    })
    GS.se <- anova(e.lmm3, effects = c("pool.se" = "0.163520470*X11 + 0.013263174*X21 + 0.816759181*X5 + 0.006457176*X11:X21=0"))
    ## diag(1/vcov(e.lmm3))[-1]/sum(1/diag(vcov(e.lmm3))[-1])

    ## no uncertainty about the weights
    expect_equal(as.double(model.tables(GS.se)), as.double(test.se), tol = 1e-6)
    expect_equal(as.double(model.tables(GS.se)$estimate), as.double(test.se2$estimate), tol = 1e-6)
    expect_equal(as.double(test.se2), c(0.2862191, 0.2762587, 90.2378692, -0.2625972, 0.8350354, 0.3029450), tol = 1e-3)
    
    ## pool.gls
    test.gls <- model.tables(e.anova, method = "pool.gls")    
    expect_equal(as.double(test.gls), c(0.2432747, 0.2586498, 90.4081859, -0.2705465, 0.7570960,  0.3494384), tol = 1e-3)

    test.rubin <- model.tables(e.anova, method = "pool.rubin")
    expect_equal(test.rubin$estimate,test.average$estimate, tol = 1e-6)

})


## * Wald test (multiple models via rbind)

test_that("rbind.anova_lmm", {

    ## same clusters
    A1 <- anova(e.lmm1, ci = TRUE, effect = c("X11=0"))
    A2 <- anova(e.lmm1, ci = TRUE, effect = c("X21=0"))
    A3 <- anova(e.lmm1, ci = TRUE, effect = c("X3=0"))
    A123.bis <- rbind(A1,A2,A3)
    A123 <- anova(e.lmm1, ci = TRUE, effect = c("X11=0","X21=0","X3=0"), robust = TRUE)

    expect_equal(coef(A123), coef(A123.bis), tol = 1e-5)
    expect_equal(vcov(A123), vcov(A123.bis), tol = 1e-5)

    ## sqrt(diag(vcov(A123)))
    ## sqrt(diag(vcov(A123.bis)))
    ## A123.bis$univariate ## rbind.Wald_lmm(A1,A2,A3)
    GS <- model.tables(A123)
    test <- model.tables(A123.bis)
    expect_equal(A123.bis$univariate[,c("estimate","se","statistic")],
                 A123$univariate[,c("estimate","se","statistic")],
                 tol = 1e-5)

    ## different clusters
    A1 <- anova(e.lmm1, ci = TRUE, effect = c("X11=0","X21=0","X3=0"))
    A2 <- anova(e.lmm2, ci = TRUE, effect = c("X11=0","X21=0","X3=0"))
    A12.ter <- rbind(A1,A2)
    A123 <- anova(e.lmm1, ci = TRUE, effect = c("X11=0","X21=0","X3=0"))

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
    e.lmm <- mlmm(chl~bmi, by = ".imp", data = df.NNA,
                  effects = "bmi=0", trace = FALSE)
    test <- model.tables(e.lmm, method = "pool.rubin")

    expect_equal(as.double(GS[GS$term=="bmi",c("estimate","std.error")]),
                 as.double(test[,c("estimate","se")]), tol = 1e-6)
    expect_equal(as.double(GS[GS$term=="bmi",c("df","p.value")]),
                 as.double(test[,c("df","p.value")]), tol = 1e-2)

})

##----------------------------------------------------------------------
### test-auto-anova.R ends here
