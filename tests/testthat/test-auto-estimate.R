### test-auto-estimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 31 2022 (11:36) 
## Version: 
## Last-Updated: maj  7 2024 (10:35) 
##           By: Brice Ozenne
##     Update #: 81
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
    library(lme4)

    library(LMMstar)
}

context("Check delta method (estimate function) for mixed model")
LMMstar.options(optimizer = "FS", method.numDeriv = "simple", precompute.moments = TRUE,
                columns.confint = c("estimate","se","df","lower","upper","p.value"))

## * Compare change with complete data
## ** Simulate data
set.seed(10)
d <- sampleRem(1e2, n.time = 2)
d$dX <- d$X2 - d$X1
d$dY <- d$Y2 - d$Y1

## ** ANCOVA
e.ANCOVA1 <- lm(Y2~Y1+X1, data = d)
e.ANCOVA2 <- lm(dY~Y1+X1, data = d)

test_that("delta method for association based on residual variance", {

    dL2 <- reshape2::melt(d, id.vars = c("id","Y1","X1"),  measure.vars = c("Y1","Y2"))

    ## ANCOVA1
    e.lmmANCOVA1 <- lmm(value ~ variable + variable:X1, data = dL2, repetition = ~variable|id)
    e.coefANCOVA1 <- coef(e.lmmANCOVA1, effects  = "all")
    e.OmegaANCOVA1 <- sigma(e.lmmANCOVA1)
    e.vcovANCOVA1 <- vcov(e.lmmANCOVA1, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")

    ## sigma(e.ANCOVA1)^2 ## 2.201905
    ## e.OmegaANCOVA1[2,2]*(1-e.coefANCOVA1["rho(Y1,Y2)"]^2) ## 2.179437 
    ## prod(e.coefANCOVA1[c("sigma","k.Y2")])^2*(1-e.coefANCOVA1["rho(Y1,Y2)"]^2) ## 2.179437 
    
    ## ## check estimate
    ## (with Y1)
    expect_equal(unname(e.OmegaANCOVA1["Y1","Y2"]/e.OmegaANCOVA1["Y1","Y1"]), unname(coef(e.ANCOVA1)["Y1"]), tol = 1e-5)
    expect_equal(unname(e.coefANCOVA1["k.Y2"]*e.coefANCOVA1["rho(Y1,Y2)"]), unname(coef(e.ANCOVA1)["Y1"]), tol = 1e-5)

    ## (with X1)
    expect_equal(unname(e.coefANCOVA1["variableY2:X1"]-e.coefANCOVA1["k.Y2"]*e.coefANCOVA1["rho(Y1,Y2)"]*e.coefANCOVA1["variableY1:X1"]),
                 unname(coef(e.ANCOVA1)["X1"]), tol = 1e-5)

    ## ## check variance
    ## (with Y1)
    e.varANCOVA1 <- e.coefANCOVA1["rho(Y1,Y2)"]^2 * e.vcovANCOVA1["k.Y2","k.Y2"] + e.coefANCOVA1["k.Y2"]^2 * e.vcovANCOVA1["rho(Y1,Y2)","rho(Y1,Y2)"] + 2*e.coefANCOVA1["k.Y2"]*e.coefANCOVA1["rho(Y1,Y2)"]*e.vcovANCOVA1["k.Y2","rho(Y1,Y2)"]
    ## expect_equal(unname(sqrt(e.varANCOVA1)), unname(sqrt(diag(vcov(e.ANCOVA1))["Y1"])), tol = 1e-1)
    expect_equal(unname(sqrt(e.varANCOVA1)), 0.04586995, tol = 1e-5) ## too small variance compared to lm (0.04610579)
    
    e.deltaANCOVA1 <- estimate(e.lmmANCOVA1, function(p){
        c(Y1 = p["rho(Y1,Y2)"]*p["k.Y2"],
          X1 = p["variableY2:X1"]-p["k.Y2"]*p["rho(Y1,Y2)"]*p["variableY1:X1"])
    })
    e.deltaANCOVA1.bis <- estimate(e.lmmANCOVA1, function(p){
        Omega <- sigma(e.lmmANCOVA1, p = p)
        c(Y1 = Omega["Y1","Y2"]/Omega["Y1","Y1"],
          X1 = p["variableY2:X1"]-(Omega["Y1","Y2"]/Omega["Y1","Y1"])*p["variableY1:X1"])
    })
    
    
    expect_equal(e.deltaANCOVA1["Y1","estimate"], unname(e.OmegaANCOVA1["Y1","Y2"]/e.OmegaANCOVA1["Y1","Y1"]), tol = 1e-10)
    expect_equal(e.deltaANCOVA1["Y1","se"], unname(sqrt(e.varANCOVA1)), tol = 1e-10)

    test <- data.frame("estimate" = c(0.92166234, -0.16568066), 
                       "se" = c(0.04586995, 0.31929291), 
                       "df" = c(35.19182165, 69.08439034), 
                       "lower" = c(0.82855953, -0.80263871), 
                       "upper" = c(1.01476516, 0.47127739), 
                       "p.value" = c(0, 0.60548972))

    expect_equal(as.double(unlist(e.deltaANCOVA1)), as.double(unlist(test)), tol = 1e-5)

    ## ANCOVA2
    dL22 <- dL2
    dL22$Y1[dL22$variable=="Y1"] <- 0
   
    e.lmmANCOVA2 <- lmm(value ~ variable+variable:X1+Y1, data = dL22, repetition = ~variable|id, type.information = "expected")
    ## summary(lm(value ~ variable+variable:X1+Y1, data = dL22))
    ## summary(gls(value ~ variable+variable:X1+Y1, data = dL22, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|variable)))
    e.coefANCOVA2 <- coef(e.lmmANCOVA2, effects  = "all")
    e.OmegaANCOVA2 <- sigma(e.lmmANCOVA2)
    e.vcovANCOVA2 <- vcov(e.lmmANCOVA2, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")

    ## check estimate
    expect_equal(unname(e.coefANCOVA2["Y1"]+e.coefANCOVA2["k.Y2"]*e.coefANCOVA2["rho(Y1,Y2)"]-1), unname(coef(e.ANCOVA2)["Y1"]), tol = 1e-5)
    expect_equal(unname(e.coefANCOVA2["variableY2:X1"]-e.coefANCOVA2["k.Y2"]*e.coefANCOVA2["rho(Y1,Y2)"]*e.coefANCOVA2["variableY1:X1"]), unname(coef(e.ANCOVA2)["X1"]), tol = 1e-5)

    e.deltaANCOVA2 <- estimate(e.lmmANCOVA2, function(p){
        c(Y1 = as.double(p["Y1"]+p["k.Y2"]*p["rho(Y1,Y2)"]-1),
          X1 = as.double(p["variableY2:X1"]-p["k.Y2"]*p["rho(Y1,Y2)"]*p["variableY1:X1"]))
    })
    ## do not match standard error
    ## summary(e.ANCOVA2)$coef[c("Y1","X1"),]
    ## e.deltaANCOVA2

    test <- data.frame("estimate" = c(-0.07833768, -0.16568061), 
                       "se" = c(0.06440727, 0.34110244), 
                       "df" = c(380.89431576, 127.63222867), 
                       "lower" = c(-0.204976, -0.84062863), 
                       "upper" = c(0.04830065, 0.50926742), 
                       "p.value" = c(0.22462908, 0.62799783))
    expect_equal(as.double(unlist(test)), as.double(unlist(e.deltaANCOVA2)), tol = 1e-3)
})

## * Association between the changes
e.lm <- lm(dY~dX, data = d)
summary(e.lm)$coef["dX",]
##  Estimate Std. Error    t value   Pr(>|t|) 
## 0.2681921  0.2790260  0.9611725  0.3388310 


test_that("delta method for association based on residual variance", {

    dL2 <- reshape2::melt(d, id.vars = c("id","X5"),  measure.vars = c("dX","dY"))

    ## bivariate mixed model estimating the association between the changes
    e.lmm2 <- lmm(value ~ variable, data = dL2, repetition = ~variable|id)

    e.coef2 <- coef(e.lmm2, effects  = "all")
    e.Omega2 <- sigma(e.lmm2)
    e.vcov22 <- vcov(e.lmm2, effects = "all", type.information = "observed", transform.sigma = "none", transform.k = "none", transform.rho = "none")

    ## check estimate
    expect_equal(unname(e.Omega2["dY","dX"]/e.Omega2["dX","dX"]), unname(coef(e.lm)["dX"]), tol = 1e-5)
    expect_equal(unname(e.coef2["rho(dX,dY)"]*e.coef2["k.dY"]), unname(coef(e.lm)["dX"]), tol = 1e-5)

    ## check variance
    e.var2 <- e.coef2["rho(dX,dY)"]^2 * e.vcov22["k.dY","k.dY"] + e.coef2["k.dY"]^2 * e.vcov22["rho(dX,dY)","rho(dX,dY)"] + 2*e.coef2["k.dY"]*e.coef2["rho(dX,dY)"]*e.vcov22["k.dY","rho(dX,dY)"]
    ## expect_equal(unname(sqrt(e.var2)), unname(sqrt(diag(vcov(e.lm))["dX"])), tol = 1e-1)
    expect_equal(unname(sqrt(e.var2)), 0.2776132, tol = 1e-5) ## too small variance compared to lm (0.279)
    
    e.delta2 <- estimate(e.lmm2, function(p){
        unname(p["rho(dX,dY)"]*p["k.dY"])
    })
    expect_equal(e.delta2[,"estimate"], unname(e.Omega2["dY","dX"]/e.Omega2["dX","dX"]), tol = 1e-10)
    expect_equal(e.delta2[,"se"], unname(sqrt(e.var2)), tol = 1e-10)

    test <- data.frame("estimate" = c(0.26819079), 
                       "se" = c(0.27761318), 
                       "df" = c(45.32849395), 
                       "lower" = c(-0.29083904), 
                       "upper" = c(0.82722062), 
                       "p.value" = c(0.33914069))
    expect_equal(as.double(unlist(e.delta2)), as.double(unlist(test)), tol = 1e-5)
    
    ## quadrivariate mixed model estimating the association between the changes
    dL4 <- reshape2::melt(d, id.vars = c("id","X5"),  measure.vars = c("X1","X2","Y1","Y2"))
    e.lmm4 <- lmm(value ~ variable, data = dL4, repetition = ~variable|id)

    Omega4 <- sigma(e.lmm4)
    C <- rbind(c(1,-1,0,0), c(0,0,1,-1))
    Omega4.diff <- C %*% Omega4 %*% t(C)
    
    expect_equal(as.double(Omega4.diff), as.double(e.Omega2), tol = 1e-5)

    e.delta4 <- estimate(e.lmm4, function(p){ ## p <- coef(e.lmm4, effects = "all")
        iOmega <- C %*% sigma(e.lmm4, p = p) %*% t(C)
        iOmega[1,2]/iOmega[1,1]
    })
    test <- data.frame("estimate" = c(0.26819312), 
                       "se" = c(0.27761261), 
                       "df" = c(16.66677284), 
                       "lower" = c(-0.31841161), 
                       "upper" = c(0.85479785), 
                       "p.value" = c(0.34781928))
    expect_equal(as.double(unlist(e.delta4)), as.double(unlist(test)), tol = 1e-3)

})

## * Random effects
set.seed(10)
dL <- sampleRem(1e2, n.times = 3, format = "long")
dL$gender <- dL$id %% 2

test_that("delta method for random effects", {
    eRI.lmm <- lmm(Y ~ X1+X2+X3 + (1|id), repetition =~visit|id, data = dL)
    eRI.ranef <- ranef(eRI.lmm, effects = "mean", se = TRUE)
    eRI.tauA <- ranef(eRI.lmm, effects = "variance", scale = "absolute", se = TRUE)
    eRI.tauR <- ranef(eRI.lmm, effects = "variance", scale = "relative", se = TRUE)

    GS.ranef <- data.frame("id" = c(1, 2, 3, 4, 5, 6), 
                           "estimate" = c(0.81086415, 1.33742122, -2.58986262, -0.41766805, 6.63078022, 0.15015072), 
                           "se" = c(0.40338857, 0.40365815, 0.64672426, 0.59263204, 0.65224109, 0.44703011), 
                           "df" = c(95.8304692, 95.50855669, 95.27569894, 95.99610407, 91.49541552, 96.01409382), 
                           "lower" = c(0.01012608, 0.53611351, -3.8737247, -1.59403396, 5.33527792, -0.73719536), 
                           "upper" = c(1.61160221, 2.13872893, -1.30600054, 0.75869785, 7.92628253, 1.03749681))

    expect_equivalent(head(eRI.ranef), GS.ranef, tol = 1e-5)


    GS.tauA <- data.frame("type" = c("total", "id", "residual"), 
                          "estimate" = c(11.643744,  9.413364,  2.230380), 
                          "se" = c(1.473533, 1.467895, 0.223038), 
                          "lower" = c(9.085982, 6.934429, 1.833406), 
                          "upper" = c(14.921532, 12.778474,  2.713307))
    expect_equivalent(eRI.tauA, GS.tauA, tol = 1e-5)

    GS.tauR <- data.frame("type" = c("total", "id", "residual"), 
                          "estimate" = c(1.0000000, 0.8084482, 0.1915518), 
                          "se" = c(0.00000000, 0.02934011, 0.02934011), 
                          "lower" = c(1.0000000, 0.7529402, 0.1418754), 
                          "upper" = c(1.0000000, 0.8680484, 0.2586219))
    expect_equivalent(eRI.tauR, GS.tauR, tol = 1e-5)
})

##----------------------------------------------------------------------
### test-auto-estimate.R ends here
