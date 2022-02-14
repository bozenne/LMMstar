### test-auto-estimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 31 2022 (11:36) 
## Version: 
## Last-Updated: feb 14 2022 (09:28) 
##           By: Brice Ozenne
##     Update #: 33
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

context("Check delta method (estimate function) for mixed model")

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

    dL2 <- melt(d, id.vars = c("id","Y1","X1"),  measure.vars = c("Y1","Y2"))

    ## ANCOVA1
    e.lmmANCOVA1 <- lmm(value ~ variable + variable:X1, data = dL2, repetition = ~variable|id)
    e.coefANCOVA1 <- coef(e.lmmANCOVA1, effects  = "all")
    e.OmegaANCOVA1 <- getVarCov(e.lmmANCOVA1)
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
    expect_equal(unname(sqrt(e.varANCOVA1)), 0.04586995, tol = 1e-6) ## too small variance compared to lm (0.04610579)
    
    e.deltaANCOVA1 <- estimate(e.lmmANCOVA1, function(p){
        c(Y1 = p["rho(Y1,Y2)"]*p["k.Y2"],
          X1 = p["variableY2:X1"]-p["k.Y2"]*p["rho(Y1,Y2)"]*p["variableY1:X1"])
    })
    
    expect_equal(e.deltaANCOVA1["Y1","estimate"], unname(e.OmegaANCOVA1["Y1","Y2"]/e.OmegaANCOVA1["Y1","Y1"]), tol = 1e-10)
    expect_equal(e.deltaANCOVA1["Y1","se"], unname(sqrt(e.varANCOVA1)), tol = 1e-10)

    test <- data.frame("estimate" = c(0.92166234, -0.16568066), 
                       "se" = c(0.04586995, 0.31929291), 
                       "df" = c(35.19182165, 69.08439034), 
                       "lower" = c(0.82855953, -0.80263871), 
                       "upper" = c(1.01476516, 0.47127739), 
                       "p.value" = c(0, 0.60548972))

    expect_equal(as.double(unlist(e.deltaANCOVA1)), as.double(unlist(test)), tol = 1e-6)

    ## ANCOVA2
    dL22 <- dL2
    dL22$Y1[dL22$variable=="Y1"] <- 0
   
    e.lmmANCOVA2 <- lmm(value ~ variable+variable:X1+Y1, data = dL22, repetition = ~variable|id, type.information = "expected")
    ## summary(lm(value ~ variable+variable:X1+Y1, data = dL22))
    ## summary(gls(value ~ variable+variable:X1+Y1, data = dL22, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|variable)))
    e.coefANCOVA2 <- coef(e.lmmANCOVA2, effects  = "all")
    e.OmegaANCOVA2 <- getVarCov(e.lmmANCOVA2)
    e.vcovANCOVA2 <- vcov(e.lmmANCOVA2, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")

    ## check estimate
    expect_equal(unname(e.coefANCOVA2["Y1"]+e.coefANCOVA2["k.Y2"]*e.coefANCOVA2["rho(Y1,Y2)"]-1), unname(coef(e.ANCOVA2)["Y1"]), tol = 1e-5)
    expect_equal(unname(e.coefANCOVA2["variableY2:X1"]-e.coefANCOVA2["k.Y2"]*e.coefANCOVA2["rho(Y1,Y2)"]*e.coefANCOVA2["variableY1:X1"]), unname(coef(e.ANCOVA2)["X1"]), tol = 1e-5)

    e.deltaANCOVA2 <- estimate(e.lmmANCOVA2, function(p){
        c(Y1 = as.double(p["Y1"]+p["k.Y2"]*p["rho(Y1,Y2)"]-1),
          X1 = as.double(p["variableY2:X1"]-p["k.Y2"]*p["rho(Y1,Y2)"]*p["variableY1:X1"]))
    })
    ## do not match standard error
    summary(e.ANCOVA2)$coef[c("Y1","X1"),]
    ## e.deltaANCOVA2

    test <- data.frame("estimate" = c(-0.07833768, -0.16568061), 
                       "se" = c(0.06440726, 0.34110244), 
                       "df" = c(386.69874806, 129.57726267), 
                       "lower" = c(-0.20496993, -0.84053168), 
                       "upper" = c(0.04829458, 0.50917047), 
                       "p.value" = c(0.22461774, 0.62798535))
    expect_equal(as.double(unlist(test)), as.double(unlist(e.deltaANCOVA2)), tol = 1e-6)
})


## * Association between the changes
e.lm <- lm(dY~dX, data = d)
summary(e.lm)$coef["dX",]
##  Estimate Std. Error    t value   Pr(>|t|) 
## 0.2681921  0.2790260  0.9611725  0.3388310 


test_that("delta method for association based on residual variance", {

    dL2 <- melt(d, id.vars = c("id","X5"),  measure.vars = c("dX","dY"))


    ## bivariate mixed model estimating the association between the changes
    e.lmm2 <- lmm(value ~ variable, data = dL2, repetition = ~variable|id)

    e.coef2 <- coef(e.lmm2, effects  = "all")
    e.Omega2 <- getVarCov(e.lmm2)
    e.vcov22 <- vcov(e.lmm2, effects = "all", type.information = "observed", transform.sigma = "none", transform.k = "none", transform.rho = "none")

    ## check estimate
    expect_equal(unname(e.Omega2["dY","dX"]/e.Omega2["dX","dX"]), unname(coef(e.lm)["dX"]), tol = 1e-5)
    expect_equal(unname(e.coef2["rho(dX,dY)"]*e.coef2["k.dY"]), unname(coef(e.lm)["dX"]), tol = 1e-5)

    ## check variance
    e.var2 <- e.coef2["rho(dX,dY)"]^2 * e.vcov22["k.dY","k.dY"] + e.coef2["k.dY"]^2 * e.vcov22["rho(dX,dY)","rho(dX,dY)"] + 2*e.coef2["k.dY"]*e.coef2["rho(dX,dY)"]*e.vcov22["k.dY","rho(dX,dY)"]
    ## expect_equal(unname(sqrt(e.var2)), unname(sqrt(diag(vcov(e.lm))["dX"])), tol = 1e-1)
    expect_equal(unname(sqrt(e.var2)), 0.2776132, tol = 1e-6) ## too small variance compared to lm (0.279)
    
    e.delta2 <- estimate(e.lmm2, function(p){
        unname(p["rho(dX,dY)"]*p["k.dY"])
    })
    expect_equal(e.delta2[,"estimate"], unname(e.Omega2["dY","dX"]/e.Omega2["dX","dX"]), tol = 1e-10)
    expect_equal(e.delta2[,"se"], unname(sqrt(e.var2)), tol = 1e-10)

    test <- data.frame("estimate" = c(0.26819079), 
                       "se" = c(0.27761318), 
                       "df" = c(45.3428542), 
                       "lower" = c(-0.29083419), 
                       "upper" = c(0.82721577), 
                       "p.value" = c(0.33913908))
    expect_equal(as.double(unlist(e.delta2)), as.double(unlist(test)), tol = 1e-6)
    
    ## quadrivariate mixed model estimating the association between the changes
    dL4 <- melt(d, id.vars = c("id","X5"),  measure.vars = c("X1","X2","Y1","Y2"))
    e.lmm4 <- lmm(value ~ variable, data = dL4, repetition = ~variable|id)

    Omega4 <- getVarCov(e.lmm4)
    C <- rbind(c(1,-1,0,0), c(0,0,1,-1))
    Omega4.diff <- C %*% Omega4 %*% t(C)
    
    expect_equal(as.double(Omega4.diff), as.double(e.Omega2), tol = 1e-5)

    e.delta4 <- estimate(e.lmm4, function(p){ ## p <- coef(e.lmm4, effects = "all")
        iOmega <- C %*% getVarCov(e.lmm4, p = p) %*% t(C)
        iOmega[1,2]/iOmega[1,1]
    })
    test <- data.frame("estimate" = c(0.26819312), 
                       "se" = c(0.27761253), 
                       "df" = c(16.65819987), 
                       "lower" = c(-0.31843493), 
                       "upper" = c(0.85482116), 
                       "p.value" = c(0.34782614))
    expect_equal(as.double(unlist(e.delta4)), as.double(unlist(test)), tol = 1e-6)

})

## ## For Paul
## set.seed(10)
## dW <- sampleRem(25, n.time = 3)
## dW$dY2 <- dW$Y2 - dW$Y1
## dW$dY3 <- dW$Y3 - dW$Y1

## dL <- melt(dW, id.vars = c("id","X1","Y1"), measure.vars = c("dY2","dY3"))

## e.lm <- lm(dY3~Y1+X1, data = dW)


## e.lmm <- lmm(value ~ variable + variable:Y1 + variable:X1, repetition =~variable|id,  data = dL)
## library(nlme)
## e.lme <- lme(value ~ variable + variable:Y1 + variable:X1, random =~1|id, weights=varIdent(form=~1|variable),
##              correlation=corSymm(form=~1|id),  data = dL)

## summary(e.lm)$coef
## print(model.tables(e.lmm), digit = 7)

## summary(e.lme)$tTable


##----------------------------------------------------------------------
### test-auto-estimate.R ends here
