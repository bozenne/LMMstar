### test-auto-estimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 31 2022 (11:36) 
## Version: 
## Last-Updated: jul 18 2025 (11:32) 
##           By: Brice Ozenne
##     Update #: 89
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
d$baseline <- d$Y1

## ** ANCOVA
e.ANCOVA1 <- lm(Y2~Y1+X1, data = d)
e.ANCOVA2 <- lm(dY~Y1+X1, data = d)

test_that("delta method for association based on residual variance", {

    dL2 <- reshape(d, direction = "long", idvar = c("id","baseline"),
                   timevar = "time", times = c("Y1","Y2"), varying = c("Y1","Y2"), 
                   v.names = "Y")
    dL2$baseline[dL2$time=="Y1"] <- 0
    rownames(dL2) <- NULL

    ## ANCOVA1
    e.lmmANCOVA1 <- lmm(Y ~ time + time:X1, data = dL2, repetition = ~time|id)
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
    expect_equal(unname(e.coefANCOVA1["timeY2:X1"]-e.coefANCOVA1["k.Y2"]*e.coefANCOVA1["rho(Y1,Y2)"]*e.coefANCOVA1["timeY1:X1"]),
                 unname(coef(e.ANCOVA1)["X1"]), tol = 1e-5)

    ## ## check variance
    ## (with Y1)
    e.varANCOVA1 <- e.coefANCOVA1["rho(Y1,Y2)"]^2 * e.vcovANCOVA1["k.Y2","k.Y2"] + e.coefANCOVA1["k.Y2"]^2 * e.vcovANCOVA1["rho(Y1,Y2)","rho(Y1,Y2)"] + 2*e.coefANCOVA1["k.Y2"]*e.coefANCOVA1["rho(Y1,Y2)"]*e.vcovANCOVA1["k.Y2","rho(Y1,Y2)"]
    ## expect_equal(unname(sqrt(e.varANCOVA1)), unname(sqrt(diag(vcov(e.ANCOVA1))["Y1"])), tol = 1e-1)
    expect_equal(unname(sqrt(e.varANCOVA1)), 0.04586995, tol = 1e-5) ## too small variance compared to lm (0.04610579)
    
    e.deltaANCOVA1 <- estimate(e.lmmANCOVA1, function(p){
        c(Y1 = p["rho(Y1,Y2)"]*p["k.Y2"],
          X1 = p["timeY2:X1"]-p["k.Y2"]*p["rho(Y1,Y2)"]*p["timeY1:X1"])
    })
    e.deltaANCOVA1.bis <- estimate(e.lmmANCOVA1, function(p){
        Omega <- sigma(e.lmmANCOVA1, p = p)
        c(Y1 = Omega["Y1","Y2"]/Omega["Y1","Y1"],
          X1 = p["timeY2:X1"]-(Omega["Y1","Y2"]/Omega["Y1","Y1"])*p["timeY1:X1"])
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
    e.lmmANCOVA2 <- lmm(Y ~ time+time:X1+baseline, data = dL2, repetition = ~time|id, type.information = "expected")
    ## summary(lm(value ~ variable+variable:X1+Y1, data = dL22))
    ## summary(gls(value ~ variable+variable:X1+Y1, data = dL22, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|variable)))
    e.coefANCOVA2 <- coef(e.lmmANCOVA2, effects  = "all")
    e.OmegaANCOVA2 <- sigma(e.lmmANCOVA2)
    e.vcovANCOVA2 <- vcov(e.lmmANCOVA2, effects = "all", transform.sigma = "none", transform.k = "none", transform.rho = "none")

    ## check estimate
    expect_equal(unname(e.coefANCOVA2["baseline"]+e.coefANCOVA2["k.Y2"]*e.coefANCOVA2["rho(Y1,Y2)"]-1), unname(coef(e.ANCOVA2)["Y1"]), tol = 1e-5)
    expect_equal(unname(e.coefANCOVA2["timeY2:X1"]-e.coefANCOVA2["k.Y2"]*e.coefANCOVA2["rho(Y1,Y2)"]*e.coefANCOVA2["timeY1:X1"]), unname(coef(e.ANCOVA2)["X1"]), tol = 1e-5)

    e.deltaANCOVA2 <- estimate(e.lmmANCOVA2, function(p){
        c(Y1 = as.double(p["baseline"]+p["k.Y2"]*p["rho(Y1,Y2)"]-1),
          X1 = as.double(p["timeY2:X1"]-p["k.Y2"]*p["rho(Y1,Y2)"]*p["timeY1:X1"]))
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

    dL2 <- reshape(d, direction = "long", idvar = "id",
                   timevar = "time", times = c("dX","dY"), varying = c("dX","dY"), 
                   v.names = "Y")

    ## bivariate mixed model estimating the association between the changes
    e.lmm2 <- lmm(Y ~ time, data = dL2, repetition = ~time|id)

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
    dL4 <- reshape(d, direction = "long", idvar = "id",
                   timevar = "variable", times = c("X1","X2","Y1","Y2"), varying = c("X1","X2","Y1","Y2"), 
                   v.names = "value")
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

## * Mediation analysis
set.seed(10)
dL <- sampleRem(1e2, n.times = 1, format = "long")

test_that("estimate for rbind.anova",{

    m.lvm <- lvm(Y ~ X5 + X1, X5 ~ X1)
    e.lvm <- estimate(m.lvm, data = dL)
    GS <- effects(e.lvm, Y~X1)
    GS.robust <- estimate(e.lvm, f = function(p){
        DE <- as.double(p["Y~X1"])
        IE <- as.double(p["Y~X5"]*p["X5~X1"])
        out <- c(direct = DE,
                 indirect = IE,
                 total = DE + IE,
                 proportion = 1/(1 + DE/IE))
    }, robust = TRUE)

    e.lmm1 <- lmm(Y ~ X5 + X1, repetition = ~visit|id, data = dL, method = "ML")
    e.aov1 <- anova(e.lmm1, effects = c("X1=0","X5=0"), simplify = FALSE)
    e.aov1.robust <- anova(e.lmm1, effects = c("X1=0","X5=0"), simplify = FALSE, robust = TRUE)
    e.lmm2 <- lmm(X5 ~ X1, repetition = ~visit|id, data = dL, method = "ML")
    e.aov2 <- anova(e.lmm2, effects = "X1=0", simplify = FALSE)
    e.aov2.robust <- anova(e.lmm2, effects = "X1=0", simplify = FALSE, robust = TRUE)
    e.maov <- rbind(e.aov1, e.aov2, name = c("Y","X5"))
    e.maov.robust <- rbind(e.aov1.robust, e.aov2.robust, name = c("Y","X5"))

    expect_equivalent(diag(vcov(e.lvm)[c("Y~X1","Y~X5","X5~X1"),c("Y~X1","Y~X5","X5~X1")]),diag(vcov(e.maov)), tol = 1e-6)
    expect_equivalent(crossprod(iid(e.lvm))[c("Y~X1","Y~X5","X5~X1"),c("Y~X1","Y~X5","X5~X1")], vcov(e.maov.robust), tol = 1e-6)

    test <- estimate(e.maov.robust, f = function(object){
        p <- coef(object)
        DE <- as.double(p["Y: X1=0"])
        IE <- as.double(p["Y: X5=0"]*p["X5: X1=0"])
        out <- c(direct = DE,
                 indirect = IE,
                 total = DE + IE,
                 proportion = 1/(1 + DE/IE))
    })
    test.robust <- estimate(e.maov.robust, f = function(object){
        p <- coef(object)
        DE <- as.double(p["Y: X1=0"])
        IE <- as.double(p["Y: X5=0"]*p["X5: X1=0"])
        out <- c(direct = DE,
                 indirect = IE,
                 total = DE + IE,
                 proportion = 1/(1 + DE/IE))
    }, df = FALSE)
    expect_equal(as.double(GS$coef)[1:3], as.double(test.robust$estimate)[c(3,1,2)], tol = 1e-5)
    expect_equivalent(as.data.frame(GS.robust$coefmat), test.robust[,c("estimate","se","lower","upper","p.value")], tol = 1e-5)
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
