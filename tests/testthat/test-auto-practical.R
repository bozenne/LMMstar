### test-auto-practical.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  7 2021 (17:03) 
## Version: 
## Last-Updated: maj  7 2024 (15:18) 
##           By: Brice Ozenne
##     Update #: 131
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
    library(numDeriv)
    library(lava)
    library(multcomp)
    library(nlme)
    library(qqtest)
    library(ggplot2)
    library(lme4)
    library(lmerTest)

    library(LMMstar)
}

context("Check lmm on the example from the course with Julie")
LMMstar.options(optimizer = "FS", method.numDeriv = "simple", precompute.moments = TRUE,
                columns.confint = c("estimate","se","df","lower","upper","p.value"))

## * Practical 1
test_that("practical 1 - gastricbypass",{

    ## ** data
    data(gastricbypassL, package = "LMMstar")

    ## ** summarize
    summarize(glucagonAUC~time, data = gastricbypassL, na.rm=TRUE)
    summarize(glucagonAUC~time|id, data = gastricbypassL, na.rm=TRUE)
    
    ## ** compound symmetry
    eCS.gls <- gls(glucagonAUC~visit,
                   data=gastricbypassL,
                   correlation=corCompSymm(form=~visit|id),
                   na.action=na.exclude,
                   control=glsControl(opt="optim"))

    eCS.lmm <- lmm(glucagonAUC~visit,
                   data=gastricbypassL,
                   repetition = ~visit|id,
                   structure = "CS")
    
    expect_equal(as.double(logLik(eCS.gls)), as.double(logLik(eCS.lmm)), tol = 1e-6)

    ## ** unstructured with missing data
    eUN.gls <- gls(glucagonAUC~visit,
                   data=gastricbypassL,
                   correlation=corSymm(form=~as.numeric(visit)|id),
                   weights=varIdent(form=~1|time),
                   na.action=na.exclude,
                   control=glsControl(opt="optim"),
                   method = "REML")

    eUN.lmm <- lmm(glucagonAUC~visit,
                   data=gastricbypassL,
                   repetition = ~visit|id,
                   structure = "UN",
                   method.fit = "REML", trace = 0)

    ## check moments
    expect_equal(as.double(logLik(eUN.gls)), as.double(logLik(eUN.lmm)), tol = 1e-6)
    
    ## GS <- numDeriv::jacobian(func = function(p){logLik(eUN.lmm, p = p, transform.sigma = "none", transform.k = "none", transform.rho = "none")}, x = coef(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all"))
    GS <- rbind(c("(Intercept)" = 0, "visit2" = 0, "visit3" = 0, "visit4" = 0, "sigma" = 2e-08, "k.2" = 7.85e-06, "k.3" = 5.37e-06, "k.4" = -2.109e-05, "rho(1,2)" = -0.00011073, "rho(1,3)" = -7.107e-05, "rho(1,4)" = 9.35e-06, "rho(2,3)" = 7.852e-05, "rho(2,4)" = -2.407e-05, "rho(3,4)" = 4.945e-05))
    expect_equal(as.double(GS), as.double(score(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all")), tol = 1e-6)

    ## GS <- -numDeriv::jacobian(func = function(p){score(eUN.lmm, p = p, transform.sigma = "none", transform.k = "none", transform.rho = "none",effects = "all")}, x = coef(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all"))
    GS <- matrix(c(0.2000435, 0.10068008, -0.00122474, 0.0667917, 0, 0.01960874, -0.00398498, 0.00220684, 0.015223, 0.02694338, -0.00181844, -0.02394494, -0.00163313, 0.00522044, 0.10068008, 0.39828548, 0.05620903, -0.04323161, 0, 0.06528472, -0.03616399, 0.07545552, 0.34097272, 0.31291727, -0.03538287, -0.31640674, 0.0569075, -0.1180618, -0.00122474, 0.05620903, 0.05070213, -0.04276096, 0, 0.00060657, -0.01939405, 0.0573912, 0.24479209, 0.18869963, -0.02474055, -0.19994584, 0.05242174, -0.1138332, 0.0667917, -0.04323161, -0.04276096, 0.08432189, 0, 0.00029111, 0.01617431, -0.04825539, -0.20558724, -0.15785621, 0.02076676, 0.16745253, -0.04422625, 0.09610526, 0, 0, 0, 0, 0.61232751, 2.68368828, 1.30246228, 1.80652473, -7.88842064, -3.41692146, 2.72607489, 2.53245722, -1.48188657, -3.02859195, 0.01960874, 0.06528472, 0.00060657, 0.00029111, 2.68368828, 119.6435905, -1.30957174, 2.67549682, -67.36667611, 2.7100885, -1.4963409, 22.4424789, -13.18941762, -0.05361995, -0.00398498, -0.03616399, -0.01939405, 0.01617431, 1.30246228, -1.30957174, 17.62474246, -5.70108564, 2.50057785, -13.49661757, -1.26046703, 9.8287541, 1.16522916, -13.33546581, 0.00220684, 0.07545552, 0.0573912, -0.04825539, 1.80652473, 2.67549682, -5.70108564, 29.08095721, -4.89725928, -3.62086875, 17.33682661, 3.88337218, -10.663066, -16.12586692, 0.015223, 0.34097272, 0.24479209, -0.20558724, -7.88842064, -67.36667611, 2.50057785, -4.89725928, 494.48230881, 170.86268327, -126.83334603, -160.48986839, 120.86002165, 52.61300133, 0.02694338, 0.31291727, 0.18869963, -0.15785621, -3.41692146, 2.7100885, -13.49661757, -3.62086875, 170.86268327, 195.7773562, -137.85293967, -150.42933562, 105.91037437, 88.19350765, -0.00181844, -0.03538287, -0.02474055, 0.02076676, 2.72607489, -1.4963409, -1.26046703, 17.33682661, -126.83334603, -137.85293967, 182.19196649, 101.08267288, -138.37463781, -84.00686145, -0.02394494, -0.31640674, -0.19994584, 0.16745253, 2.53245722, 22.4424789, 9.8287541, 3.88337218, -160.48986839, -150.42933562, 101.08267288, 155.69993421, -103.48124038, -60.02043103, -0.00163313, 0.0569075, 0.05242174, -0.04422625, -1.48188657, -13.18941762, 1.16522916, -10.663066, 120.86002165, 105.91037437, -138.37463781, -103.48124038, 144.14781848, 59.2350391, 0.00522044, -0.1180618, -0.1138332, 0.09610526, -3.02859195, -0.05361995, -13.33546581, -16.12586692, 52.61300133, 88.19350765, -84.00686145, -60.02043103, 59.2350391, 99.27770758), 
                 nrow = 14, 
                 ncol = 14, 
                 dimnames = list(c("(Intercept)", "visit2", "visit3", "visit4", "sigma", "k.2", "k.3", "k.4", "rho(1,2)", "rho(1,3)", "rho(1,4)", "rho(2,3)", "rho(2,4)", "rho(3,4)"),c("(Intercept)", "visit2", "visit3", "visit4", "sigma", "k.2", "k.3", "k.4", "rho(1,2)", "rho(1,3)", "rho(1,4)", "rho(2,3)", "rho(2,4)", "rho(3,4)")) 
                 )
    expect_equal(as.double(GS), as.double(information(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all")), tol = 1e-6)

    ## ** extract information
    capture.output(model.tables(eUN.lmm))
    capture.output(summary(eUN.lmm, columns = "statistic"))
    capture.output(confint(eUN.lmm, effects = "all", backtransform = TRUE))
    expect_equal(unname(coef(eUN.lmm, effects = "variance", transform.k = "sd")),
                 unname(coef(eUN.lmm, effects = "variance")[1]*c(1,coef(eUN.lmm, effects = "variance")[-1])),
                 tol = 1e-6)
    autoplot(eUN.lmm)

    apply(residuals(eUN.lmm, type = "normalized", format = "wide")[,-1],2,sd,na.rm=TRUE)
    cor(residuals(eUN.lmm, type = "normalized", format = "wide")[,-1], use = "pairwise")
    qqtest::qqtest(na.omit(residuals(eUN.lmm, type = "normalized")))

    eUN.lmm_anova <- anova(eUN.lmm, effects = c("visit3-visit2=0"), ci = TRUE)
    capture.output(eUN.lmm_anova)
    capture.output(summary(eUN.lmm_anova))
    expect_equal(eUN.lmm_anova$multivariate$df.denom, 19.20982, tol = 1e-1) ## Richardson
})

## * Practical 2
test_that("practical 2 - ncgs",{
    
    data(ncgsL, package = "LMMstar")
    ncgsL$visit <- as.numeric(ncgsL$visit)
    ncgsL$highdose.time <- ncgsL$time
    ncgsL$highdose.time[ncgsL$group=="placebo"] <- "0"
    ncgsL$time <- relevel(as.factor(ncgsL$time), ref="0")
    ncgsL$highdose.time <- relevel(as.factor(ncgsL$highdose.time), ref="0")
    
    ncgsL$treatment <- factor(ncgsL$group, c("none","placebo","highdose"))
    ncgsL$treatment[ncgsL$time=="0"] <- "none"
    
    ncgsL$treatment2 <- factor(ncgsL$group, c("placebo","highdose"))
    ncgsL$treatment2[ncgsL$time=="0"] <- "placebo"

    ## ** unstructured with missing data
    e.gls <- gls(cholest~time+highdose.time,
                 data=ncgsL,
                 correlation=corSymm(form=~visit|id),
                 weights=varIdent(form=~1|time),
                 na.action=na.exclude,
                 control=glsControl(opt='optim'))

    e.lmm <- lmm(cholest~time+highdose.time,
                 data=ncgsL,
                 repetition = ~time|id,
                 structure = "UN",
                 df = FALSE)
    
    e2.lmm <- suppressMessages(lmm(cholest~treatment*time,
                                   data=ncgsL,
                                   repetition = ~time|id,
                                   structure = "UN"))
    
    e3.lmm <- suppressMessages(lmm(cholest~treatment2*time,
                                   data=ncgsL,
                                   repetition = ~time|id,
                                   structure = "UN"))

    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lmm)), tol = 1e-6)
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e2.lmm)), tol = 1e-6)
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e3.lmm)), tol = 1e-6)
    
    ## ** extract information
    confint(e2.lmm, effects = "all", backtransform = TRUE)[,c("estimate","lower","upper")]
    e2.lmm_anova <- anova(e2.lmm, effects = c("treatmenthighdose:time6-treatmentplacebo:time6=0","treatmenthighdose:time12-treatmentplacebo:time12=0"), ci = TRUE)
    expect_equal(e2.lmm_anova$multivariate$df.denom, 100.0411, tol = 1e-1) ## Richardson
    autoplot(e2.lmm) 
    autoplot(e2.lmm, color = "group") 
    autoplot(e2.lmm, color = "group", ci.alpha = 0.25)
})

test_that("practical 2 - vitamin",{
    data(vitaminL, package = "LMMstar")

    vitaminL$visit.num <- as.numeric(vitaminL$visit)
    
    ## define treatment
    vitaminL$treatment <- factor(vitaminL$group, c("none","control","vitamin"))
    vitaminL$treatment[vitaminL$time<=4] <- "none"

    ## define interaction manually
    vitaminL$vita.time <- vitaminL$time
    vitaminL$vita.time[vitaminL$group=="control" | vitaminL$time<=4] <- "1"

    ## ** fit unstructured
    e.gls <- gls(weight~visit+vita.time,
                 data=vitaminL,
                 correlation=corSymm(form=~visit.num|animal),
                 weights=varIdent(form=~1|visit),
                 na.action=na.exclude,
                 control=glsControl(opt='optim'))

    e0.lmm <- suppressWarnings(lmm(weight~visit+vita.time,
                                   data=vitaminL,
                                   repetition = ~visit|animal,
                                   structure = "UN"))
    
    e.lmm <- suppressMessages(suppressWarnings(lmm(weight~treatment*visit,
                                                   data=vitaminL,
                                                   repetition = ~visit|animal,
                                                   structure = "UN")))

    expect_equal(as.double(logLik(e0.lmm)), as.double(logLik(e.lmm)), tol = 1e-6)
    expect_equal(as.double(logLik(e.lmm)), -234.28331772, tol = 1e-6)
    expect_equal(as.double(logLik(e.gls)), -232.08183621, tol = 1e-6)
    
    ## ** extract information
    autoplot(e.lmm, color = "group",ci =FALSE)
    qqtest::qqtest(residuals(e.lmm))
    confint(e.lmm)
    e.lmm_anova <- anova(e.lmm, effects = "treatmentvitamin:visit6 - treatmentcontrol:visit6 = 0", ci = TRUE)
    ## expect_equal(e.lmm_anova$all$df.denom, 0.4972829, tol = 1e-1) ## Richardson
   
})

## * Practical 3
test_that("practical 3 - swabsL",{
   
    data(swabsL, package = "LMMstar")

    ## ** unstructured
    eUN.gls <- gls(swabs ~ crowding + name,
                   data = swabsL,
                   correlation = corSymm(form = ~ as.numeric(name) | family),
                   weights = varIdent(form = ~ 1 | name)
                   )
    eUN.lmm <- lmm(swabs ~ crowding + name, data = swabsL, structure = "UN", repetition = ~name|family)
    expect_equal(as.double(logLik(eUN.lmm)),as.double(logLik(eUN.gls)), tol = 1e-6)
    ## summary(eUN.lmm)
    sigma(eUN.lmm)

    ## ** compound symmetry
    eCS.lme <- lme(swabs ~ crowding + name,
                   data = swabsL,
                   random =~ 1 | family)
    eCS.lmm <- lmm(swabs ~ crowding + name, data = swabsL, structure = "CS", repetition = ~name|family)
    expect_equal(as.double(logLik(eCS.lmm)), as.double(logLik(eCS.lme)), tol = 1e-6)

    capture.output(summary(eCS.lmm))
    sigma(eCS.lmm)
    autoplot(eCS.lmm)

    ## with interaction
    eCSI.lmm <- lmm(swabs ~ crowding * name, data = swabsL, structure = "CS", repetition = ~name|family)

    eRI.lmm <- lmm(swabs ~ crowding * name + (1|family), data = swabsL)
    eCSI.lme <- lme(swabs ~ crowding * name, random =~ 1 | family, data = swabsL)
    GS <- as.data.frame(ranef(eCSI.lme, augFrame = TRUE))

    expect_equal(as.double(GS[,1]),as.double(ranef(eRI.lmm)), tol = 1e-4)
})

## * Practical 4
test_that("practical 4 - bloodpressureL",{

    data(bloodpressureL, package = "LMMstar")
    bloodpressureL$period.num <- as.numeric(bloodpressureL$period)
    bloodpressureL$treatment.num <- as.numeric(bloodpressureL$treatment)
            
    ## ** compound symmetry
    eCS.gls <- gls(duration ~ period + treatment,
                   data = bloodpressureL,
                   correlation = corCompSymm(form=~ 1 | id))
    eCS.lmm <- lmm(duration ~ period + treatment,
                   data = bloodpressureL,
                   structure = "CS", repetition = ~ period | id)
    expect_equal(as.double(logLik(eCS.lmm)),as.double(logLik(eCS.gls)), tol = 1e-6)
    plot(eCS.lmm, color = "sequence", obs.alpha = 0.25, ci.alpha = 0.05)
    confint(eCS.gls)
    confint(eCS.lmm)

    ## ** unstructured
    eUNP.gls <- gls(duration ~ period + treatment,
                    data = bloodpressureL,
                    correlation = corSymm(form=~ period.num | id),
                    weights = varIdent(form=~ 1|period),
                    )
    eUNP.lmm <- lmm(duration ~ period + treatment,
                    data = bloodpressureL,
                    structure = "UN", repetition = ~ period | id)
    expect_equal(as.double(logLik(eUNP.lmm)),as.double(logLik(eUNP.gls)), tol = 1e-6)
    capture.output(summary(eUNP.lmm))

    eUNT.gls <- gls(duration ~ period + treatment,
                    data = bloodpressureL,
                    correlation = corSymm(form=~ treatment.num | id),
                    weights = varIdent(form=~ 1|treatment),
                    )
    eUNT.lmm <- lmm(duration ~ period + treatment,
                    data = bloodpressureL,
                    structure = "UN", repetition = ~ treatment | id)
    capture.output(summary(eUNT.lmm))
    expect_equal(as.double(logLik(eUNT.lmm)),as.double(logLik(eUNT.gls)), tol = 1e-6)

})

## * Practical 6
test_that("practical 6 - vasscoresL",{
    data(vasscoresL, package = "LMMstar")

    summarize(vas ~ treatment, data = vasscoresL, na.rm = TRUE)
    summarize(vas ~ treatment|id, data = vasscoresL, na.rm = TRUE)
    summarize(vas ~ treatment+group|id, data = vasscoresL, na.rm = TRUE)

    ## ** model on all pairs: AB, AC, BC
    fit.CS <- lmm(vas~-1+treatment, data=vasscoresL,
                  repetition=~treatment|id, structure="CS")
    capture.output(summary(fit.CS))
    
    GS <- matrix(c(2099.98852585, 1675.79243303, 1675.79243303, 1675.79243303, 2099.98852585, 1675.79243303, 1675.79243303, 1675.79243303, 2099.98852585), 
                 nrow = 3, 
                 ncol = 3, 
                 dimnames = list(c("A", "B", "C"),c("A", "B", "C")) 
                 ) 

    expect_equivalent(sigma(fit.CS), GS, tol = 1e-5)
    
    ## GS <- lmer(vas~-1+treatment+(1|id), data=vasscoresL)    
    GS <- lme(vas~-1+treatment, random =~ 1|id, data=vasscoresL, na.action = na.omit)    
    expect_equivalent(unname(coef(fit.CS)), unname(fixef(GS)), tol = 1e-5)
    ## expect_equivalent(unname(model.tables(fit.CS)$se), unname(summary(GS)$coef[,"Std. Error"]), tol = 1e-2)
    ## expect_equivalent(unname(model.tables(fit.CS)$df), unname(summary(GS)$coef[,"df"]), tol = 1e-2)
    expect_equivalent(unname(model.tables(fit.CS)$se), unname(summary(GS)$tTable[,"Std.Error"]), tol = 1e-2)
    expect_equivalent(unname(model.tables(fit.CS)$df), c(43.64236203, 43.64236203, 43.64236203), tol = 1e-2)

    ## autoplot(fit.CS)
    suppressWarnings(plot(fit.CS, obs.alpha = 0.1))

    fit.UN <- lmm(vas~-1+treatment, data=vasscoresL,
                  repetition=~treatment|id, structure="UN")
    capture.output(summary(fit.UN))

    ## GS <- lmm(vas~-1+treatment, data=vasscoresL,
    ##               repetition=~treatment|id, structure="UN",
    ##               control = list(optimizer = "FS"))

    GS <- matrix(c(2038.74222977, 1892.18702336, 1797.44635389, 1892.18702336, 2111.43949953, 1369.79020341, 1797.44635389, 1369.79020341, 2157.87873809), 
                 nrow = 3, 
                 ncol = 3, 
                 dimnames = list(c("A", "B", "C"),c("A", "B", "C")) 
                 )
    expect_equivalent(sigma(fit.UN), GS, tol = 1e-5)

    vasscoreL.imputed <- fitted(fit.UN, type = "outcome", keep.data = TRUE)
    set.seed(11)
    gg <- ggplot(vasscoreL.imputed, aes(x = treatment, y = vas, group = id))
    gg + geom_point(aes(color = impute)) + geom_line()

    ## ** model some pairs: AB, AC
    ## IMPORTANT: check extraction of the residual variance-covariance matrix when not all pairs of time are analyzed
    vasscoreL.red <- vasscoresL[vasscoresL$group %in% c("AB","AC"),]

    fit.CS.red <- lmm(vas~-1+treatment, data=vasscoreL.red,
                      repetition=~treatment|id, structure="CS")
    capture.output(summary(fit.CS.red))
    
    GS <- matrix(c(2006.9354323, 1763.57527658, 1763.57527658, 1763.57527658, 2006.9354323, 1763.57527658, 1763.57527658, 1763.57527658, 2006.9354323), 
                 nrow = 3, 
                 ncol = 3, 
                 dimnames = list(c("A", "B", "C"),c("A", "B", "C")) 
                 ) 

    expect_equivalent(sigma(fit.CS.red), GS, tol = 1e-5)

    fit.UN.red <- lmm(vas~-1+treatment, data=vasscoreL.red,
                      repetition=~treatment|id, structure="UN", control = list(optimizer = "FS")) ## FS optimizer mandatory, GLS does something weird as it is able to estimate the correlation (B,C)
    capture.output(summary(fit.UN.red))
    
    GS <- matrix(c(1996.68088251, 1926.29216631, 1573.59329163, 1926.29216631, 2231.4711806, NA, 1573.59329163, NA, 1753.88563313), 
                 nrow = 3, 
                 ncol = 3, 
                 dimnames = list(c("A", "B", "C"),c("A", "B", "C")) 
                 ) 

    expect_equivalent(sigma(fit.UN.red), GS, tol = 1e-5)
})
##----------------------------------------------------------------------
### test-auto-practical.R ends here
