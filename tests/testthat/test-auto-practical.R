### test-auto-practical.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  7 2021 (17:03) 
## Version: 
## Last-Updated: jul 18 2025 (16:32) 
##           By: Brice Ozenne
##     Update #: 140
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
    GS <- rbind(c("(Intercept)" = 0, "visit2" = 0, "visit3" = 0, "visit4" = 0, "sigma" = 0, "k.2" = 1.24e-06, "k.3" = 8.5e-07, "k.4" = -3.32e-06, "rho(1,2)" = -1.746e-05, "rho(1,3)" = -1.12e-05, "rho(1,4)" = 1.47e-06, "rho(2,3)" = 1.238e-05, "rho(2,4)" = -3.8e-06, "rho(3,4)" = 7.8e-06))
    expect_equal(as.double(GS), as.double(score(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all")), tol = 1e-6)

    ## GS <- -numDeriv::jacobian(func = function(p){score(eUN.lmm, p = p, transform.sigma = "none", transform.k = "none", transform.rho = "none",effects = "all")}, x = coef(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all"))
    GS <- matrix(c(0.20004352, 0.10067992, -0.00122488, 0.06679182, 0, 0.0196087, -0.0039849, 0.00220668, 0.0152223, 0.02694279, -0.00181838, -0.02394434, -0.00163328, 0.00522075, 0.10067992, 0.39828487, 0.0562088, -0.04323152, 0, 0.06528458, -0.03616402, 0.07545536, 0.3409704, 0.31291605, -0.0353833, -0.31640557, 0.05690788, -0.11806133, -0.00122488, 0.0562088, 0.05070214, -0.04276106, 0, 0.00060655, -0.01939419, 0.05739131, 0.24479116, 0.18869954, -0.02474101, -0.19994578, 0.05242226, -0.11383326, 0.06679182, -0.04323152, -0.04276106, 0.08432205, 0, 0.00029112, 0.01617447, -0.04825559, -0.2055869, -0.15785648, 0.02076719, 0.16745284, -0.04422678, 0.09610551, 0, 0, 0, 0, 0.61232758, 2.68368752, 1.30246127, 1.80652693, -7.88839364, -3.41690543, 2.72607278, 2.53243921, -1.48188089, -3.0286055, 0.0196087, 0.06528458, 0.00060655, 0.00029112, 2.68368752, 119.64341672, -1.30956579, 2.67549075, -67.36647678, 2.71012078, -1.49637096, 22.44235147, -13.1893562, -0.05360568, -0.0039849, -0.03616402, -0.01939419, 0.01617447, 1.30246127, -1.30956579, 17.62472425, -5.70109557, 2.50063266, -13.4965165, -1.26048346, 9.82865019, 1.16525033, -13.33551209, 0.00220668, 0.07545536, 0.05739131, -0.04825559, 1.80652693, 2.67549075, -5.70109557, 29.08100962, -4.8972862, -3.62094093, 17.33689055, 3.88342235, -10.66310143, -16.12598707, 0.0152223, 0.3409704, 0.24479116, -0.2055869, -7.88839364, -67.36647678, 2.50063266, -4.8972862, 494.47921212, 170.86106875, -126.83263341, -160.48839026, 120.859385, 52.61261782, 0.02694279, 0.31291605, 0.18869954, -0.15785648, -3.41690543, 2.71012078, -13.4965165, -3.62094093, 170.86106875, 195.77662999, -137.85288085, -150.42863432, 105.91022589, 88.19349642, -0.00181838, -0.0353833, -0.02474101, 0.02076719, 2.72607278, -1.49637096, -1.26048346, 17.33689055, -126.83263341, -137.85288085, 182.19198937, 101.08256312, -138.37455048, -84.00691074, -0.02394434, -0.31640557, -0.19994578, 0.16745284, 2.53243921, 22.44235147, 9.82865019, 3.88342235, -160.48839026, -150.42863432, 101.08256312, 155.69932625, -103.48120197, -60.02029011, -0.00163328, 0.05690788, 0.05242226, -0.04422678, -1.48188089, -13.1893562, 1.16525033, -10.66310143, 120.859385, 105.91022589, -138.37455048, -103.48120197, 144.14781723, 59.2349263, 0.00522075, -0.11806133, -0.11383326, 0.09610551, -3.0286055, -0.05360568, -13.33551209, -16.12598707, 52.61261782, 88.19349642, -84.00691074, -60.02029011, 59.2349263, 99.27828035), 
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

    eUN.lmm_anova <- anova(eUN.lmm, effects = c("visit3-visit2=0"))
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
    
    
    e.lmm <- lmm(weight~visit+vita.time,
                 data=vitaminL,
                 repetition = ~visit|animal,
                 structure = "UN",
                 control = list(transform.sigma = "none", transform.k = "none", transform.rho = "none"))
    e.lmm2 <- lmm(weight~treatment*visit,
                 data=vitaminL,
                 repetition = ~visit|animal,
                 structure = "UN",
                 control = list(transform.sigma = "none", transform.k = "none", transform.rho = "none"))

    expect_equal(as.double(logLik(e.lmm2)), as.double(logLik(e.gls)), tol = 1e-6)
    expect_equal(as.double(logLik(e.lmm)), as.double(logLik(e.gls)), tol = 1e-6)
    expect_equal(as.double(logLik(e.gls)), -232.08183621, tol = 1e-6)

    ## here using a transformation seems worse ...
    expect_equal(as.double(logLik(e0.lmm)), -233.2156, tol = 1e-6)
    expect_true(all(abs(score(e0.lmm, p = coef(e.lmm, effects = "all"), effects = "all"))<1e-4))
    
    ## ** extract information
    autoplot(e.lmm, color = "group",ci =FALSE)
    qqtest::qqtest(residuals(e.lmm))
    confint(e.lmm)
    e.lmm_anova <- anova(e.lmm2, effects = "treatmentvitamin:visit6 - treatmentcontrol:visit6 = 0")
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
