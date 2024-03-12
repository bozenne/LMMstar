### test-auto-practical.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  7 2021 (17:03) 
## Version: 
## Last-Updated: mar 12 2024 (10:26) 
##           By: Brice Ozenne
##     Update #: 126
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
    gastricbypassL$time <- factor(gastricbypassL$time,
                                  levels = c("3monthsBefore", "1weekBefore", "1weekAfter", "3monthsAfter"),
                                  labels = c("m3B","w1B","w1A","m3A"))
    gastricbypassL$visit <- as.numeric(gastricbypassL$visit)

    ## ** summarize
    summarize(glucagonAUC~time, data = gastricbypassL, na.rm=TRUE)
    summarize(glucagonAUC~time|id, data = gastricbypassL, na.rm=TRUE)
    
    ## ** compound symmetry
    eCS.gls <- gls(glucagonAUC~time,
                   data=gastricbypassL,
                   correlation=corCompSymm(form=~visit|id),
                   na.action=na.exclude,
                   control=glsControl(opt="optim"))

    eCS.lmm <- lmm(glucagonAUC~time,
                   data=gastricbypassL,
                   repetition = ~time|id,
                   structure = "CS")
    
    expect_equal(as.double(logLik(eCS.gls)), as.double(logLik(eCS.lmm)), tol = 1e-6)

    ## ** unstructured with missing data
    eUN.gls <- gls(glucagonAUC~time,
                   data=gastricbypassL,
                   correlation=corSymm(form=~visit|id),
                   weights=varIdent(form=~1|time),
                   na.action=na.exclude,
                   control=glsControl(opt="optim"),
                   method = "REML")

    eUN.lmm <- lmm(glucagonAUC~time,
                   data=gastricbypassL,
                   repetition = ~time|id,
                   structure = "UN",
                   method.fit = "REML", trace = 0)

    ## check moments
    expect_equal(as.double(logLik(eUN.gls)), as.double(logLik(eUN.lmm)), tol = 1e-6)
    
    ## GS <- numDeriv::jacobian(func = function(p){logLik(eUN.lmm, p = p, transform.sigma = "none", transform.k = "none", transform.rho = "none")}, x = coef(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all"))
    GS <- rbind(c(0, 0, 0, 0, 0, 3e-08, 3e-08, -1.2e-07, -5.6e-07, -3.9e-07, 3e-08, 4.1e-07, -1.2e-07, 2.7e-07))
    expect_equal(as.double(GS), as.double(score(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all")), tol = 1e-6)

    ## GS <- -numDeriv::jacobian(func = function(p){score(eUN.lmm, p = p, transform.sigma = "none", transform.k = "none", transform.rho = "none",effects = "all")}, x = coef(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all"))
    GS <- cbind(c(3.75e-06, 1.79e-06, -3e-08, 1.43e-06, 0, 8.477e-05, -1.843e-05, 9.63e-06, 6.084e-05, 0.00011196, -7.35e-06, -9.894e-05, -7.91e-06, 2.433e-05), 
                c(1.79e-06, 6.88e-06, 1.06e-06, -8.4e-07, 0, 0.00027411, -0.00016572, 0.00035807, 0.00141665, 0.00130009, -0.00014701, -0.00131459, 0.00023644, -0.00049052), 
                c(-3e-08, 1.06e-06, 1.04e-06, -9.1e-07, 0, 2.78e-06, -9.7e-05, 0.00029724, 0.00111002, 0.00085567, -0.00011219, -0.00090667, 0.00023771, -0.00051618), 
                c(1.43e-06, -8.4e-07, -9.1e-07, 1.86e-06, 0, 1.38e-06, 8.377e-05, -0.00025882, -0.0009654, -0.00074127, 9.752e-05, 0.00078633, -0.00020768, 0.00045129), 
                c(0, 0, 0, 0, 1.035e-05, 0.01115006, 0.00590608, 0.00848313, -0.03243156, -0.01404792, 0.01120771, 0.01041161, -0.00609246, -0.01245153), 
                c(8.477e-05, 0.00027411, 2.78e-06, 1.38e-06, 0.01115006, 122.18555589, -1.45964668, 3.08816828, -68.07837646, 2.73876754, -1.51219029, 22.6795013, -13.32873148, -0.05416957), 
                c(-1.843e-05, -0.00016572, -9.7e-05, 8.377e-05, 0.00590608, -1.45964668, 21.44038778, -7.18200571, 2.75807973, -14.88594103, -1.39025149, 10.84047319, 1.28521521, -14.70839091), 
                c(9.63e-06, 0.00035807, 0.00029724, -0.00025882, 0.00848313, 3.08816828, -7.18200571, 37.93794573, -5.59355628, -4.13575786, 19.80175087, 4.43555263, -12.17911963, -18.41870059), 
                c(6.084e-05, 0.00141665, 0.00111002, -0.0009654, -0.03243156, -68.07837647, 2.75807974, -5.59355628, 494.47865121, 170.86077633, -126.83250435, -160.48812254, 120.8592697, 52.61254837), 
                c(0.00011196, 0.00130009, 0.00085567, -0.00074127, -0.01404792, 2.73876752, -14.88594103, -4.13575786, 170.86077634, 195.77649853, -137.8528703, -150.42850734, 105.91019906, 88.19349447), 
                c(-7.35e-06, -0.00014701, -0.00011219, 9.752e-05, 0.01120771, -1.51219029, -1.39025149, 19.80175087, -126.83250437, -137.85287029, 182.19199363, 101.08254329, -138.37453473, -84.00691973), 
                c(-9.894e-05, -0.00131459, -0.00090667, 0.00078633, 0.01041161, 22.67950128, 10.8404732, 4.43555263, -160.48812249, -150.42850737, 101.08254331, 155.69921616, -103.48119506, -60.02026461), 
                c(-7.91e-06, 0.00023644, 0.00023771, -0.00020768, -0.00609246, -13.3287315, 1.28521521, -12.17911962, 120.85926971, 105.91019907, -138.37453473, -103.48119506, 144.14781706, 59.2349059), 
                c(2.433e-05, -0.00049052, -0.00051618, 0.00045129, -0.01245153, -0.05416957, -14.70839091, -18.41870059, 52.61254837, 88.19349446, -84.00691974, -60.02026461, 59.23490589, 99.2783842)
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

    apply(residuals(eUN.lmm, type = "normalized", format = "wide"),2,sd,na.rm=TRUE)
    cor(residuals(eUN.lmm, type = "normalized", format = "wide")[,-1], use = "pairwise")
    qqtest::qqtest(na.omit(residuals(eUN.lmm, type = "normalized")))

    eUN.lmm_anova <- anova(eUN.lmm, effects = c("timew1A-timew1B=0"), ci = TRUE)
    capture.output(eUN.lmm_anova)
    capture.output(summary(eUN.lmm_anova))
    expect_equal(eUN.lmm_anova$multivariate$df.denom, 19.24699, tol = 1e-1) ## Richardson
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

    expect_equal(as.double(GS[,1]),as.double(ranef(eRI.lmm)$estimate), tol = 1e-4)
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

    vasscoreL.imputed <- fitted(fit.UN, type = "outcome", keep.newdata = TRUE)
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
