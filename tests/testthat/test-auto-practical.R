### test-auto-practical.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  7 2021 (17:03) 
## Version: 
## Last-Updated: May 30 2022 (23:35) 
##           By: Brice Ozenne
##     Update #: 110
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
    library(emmeans)
    library(nlme)
    library(qqtest)
    library(ggplot2)
    library(lme4)
    library(lmerTest)

    library(LMMstar)
}

context("Check lmm on the example from the course with Julie")
LMMstar.options(optimizer = "gls", method.numDeriv = "simple", precompute.moments = TRUE,
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
    
    GS <- numDeriv::jacobian(func = function(p){logLik(eUN.lmm, p = p, transform.sigma = "none", transform.k = "none", transform.rho = "none")}, x = coef(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all"))
    expect_equal(as.double(GS), as.double(score(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all")), tol = 1e-6)

    GS <- -numDeriv::jacobian(func = function(p){score(eUN.lmm, p = p, transform.sigma = "none", transform.k = "none", transform.rho = "none",effects = "all")}, x = coef(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all"))
    expect_equal(as.double(GS), as.double(information(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all")), tol = 1e-6)

    ## ** extract information
    capture.output(model.tables(eUN.lmm))
    capture.output(summary(eUN.lmm, columns = "statistic"))
    capture.output(confint(eUN.lmm, effects = "all", backtransform = TRUE))
    expect_equal(unname(coef(eUN.lmm, effects = "variance", transform.k = "sd")),
                 unname(coef(eUN.lmm, effects = "variance")[1]*c(1,coef(eUN.lmm, effects = "variance")[-1])),
                 tol = 1e-6)
    ## emmeans(eUN.lmm, specs = ~time)
    ## emmip(eUN.lmm, ~time)
    autoplot(eUN.lmm)

    apply(residuals(eUN.lmm, type = "normalized", format = "wide"),2,sd,na.rm=TRUE)
    cor(residuals(eUN.lmm, type = "normalized", format = "wide")[,-1], use = "pairwise")
    qqtest::qqtest(na.omit(residuals(eUN.lmm, type = "normalized")))

    eUN.lmm_anova <- anova(eUN.lmm, effects = c("timew1A-timew1B=0"), ci = TRUE)
    capture.output(eUN.lmm_anova)
    capture.output(summary(eUN.lmm_anova))
    expect_equal(eUN.lmm_anova$all$df.denom, 19.24699, tol = 1e-1) ## Richardson
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
    expect_equal(e2.lmm_anova$all$df.denom, 100.0411, tol = 1e-1) ## Richardson
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

    e0.lmm <- lmm(weight~visit+vita.time,
                  data=vitaminL,
                  repetition = ~visit|animal,
                  structure = "UN", control = list(trace = 3))

    
    e.lmm <- suppressWarnings(lmm(weight~treatment*visit,
                                  data=vitaminL,
                                  repetition = ~visit|animal,
                                  structure = "UN"))

    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lmm)), tol = 1e-6)
    
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
    ## emmip(eCS.lmm, crowding~name)

    ## with interaction
    eCSI.lmm <- lmm(swabs ~ crowding * name, data = swabsL, structure = "CS", repetition = ~name|family)
    emmip(eCSI.lmm, crowding~name)

    eCSI.lme <- lme(swabs ~ crowding * name, random =~ 1 | family, data = swabsL)
    GS <- as.data.frame(ranef(eCSI.lme, augFrame = TRUE))

    expect_equal(as.double(GS[,1]),as.double(coef(eCSI.lmm, effects = "ranef")), tol = 1e-4)
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
    emmip(eCS.lmm, treatment~period)

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

    expect_equal(sigma(fit.CS), GS, tol = 1e-5)
    
    GS <- lmer(vas~-1+treatment+(1|id), data=vasscoresL)    
    expect_equal(unname(coef(fit.CS)), unname(fixef(GS)), tol = 1e-5)
    expect_equal(unname(model.tables(fit.CS)$se), unname(summary(GS)$coef[,"Std. Error"]), tol = 1e-2)
    expect_equal(unname(model.tables(fit.CS)$df), unname(summary(GS)$coef[,"df"]), tol = 1e-2)

    ## autoplot(fit.CS)
    suppressWarnings(autoplot(fit.CS, obs.alpha = 0.1))
    dummy.coef(fit.CS)
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
    expect_equal(sigma(fit.UN), GS, tol = 1e-5)

    vasscoreL.imputed <- fitted(fit.UN, impute = TRUE, keep.newdata = TRUE)
    set.seed(11)
    gg <- ggplot(vasscoreL.imputed, aes(x=treatment, y = vas, group = id))
    gg + geom_jitter(aes(color = imputed), height = 1, width = 0) + geom_line()

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

    expect_equal(sigma(fit.CS.red), GS, tol = 1e-5)

    fit.UN.red <- lmm(vas~-1+treatment, data=vasscoreL.red,
                      repetition=~treatment|id, structure="UN", control = list(optimizer = "FS")) ## FS optimizer mandatory, GLS does something weird as it is able to estimate the correlation (B,C)
    capture.output(summary(fit.UN.red))
    
    GS <- matrix(c(1996.68088251, 1926.29216631, 1573.59329163, 1926.29216631, 2231.4711806, NA, 1573.59329163, NA, 1753.88563313), 
                 nrow = 3, 
                 ncol = 3, 
                 dimnames = list(c("A", "B", "C"),c("A", "B", "C")) 
                 ) 

    expect_equal(sigma(fit.UN.red), GS, tol = 1e-5)
})
##----------------------------------------------------------------------
### test-auto-practical.R ends here
