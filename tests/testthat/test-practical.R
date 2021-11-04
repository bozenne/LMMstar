### test-practical.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  7 2021 (17:03) 
## Version: 
## Last-Updated: nov  4 2021 (16:49) 
##           By: Brice Ozenne
##     Update #: 64
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

    library(LMMstar)
}

context("Check lmm on the example from the course with Julie")
LMMstar.options(optimizer = "gls", method.numDeriv = "simple", precompute.moments = TRUE,
                columns.confint = c("estimate","se","df","lower","upper","p.value"))

test.practical <- FALSE

## * Practical 1
test_that("practical 1 - gastricbypass",{
    if(test.practical){
    data(gastricbypassL, package = "LMMstar")
    gastricbypassL$time <- factor(gastricbypassL$time,
                                  levels = c("3 months before surgery", "1 week before surgery", "1 week after surgery", "3 months after surgery"),
                                  labels = c("m3B","w1B","w1A","m3A"))
    gastricbypassL$visit <- as.numeric(gastricbypassL$visit)

    ## ** compound symmetry
    eCS.gls <- gls(glucagon~time,
                   data=gastricbypassL,
                   correlation=corCompSymm(form=~visit|id),
                   na.action=na.exclude,
                   control=glsControl(opt="optim"))

    eCS.lmm <- lmm(glucagon~time,
                   control=glsControl(opt="optim"),
                   data=gastricbypassL,
                   repetition = ~time|id,
                   structure = "CS")
    
    expect_equal(as.double(logLik(eCS.gls)), as.double(logLik(eCS.lmm)), tol = 1e-6)
    
    ## ** unstructured without missing data
    test.nna <- tapply(is.na(gastricbypassL$glucagon),gastricbypassL$id, sum)==0
    gastricbypassL.full <- gastricbypassL[gastricbypassL$id %in% names(which(test.nna)),]
    
    eUN.gls <- gls(glucagon~time,
                   data=gastricbypassL.full,
                   correlation=corSymm(form=~visit|id),
                   weights=varIdent(form=~1|time),
                   control=glsControl(opt="optim"),
                   method = "REML")

    eUN.lmm <- lmm(glucagon~time,
                   data=gastricbypassL.full,
                   control=glsControl(opt="optim"),
                   repetition = ~time|id,
                   structure = "UN",
                   method.fit = "REML", trace = 0)
    
    expect_equal(as.double(logLik(eUN.gls)), as.double(logLik(eUN.lmm)), tol = 1e-6)

    ## ** unstructured with missing data
    eUN.gls <- gls(glucagon~time,
                   data=gastricbypassL,
                   correlation=corSymm(form=~visit|id),
                   weights=varIdent(form=~1|time),
                   na.action=na.exclude,
                   control=glsControl(opt="optim"),
                   method = "REML")

    eUN.lmm <- lmm(glucagon~time,
                   data=gastricbypassL,
                   control=glsControl(opt="optim"),
                   repetition = ~time|id,
                   structure = "UN",
                   method.fit = "REML", trace = 0)

    ## check moments
    expect_equal(as.double(logLik(eUN.gls)), as.double(logLik(eUN.lmm)), tol = 1e-6)
    
    GS <- jacobian(func = function(p){logLik(eUN.lmm, p = p, transform.sigma = "none", transform.k = "none", transform.rho = "none")}, x = coef(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all"))
    expect_equal(as.double(GS), as.double(score(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all")), tol = 1e-6)

    GS <- -jacobian(func = function(p){score(eUN.lmm, p = p, transform.sigma = "none", transform.k = "none", transform.rho = "none",effects = "all")}, x = coef(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all"))
    expect_equal(as.double(GS), as.double(information(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none", effects = "all")), tol = 1e-6)

    ## ** extract information
    confint(eUN.lmm, effects = "all", backtransform = TRUE)[,c("estimate","lower","upper")]
    coef(eUN.lmm, transform.k = "sd")
    emmeans(eUN.lmm, specs = ~time)
    emmip(eUN.lmm, ~time)
    autoplot(eUN.lmm)
    apply(residuals(eUN.lmm, type = "normalized", format = "wide"),2,sd,na.rm=TRUE)
    cor(residuals(eUN.lmm, type = "normalized", format = "wide")[,-1], use = "pairwise")
    qqtest(na.omit(residuals(eUN.lmm, type = "normalized")))

    eUN.lmm_anova <- anova(eUN.lmm, effects = c("timew1A-timew1B=0"), ci = TRUE)
    expect_equal(eUN.lmm_anova$all$df.denom, 19.24699, tol = 1e-1) ## Richardson
    }
})

## * Practical 2
test_that("practical 2 - ncgs",{
    if(test.practical){
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
                 control=glsControl(opt='optim'),
                 df = FALSE)
    
    e2.lmm <- suppressWarnings(lmm(cholest~treatment*time,
                                   data=ncgsL,
                                   repetition = ~time|id,
                                   structure = "UN",
                                   control=glsControl(opt='optim')))
    
    e3.lmm <- suppressWarnings(lmm(cholest~treatment2*time,
                                   data=ncgsL,
                                   repetition = ~time|id,
                                   structure = "UN",
                                   control=glsControl(opt='optim')))

    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lmm)), tol = 1e-6)
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e2.lmm)), tol = 1e-6)
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e3.lmm)), tol = 1e-6)
    
    ## ** extract information
    confint(e2.lmm, effects = "all", backtransform = TRUE)[,c("estimate","lower","upper")]
    e2.lmm_anova <- anova(e2.lmm, effects = c("treatmenthighdose:time6-treatmentplacebo:time6=0","treatmenthighdose:time12-treatmentplacebo:time12=0"), ci = TRUE)
    expect_equal(e2.lmm_anova$all$df.denom, 100.0411, tol = 1e-1) ## Richardson
    autoplot(e2.lmm, color = "group") 
    autoplot(e2.lmm, color = "group", alpha = 0.25)
    }
})

test_that("practical 2 - vitamin",{
    if(test.practical){
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
                  type.information = "expected",
                  repetition = ~visit|animal,
                  structure = "UN",
                  control=glsControl(opt='optim'))

    
    e.lmm <- suppressWarnings(lmm(weight~treatment*visit,
                                  data=vitaminL,
                                  repetition = ~visit|animal,
                                  structure = "UN",
                                  control=glsControl(opt='optim')))

    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lmm)), tol = 1e-6)
    
    ## ** extract information
    autoplot(e.lmm, color = "group",ci =FALSE)
    qqtest(residuals(e.lmm))
    confint(e.lmm)
    e.lmm_anova <- anova(e.lmm, effects = "treatmentvitamin:visit6 - treatmentcontrol:visit6 = 0", ci = TRUE)
    ## expect_equal(e.lmm_anova$all$df.denom, 0.4972829, tol = 1e-1) ## Richardson
   
    }
})

## * Practical 3
test_that("practical 3 - swabsL",{
    if(test.practical){
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
    getVarCov(eUN.lmm)

    ## ** compound symmetry
    eCS.lme <- lme(swabs ~ crowding + name,
                   data = swabsL,
                   random =~ 1 | family)
    eCS.lmm <- lmm(swabs ~ crowding + name, data = swabsL, structure = "CS", repetition = ~name|family)
    expect_equal(as.double(logLik(eCS.lmm)), as.double(logLik(eCS.lme)), tol = 1e-6)

    ## summary(eCS.lmm)
    getVarCov(eCS.lmm)
    emmip(eCS.lmm, crowding~name)

    ## with interaction
    eCSI.lmm <- lmm(swabs ~ crowding * name, data = swabsL, structure = "CS", repetition = ~name|family)
    emmip(eCSI.lmm, crowding~name)
    }
})

## * Practical 4
test_that("practical 4 - bloodpressureL",{
    if(test.practical){
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
        


        eUNT.gls <- gls(duration ~ period + treatment,
                       data = bloodpressureL,
                       correlation = corSymm(form=~ treatment.num | id),
                       weights = varIdent(form=~ 1|treatment),
                       )
        eUNT.lmm <- lmm(duration ~ period + treatment,
                         data = bloodpressureL,
                         structure = "UN", repetition = ~ treatment | id)
        expect_equal(as.double(logLik(eUNT.lmm)),as.double(logLik(eUNT.gls)), tol = 1e-6)
    }
})


##----------------------------------------------------------------------
### test-practical.R ends here
