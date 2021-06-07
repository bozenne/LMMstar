### test-practical.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  7 2021 (17:03) 
## Version: 
## Last-Updated: Jun  8 2021 (00:33) 
##           By: Brice Ozenne
##     Update #: 10
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

    library(LMMstar)
}

context("Check lmm on the example from the course with Julie")
LMMstar.options(method.numDeriv = "Richardson")

## * Practical 1
test_that("practical 1 - gastricbypass",{
    data(gastricbypassL)
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
                   structure = CS(~time|id))
    
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
                   structure = UN(~time|id),
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
                   structure = UN(~time|id),
                   method.fit = "REML", trace = 0)
    
    expect_equal(as.double(logLik(eUN.gls)), as.double(logLik(eUN.lmm)), tol = 1e-6)

    ## ** extract information
    backtransform(confint(eUN.lmm))[,c("estimate","lower","upper")]
    coef(eUN.lmm, transform.k = "sd")
    emmeans(eUN.lmm, specs = ~time)
    emmip(eUN.lmm, ~time)
    apply(residuals(eUN.lmm, type.residual = "normalized", format = "wide"),2,sd,na.rm=TRUE)
    cor(residuals(eUN.lmm, type.residual = "normalized", format = "wide")[,-1], use = "pairwise")
    qqtest(na.omit(residuals(eUN.lmm, type.residual = "normalized")))

    anova(eUN.lmm, effects = c("timew1A-timew1B=0"), ci = TRUE)
})

## * Practical 2
test_that("practical 2 - ncgs",{
    data(ncgsL)
    ncgsL$visit <- as.numeric(ncgsL$visit)
    ncgsL$highdose.time <- ncgsL$time
    ncgsL$highdose.time[ncgsL$group=="placebo"] <- "0"
    ncgsL$time <- relevel(as.factor(ncgsL$time), ref="0")
    ncgsL$highdose.time <- relevel(as.factor(ncgsL$highdose.time), ref="0")
    
    ncgsL$treatment <- factor(ncgsL$group, c("none","placebo","highdose"))
    ncgsL$treatment[ncgsL$time=="0"] <- "none"
    
    ## ** unstructured with missing data
    e.gls <- gls(cholest~time+highdose.time,
                 data=ncgsL,
                 correlation=corSymm(form=~visit|id),
                 weights=varIdent(form=~1|time),
                 na.action=na.exclude,
                 control=glsControl(opt='optim'))

    e.lmm <- lmm(cholest~time+highdose.time,
                 data=ncgsL,
                 structure = UN(~time|id),
                 control=glsControl(opt='optim'))

                 
    e2.lmm <- lmm(cholest~treatment*time,
                 data=ncgsL,
                 structure = UN(~time|id),
                 control=glsControl(opt='optim'))

    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lmm)), tol = 1e-6)
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e2.lmm)), tol = 1e-6)
    
    ## ** extract information
    backtransform(confint(e2.lmm))[,c("estimate","lower","upper")]
    anova(e2.lmm, effects = c("treatmenthighdose:time6-treatmentplacebo:time6=0","treatmenthighdose:time12-treatmentplacebo:time12=0"), ci = TRUE)
    autoplot(e2.lmm, color = "group") 
    autoplot(e2.lmm, color = "treatment", alpha = 0.5) 
})
## * Practical 3
if(FALSE){

}

## * Practical 4
if(FALSE){

}


##----------------------------------------------------------------------
### test-practical.R ends here
