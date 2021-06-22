### test-previous-bug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 23 2020 (12:33) 
## Version: 
## Last-Updated: jun 22 2021 (10:57) 
##           By: Brice Ozenne
##     Update #: 23
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
    library(data.table)
    library(nlme)

    library(LMMstar)
}

context("Previous bug")

## * from: Julie Lyng Forman <jufo@sund.ku.dk> date: Fri, 23 Oct 2020 10:05:40 +0000
data(gastricbypassL, package = "LMMstar")
test_that("summarize - gastricbypass example",{
    g.summaries <- summarize(weight~time, data=gastricbypassL, na.rm=TRUE)
})

## * from: Julie Lyng Forman <jufo@sund.ku.dk> date: Tuesday, Nov 3, 2020 12:13:07 PM
vitaminL <- data.frame("group" = c("C", "C", "C", "C", "C", "T", "T", "T", "T", "T", "C", "C", "C", "C", "C", "T", "T", "T", "T", "T", "C", "C", "C", "C", "C", "T", "T", "T", "T", "T", "C", "C", "C", "C", "C", "T", "T", "T", "T", "T", "C", "C", "C", "C", "C", "T", "T", "T", "T", "T", "C", "C", "C", "C", "C", "T", "T", "T", "T", "T"), 
                       "animal" = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), 
                       "visit" = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6), 
                       "weight" = c(455, 467, 445, 485, 480, 514, 440, 495, 520, 503, 460, 565, 530, 542, 500, 560, 480, 570, 590, 555, 510, 610, 580, 594, 550, 565, 536, 569, 610, 591, 504, 596, 597, 583, 528, 524, 484, 585, 637, 605, 436, 542, 582, 611, 562, 552, 567, 576, 671, 649, 466, 587, 619, 612, 576, 597, 569, 677, 702, 675), 
                       "time" = c("1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "4", "4", "4", "4", "4", "4", "4", "4", "4", "4", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "7", "7", "7", "7", "7", "7", "7", "7", "7", "7"), 
                       "vita.time" = c("1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "5", "5", "5", "5", "5", "1", "1", "1", "1", "1", "6", "6", "6", "6", "6", "1", "1", "1", "1", "1", "7", "7", "7", "7", "7"))
vitaminL$time <- as.factor(vitaminL$time)

test_that("getCoef - Non-positive definite approximate variance-covariance",{
    e.gls <- gls(weight ~ time + vita.time,
                 data = vitaminL,
                 correlation = corSymm(form =~ as.numeric(time)|animal),
                 weight = varIdent(form =~ 1|time),
                 na.action = na.exclude,
                 control = glsControl(opt = 'optim'))
    ## intervals(e.gls)
    getCoef(e.gls)
    getCoef(e.gls, effects = "variance")
})

## * from: Julie Lyng Forman <jufo@sund.ku.dk> date: Tuesday, 06/18/21 2:03 PM
test_that("lmm - error due to minus sign in levels of a categorical variable",{
    data(gastricbypassL, package = "LMMstar")
    gastricbypassL$time2 <- factor(gastricbypassL$time,
                                  levels = c("3 months before surgery", "1 week before surgery", "1 week after surgery", "3 months after surgery"),
                                  labels = c("-3 months", "-1 week", "+1 week", "+3 months"))

    eCS.lmm <- lmm(glucagon~time,
                   control=glsControl(opt="optim"),
                   data=gastricbypassL,
                   structure = CS(~time|id))
    e.gls <- gls(glucagon~time,
                 control=glsControl(opt="optim"),
                 data=gastricbypassL,
                 na.action = na.omit,
                 correlation = corCompSymm(form=~time|id))
    ## expect_equal(names(coef(eCS.lmm, effects = "mean")), c("(Intercept)", "time1 week before surgery", "time1 week after surgery", "time3 months after surgery"))
    expect_equal(coef(e.gls), coef(eCS.lmm, effects = "mean"), tol = 1e-5)

    ## previously a bug when asking the confindence interval for sd and only for variance effects
    eUN.lmm <- lmm(glucagon~time,
                   control=glsControl(opt="optim"),
                   data=gastricbypassL,
                   structure = UN(~time|id))

    GS <- data.frame("estimate" = c(8.23783313, 8.07975612, 8.7153242, 8.40737066), 
                     "se" = c(0.16221473, 0.16264897, 0.16418598, 0.16223679), 
                     "statistic" = c(NA, NA, NA, NA), 
                     "df" = c(4.77200085, 11.7952058, 16.73655236, 17.1791855), 
                     "lower" = c(7.81478265, 7.72469084, 8.36850623, 8.06535271), 
                     "upper" = c(8.66088361, 8.43482139, 9.06214217, 8.74938861), 
                     "null" = c(NA, NA, NA, NA), 
                     "p.value" = c(NA, NA, NA, NA))
    test <- confint(eUN.lmm, transform.k = "logsd", effects = "variance")
    
    expect_equal(test$estimate, GS$estimate, tol = 1e-5)
    expect_equal(test$se, GS$se, tol = 1e-5)
    expect_equal(test$df, GS$df, tol = 1e-1)
    expect_equal(test$lower, GS$lower, tol = 1e-2)
    expect_equal(test$upper, GS$upper, tol = 1e-2)
    
})


######################################################################
### test-previous-bug.R ends here
