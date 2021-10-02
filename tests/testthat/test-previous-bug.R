### test-previous-bug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 23 2020 (12:33) 
## Version: 
## Last-Updated: okt  2 2021 (17:10) 
##           By: Brice Ozenne
##     Update #: 78
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
LMMstar.options(optimizer = "gls", method.numDeriv = "Richardson", precompute.moments = TRUE,
                columns.confint = c("estimate","se","df","lower","upper","p.value"))

## * from: Julie Lyng Forman <jufo@sund.ku.dk> date: Fri, 23 Oct 2020 10:05:40 +0000
data(gastricbypassL, package = "LMMstar")
test_that("summarize - gastricbypass example",{
    g.summaries <- summarize(weight~time, data=gastricbypassL, na.rm=TRUE)
    expect_s3_class(g.summaries,"data.frame")
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
    getCoef(e.gls, effects = "variance")
    GS <- data.frame("estimate" = c(480.4, 54.8, 91.1, 90.06463245, 57.41601691, 99.73654382, -12.3292649, 73.96796617, 55.72691236), 
                       "std.error" = c(8.91159604, 9.07876221, 9.78750704, 15.79765582, 20.53007571, 18.68081963, 14.95122806, 18.19375468, 19.20870004), 
                       "t.value" = c(53.9072909, 6.03606513, 9.30778385, 5.70113905, 2.79667828, 5.33898115, -0.82463225, 4.06556906, 2.90112877), 
                       "p.value" = c(0, 1.8e-07, 0, 6e-07, 0.00726378, 2.18e-06, 0.41342073, 0.00016606, 0.00547614), 
                       "lower" = c(462.50922443, 36.57362433, 71.45075971, 58.34951501, 16.20017011, 62.23323352, -42.3451077, 37.44247955, 17.16383791), 
                       "upper" = c(498.29077557, 73.02637567, 110.74924029, 121.77974989, 98.63186372, 137.23985412, 17.68657789, 110.49345279, 94.28998681))
    test <- getCoef(e.gls)
    expect_equivalent(GS, test, tol = 1e-5)
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
                   repetition = ~time|id,
                   structure = "CS")
    e.gls <- gls(glucagon~time,
                 control=glsControl(opt="optim"),
                 data=gastricbypassL,
                 na.action = na.omit,
                 correlation = corCompSymm(form=~time|id))
    expect_equal(as.double(logLik(e.gls)), logLik(eCS.lmm), tol = 1e-5)
    expect_equal(coef(e.gls), coef(eCS.lmm, effects = "mean"), tol = 1e-5)

    ## previously a bug when asking the confindence interval for sd and only for variance effects
    eUN.lmm <- lmm(glucagon~time2,
                   control=glsControl(opt="optim"),
                   data=gastricbypassL,
                   repetition = ~time|id,
                   structure = "UN")

    GS <- data.frame("estimate" = c(8.23783313, 8.07975612, 8.7153242, 8.40737066), 
                     "se" = c(0.16221473, 0.16265549, 0.16425097, 0.16223679), 
                     "df" = c(4.77147491, 11.80714451, 16.73482229, 17.18382852), 
                     "lower" = c(7.81476781, 7.72471717, 8.36836618, 8.06535968), 
                     "upper" = c(8.66089845, 8.43479507, 9.06228223, 8.74938164))


    test <- confint(eUN.lmm, transform.k = "logsd", effects = "variance", backtransform = FALSE,
                    columns = c("estimate","se","df","lower","upper"))
    ## test2 <- confint(eUN.lmm, transform.k = "logsd", effects = "all", backtransform = FALSE,
    ##                 columns = c("estimate","se","df","lower","upper"))

    expect_equal(test$estimate, GS$estimate, tol = 1e-5)
    expect_equal(test$se, GS$se, tol = 1e-5)
    expect_equal(test$df, GS$df, tol = 1e-1)
    expect_equal(test$lower, GS$lower, tol = 1e-2)
    expect_equal(test$upper, GS$upper, tol = 1e-2)
    
})

## * from: Julie Lyng Forman <jufo@sund.ku.dk> date: Wednesday, 07/07/21 11:00 PM (1/3)
test_that("lmm - error when predicting due to missing values in the covariates",{
    data(gastricbypassL, package = "LMMstar")
    e.lmm  <- lmm(weight ~ glucagon, repetition = ~ time|id , structure = "CS", data = gastricbypassL)
    summary(e.lmm, print = FALSE) ## was bugging at some point
    capture.output(summary(e.lmm, hide.fit = TRUE, hide.sd = TRUE, hide.cor = TRUE)) ## was bugging at some point
    expect_equal(NROW(predict(e.lmm, newdata = gastricbypassL)),NROW(gastricbypassL))

    set.seed(10)
    gastricbypassL$Gender <- factor(as.numeric(gastricbypassL$id) %% 2, levels = 0:1, labels = c("M","F"))
    e2.lmm  <- lmm(weight ~ Gender*glucagon, repetition =Gender ~ time|id , structure = "CS", data = gastricbypassL)
    getVarCov(e2.lmm)
    summary(e2.lmm, print = FALSE)
})

## * from: Julie Lyng Forman <jufo@sund.ku.dk> date: Wednesday, 07/07/21 11:00 PM (2/3)
data("gastricbypassL", package = "LMMstar")
dfres.R <- gastricbypassL[order(gastricbypassL$id),]

test_that("lmm - studentized and normalized residuals",{
    fit.main <- lmm(weight~time, 
                    repetition=~visit|id,
                    structure="UN",
                    data=dfres.R,
                    df=TRUE)
    

    dfres.R$fitted <- predict(fit.main, newdata = dfres.R)$estimate
    dfres.R$residual <- residuals(fit.main, type = 'response')
    dfres.R$rpearson <- residuals(fit.main, type = 'pearson')
    dfres.R$rstudent <- residuals(fit.main, type = 'studentized')
    dfres.R$normalized <- residuals(fit.main, type = 'normalized')
    dfres.R$rscaled <- residuals(fit.main, type = 'scaled')

    dfres.SAS <- data.frame("id" = c( 1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4,  4,  5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8,  9,  9,  9,  9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 20, 20, 20, 20), 
                            "visit" = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4), 
                            "time" = c("-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month", "-3 month", "-1 week", "+1 week", "+3 month"), 
                            "weight" = c(127.2, 120.7, 115.5, 108.1, 165.2, 153.4, 149.2, 132.0, 109.7, 101.6,  97.7,  87.1, 146.2, 142.4, 136.7, 123.0, 113.1, 105.6,  99.9,  87.7, 158.8, 143.6, 134.6, 108.7, 115.4, 108.8, 103.0,  88.9, 123.8, 115.1, 112.2,  97.6, 105.8, 102.0,  98.0,  86.6, 118.0, 105.3,  99.3,  90.9, 137.7, 131.9, 126.3, 105.4, 129.0, 120.7, 116.6, 100.0, 117.4, 108.5, 104.5,  95.4, 118.0, 113.4, 108.3,  93.6, 115.0, 109.7, 103.5,  94.1, 131.2, 127.3, 119.5, 103.5, 122.4, 113.9, 109.0,  99.4, 151.6, 143.0, 135.3, 118.5, 173.0, 162.2, 155.0, 148.0, 100.9,  95.7,  89.9,  78.8), 
                            "glucagon_auc" = c( 5032.50,  4942.50, 20421.00,  9249.45, 12142.50, 14083.50, 10945.50,  7612.50, 10321.35,  6202.50, 20121.00, 17704.50,  6693.00,  6631.50, 13090.50,  4551.00,  7090.50, NA, 19155.00, 12345.00, 10386.00,  7609.50, 11778.00,  8014.80, 10258.80,  8424.60, 29979.75, 11837.70, 16797.75, 16300.50, 11040.00,  6163.50,  2500.50,  2376.00, 13657.50, 11286.00, 11862.00,  8212.50, 22875.00, 12339.00,  5199.00,  5502.00,  7906.50,  6543.00,  5146.50,  5217.00, 12897.00, 11499.00,  3642.00,  5911.50, 26555.10, 23245.50,  4372.95,  4520.10, 12903.90, 12536.10,  5410.50,  7833.00, NA, 18148.50,  7015.50,  5497.50, 16290.00, 10536.00,  6673.50,  4857.00, 17560.50,  8434.95,  5485.50,  5010.00, 16269.00,  7441.50,  6879.00,  7953.00, 12036.00, 10362.00, 14299.50,  8739.00, 26638.50, 11410.50), 
                            "Pred" = c(128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365, 128.970, 121.240, 115.700, 102.365), 
                            "StdErrPred" = c(4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365, 4.532370, 4.228446, 4.086486, 3.813365), 
                            "DF" = c(19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19), 
                            "Alpha" = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05), 
                            "Lower" = c(119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354, 119.48364, 112.38976, 107.14689,  94.38354), 
                            "Upper" = c(138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465, 138.4564, 130.0902, 124.2531, 110.3465), 
                            "Resid" = c( -1.770,  -0.540,  -0.200,   5.735,  36.230,  32.160,  33.500,  29.635, -19.270, -19.640, -18.000, -15.265,  17.230,  21.160,  21.000,  20.635, -15.870, -15.640, -15.800, -14.665,  29.830,  22.360,  18.900,   6.335, -13.570, -12.440, -12.700, -13.465,  -5.170,  -6.140,  -3.500,  -4.765, -23.170, -19.240, -17.700, -15.765, -10.970, -15.940, -16.400, -11.465,   8.730,  10.660,  10.600,   3.035,   0.030,  -0.540,   0.900,  -2.365, -11.570, -12.740, -11.200,  -6.965, -10.970,  -7.840,  -7.400,  -8.765, -13.970, -11.540, -12.200,  -8.265,   2.230,   6.060,   3.800,   1.135,  -6.570,  -7.340,  -6.700,  -2.965,  22.630,  21.760,  19.600,  16.135,  44.030,  40.960,  39.300,  45.635, -28.070, -25.540, -25.800, -23.565), 
                            "ScaledResid" = c(-0.08732387,  0.40461942,  0.21188042,  1.29889096,  1.78742581, -0.47809719,  1.88053634, -0.45788462, -0.95069543, -0.68290924,  0.71119677,  0.09833756,  0.85005097,  1.94137195,  0.34442823,  0.17580169, -0.78295466, -0.36471018, -0.51971403,  0.02661820,  1.47167849, -1.91658820, -1.82531105, -1.88501707, -0.66948298,  0.03323736, -0.52919839, -0.37617999, -0.25506463, -0.50516419,  1.80013244, -0.74212007, -1.14310394,  0.79644955,  0.56441157, -0.18509653, -0.54121063, -2.14849561, -0.64644793,  1.09466009,  0.43069907,  0.96100412,  0.18960696, -1.64049655,  0.00148007, -0.20990181,  1.05413810, -0.96048682, -0.57081194, -0.76048611,  0.82886609,  0.59788227, -0.54121063,  0.84640117,  0.06265796, -0.59452875, -0.68921718,  0.50256394, -0.82512583,  0.75959274,  0.11001820,  1.47930296, -1.57887521, -0.31025415, -0.32413435, -0.47089165,  0.29795027,  0.68502551,  1.11646276,  0.31965248, -1.01479308, -0.12758494,  2.17224284,  0.11269252, -0.11671776,  2.46878558, -1.38484798,  0.13994872, -0.88962184,  0.07405488), 
                            "ScaledDep" = c( 6.27547785,  1.20149294, -0.61972831,  1.14820414,  8.15022752,  0.31877633,  1.04892761, -0.60857144,  5.41210629,  0.11396429, -0.12041196, -0.05234926,  7.21285268,  2.73824547, -0.48718050,  0.02511487,  5.57984705,  0.43216335, -1.35132276, -0.12406862,  7.83448021, -1.11971468, -2.65691978, -2.03570389,  5.69331874,  0.83011088, -1.36080712, -0.52686681,  6.10773709,  0.29170934,  0.96852371, -0.89280689,  5.21969777,  1.59332308, -0.26719716, -0.33578335,  5.82159109, -1.35162208, -1.47805666,  0.94397327,  6.79350078,  1.75787765, -0.64200178, -1.79118337,  6.36428178,  0.58697172,  0.22252937, -1.11117364,  5.79198978,  0.03638742, -0.00274264,  0.44719545,  5.82159109,  1.64327469, -0.76895077, -0.74521558,  5.67358453,  1.29943746, -1.65673456,  0.60890592,  6.47281992,  2.27617649, -2.41048394, -0.46094097,  6.03866736,  0.32598188, -0.53365846,  0.53433869,  7.47926448,  1.11652601, -1.84640182, -0.27827176,  8.53504456,  0.90956605, -0.94832649,  2.31809876,  4.97795373,  0.93682225, -1.72123058, -0.07663194), 
                            "StudentResid" = c(-0.08959240, -0.02929788, -0.01122802,  0.34502317,  1.83386018,  1.74485131,  1.88069333,  1.78287040, -0.97539293, -1.06557462, -1.01052179, -0.91835723,  0.87213389,  1.14804271,  1.17894209,  1.24142165, -0.80329454, -0.84855331, -0.88701357, -0.88226065,  1.50991027,  1.21314910,  1.06104788,  0.38111976, -0.68687504, -0.67493626, -0.71297926, -0.81006749, -0.26169078, -0.33312771, -0.19649035, -0.28666703, -1.17279990, -1.04387248, -0.99367976, -0.94843772, -0.55527039, -0.86482991, -0.92069763, -0.68974554,  0.44188792,  0.57836178,  0.59508505,  0.18258855,  0.00151852, -0.02929788,  0.05052609, -0.14228070, -0.58564069, -0.69121286, -0.62876911, -0.41902117, -0.55527039, -0.42536176, -0.41543674, -0.52731092, -0.70712191, -0.62610647, -0.68490921, -0.49723043,  0.11287630,  0.32878728,  0.21333238,  0.06828270, -0.33255483, -0.39823410, -0.37613866, -0.17837728,  1.14546663,  1.18059591,  1.10034595,  0.97069728,  2.22867413,  2.22229818,  2.20630590,  2.74544595, -1.42082405, -1.38568104, -1.44841456, -1.41769330), 
                            "PearsonResid" = c(-0.08732387, -0.02855604, -0.01094372,  0.33628699,  1.78742581,  1.70067068,  1.83307311,  1.73772711, -0.95069543, -1.03859366, -0.98493480, -0.89510391,  0.85005097,  1.11897362,  1.14909060,  1.20998815, -0.78295466, -0.82706746, -0.86455388, -0.85992131,  1.47167849,  1.18243148,  1.03418154,  0.37146959, -0.66948298, -0.65784650, -0.69492622, -0.78955612, -0.25506463, -0.32469272, -0.19151510, -0.27940846, -1.14310394, -1.01744104, -0.96851922, -0.92442274, -0.54121063, -0.84293193, -0.89738504, -0.67228079,  0.43069907,  0.56371733,  0.58001716,  0.17796530,  0.00148007, -0.02855604,  0.04924674, -0.13867807, -0.57081194, -0.67371096, -0.61284832, -0.40841132, -0.54121063, -0.41459136, -0.40491764, -0.51395911, -0.68921718, -0.61025310, -0.66756692, -0.48464028,  0.11001820,  0.32046220,  0.20793068,  0.06655375, -0.32413435, -0.38815058, -0.36661462, -0.17386067,  1.11646276,  1.15070255,  1.07248456,  0.94611868,  2.17224284,  2.16602833,  2.15044099,  2.67592970, -1.38484798, -1.35059482, -1.41173989, -1.38179650))

    expect_equal(dfres.R$residual, dfres.SAS$Resid, tol = 1e-8)
    expect_equal(dfres.R$rpearson, dfres.SAS$PearsonResid, tol = 1e-4)
    expect_equal(dfres.R$rstudent, dfres.SAS$StudentResid, tol = 1e-4)
    expect_equal(dfres.R$rscaled, dfres.SAS$ScaledResid, tol = 1e-4)
})

test_that("lmm - predicted values",{
    fit.main <- lmm(weight~time, 
                    repetition=~visit|id,
                    structure="UN",
                    data=dfres.R,
                    df=TRUE)
    set.seed(11)
    fit.main2 <- lmm(weight~time, 
                     repetition=~visit|id,
                     structure="UN",
                     data=dfres.R[sample.int(NROW(dfres.R),replace = FALSE),,drop=FALSE],
                     df=TRUE)
    ## check sensitivity to ordering of the values
    expect_equal(logLik(fit.main2),logLik(fit.main))

    ## error due to wrong factor
    expect_error(predict(fit.main, newdata = data.frame(time = "-1 week"), se = FALSE))
    ## valid prediction
    expect_equal(predict(fit.main, newdata = data.frame(time = "1 week before surgery"), se = FALSE)[[1]],
                 sum(coef(fit.main)[1:2]))
    expect_equal(predict(fit.main, newdata = data.frame(time = "1 week before surgery"), se = "estimation"),
                 data.frame("estimate" = c(121.24), 
                            "se" = c(4.22845441), 
                            "df" = c(18.99992887), 
                            "lower" = c(112.39029934), 
                            "upper" = c(130.08970066)),
                 tol = 1e-3)
    expect_equal(predict(fit.main, newdata = data.frame(time = "1 week before surgery", visit = 1:4, id = c(1,1,2,2)), se = "total"),
                 data.frame("estimate" = c(121.24, 121.24, 121.24, 121.24), 
                            "se" = c(20.70577588, 19.37721241, 18.75815531, 17.57030289), 
                            "df" = c(18.99992887, 18.99992887, 18.99992887, 18.99992887), 
                            "lower" = c(77.90230205, 80.68301804, 81.97871976, 84.4649241), 
                            "upper" = c(164.57769795, 161.79698196, 160.50128024, 158.0150759)),
                 tol = 1e-3)

    data("gastricbypassW", package = "LMMstar")
    GS <- predict(lm(weight2 ~ weight1, data = gastricbypassW), newdata = data.frame(weight1 = 50), se = TRUE)
    test <- predict(fit.main, newdata = data.frame(time = c("3 months before surgery","1 week before surgery"), visit = 1:2, weight = c(50,NA), id = c(1,1)),
                    type = "dynamic", keep.newdata = FALSE)
    expect_equivalent(test$estimate, GS$fit, tol = 1e-3)
    expect_equivalent(test,
                      data.frame("estimate" = c(48.3228695), 
                                 "se" = c(17.10624779), 
                                 "df" = c(Inf), 
                                 "lower" = c(14.79523992), 
                                 "upper" = c(81.85049907)),
                      tol = 1e-3)
    
})

## * from: Julie Lyng Forman <jufo@sund.ku.dk> date: Wednesday, 07/07/21 11:00 PM (3/3)
test_that("lmm - constrain model, cluster with 1 or several observations",{

    ## data
    data("calciumL", package = "LMMstar")
    calciumL$time <- 0.5 * (calciumL$visit-1)
    calciumL$timefac <- factor(calciumL$time)
    
    ## treatment variable
    calciumL$treat <- factor(calciumL$grp, c('N','P','C'))
    calciumL$treat[calciumL$time=="0"] <- "N"
    table(calciumL$time, calciumL$treat)
    table(calciumL$grp, calciumL$treat)

    calciumL$treat2 <- calciumL$treat
    calciumL$treat2[calciumL$treat=="N"] <- "C"
    calciumL$treat2 <- droplevels(calciumL$treat2)
    
    ## constrained time-treatment interaction
    calciumL$treat.time <- calciumL$timefac
    calciumL$treat.time[calciumL$grp=='P'] <- "0"

    ## Set reference points for time and treatment factors:
    calciumL$timfac <- relevel(calciumL$timefac, ref="0")
    calciumL$treat <- relevel(calciumL$treat, ref="N")
    calciumL$treat2 <- relevel(calciumL$treat2, ref="C")
    calciumL$treat.time <- relevel(calciumL$treat.time, ref="0")

    ## Fit the constrained linear mixed model:
    fit.clmm <- suppressWarnings(lmm(bmd~treat*timefac,
                                     repetition=~visit|girl,
                                     structure="UN",
                                     data=calciumL,
                                     df=FALSE))
    fit.clmm.bis <- suppressWarnings(lmm(bmd~0+treat:timefac,
                                         repetition=~visit|girl,
                                         structure="UN",
                                         data=calciumL,
                                         df=FALSE))
    ## coef(fit.clmm)
    ## coef(fit.clmm.bis)
    ## summary(fit.clmm)
    ## summary(fit.clmm.bis)
    ## autoplot(fit.clmm)

    fit.clmm2 <- suppressWarnings(lmm(bmd~treat2*timefac,
                                      repetition=~visit|girl,
                                      structure="UN",
                                      data=calciumL,
                                      df=FALSE))
    ## coef(fit.clmm2)
    ## summary(fit.clmm2)
    
    expect_equal(logLik(fit.clmm),logLik(fit.clmm.bis), tol = 1e-5)
    expect_equal(logLik(fit.clmm),logLik(fit.clmm2), tol = 1e-5)

})

## * from: Brice 08/23/21 16:45:12
test_that("glht - number of parameters",{
    
    set.seed(10)
    dL <- sampleRem(100, n.times = 3, format = "long")
 
    ## fit Linear Mixed Model
    eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)

    LMMstar.options(effects = c("mean","variance","correlation"))
    Mc <- matrix(0, nrow = 1, ncol = length(coef(eUN.lmm)), dimnames = list(NULL, names(coef(eUN.lmm))))
    Mc[,2] <- 1

    CI.glht <- multcomp::glht(eUN.lmm, linfct = Mc, rhs = 0, df = 10,
                              coef. = function(iX){coef.lmm(iX, effects = "all")},
                              vcov. = function(iX){vcov.lmm(iX, effects = "all")})
    expect_equal(NCOL(CI.glht$linfct),length(coef(eUN.lmm)))

    LMMstar.options(effects = c("mean"))
    Mc <- matrix(0, nrow = 1, ncol = length(coef(eUN.lmm)), dimnames = list(NULL, names(coef(eUN.lmm))))
    Mc[,2] <- 1

    CI.glht <- multcomp::glht(eUN.lmm, linfct = Mc, rhs = 0, df = 10)
    expect_equal(NCOL(CI.glht$linfct),length(coef(eUN.lmm)))

})

## * from: Brice 09/22/21 11:20 
test_that("lmm - estimation with missing data",{
    data(armd.wide, package = "nlmeU")

    armd.long <- reshape2::melt(armd.wide, 
                                id.var = c("subject","treat.f","lesion","miss.pat"),
                                measure.vars = c("visual0","visual4","visual12","visual24","visual52"),
                                variable.name = "week", 
                                value.name = "visual")
    armd.long <- armd.long[order(armd.long$subject),]
    armd.long$week <- factor(armd.long$week, 
                             level = c("visual0","visual4","visual12","visual24","visual52"), 
                             labels = c(0,4,12,24,52))
    rownames(armd.long) <- NULL

    e.lmm <- lmm(visual ~ week + week:treat.f,
                 repetition = ~ week | subject,
                 structure = "UN",
                 data = armd.long, df = FALSE)
    ## e.gls <- gls(visual ~ week + week:treat.f,
    ##              correlation = corSymm(form =~ as.numeric(week) | subject),
    ##              weights = varIdent(form =~1|week),
    ##              na.action = na.omit,
    ##              data = armd.long)

    expect_equal(as.double(logLik(e.lmm)),as.double(logLik(e.lmm$gls[[1]])))

    e2.lmm <- lmm(visual ~ 0 + week + week:treat.f,
                  repetition = treat.f ~ week | subject,
                  structure = "UN",
                  data = armd.long, df = FALSE)
    expect_equal(logLik(e2.lmm),sum(sapply(e2.lmm$gls,logLik)), tol = 1e-3)
    
    ## LMMstar.options(optimizer = "FS")
    ## e3.lmm <- lmm(visual ~ week + week:treat.f,
    ##               repetition = ~ week | subject,
    ##               structure = "UN",
    ##               data = armd.long)
    ## expect_equal(as.double(logLik(e3.lmm)),as.double(logLik(e.lmm)), tol = 1e-2)


})

## * from: Brice 09/22/21 11:20 
test_that("residuals.lmm - missing data",{
    data(ckdL, package = "LMMstar")

    ckdL$treat <- ckdL$allocation
    ckdL$treat[ckdL$time==0] <- "A"
    ckdL$treat <- relevel(ckdL$treat, ref="A")

    fit.main <- suppressWarnings(lmm(aix~time+time:treat, 
                                     repetition=~visit|id, 
                                     structure='UN',
                                     df=TRUE,
                                     data=ckdL))
    ## was leading to an error
    res <- residuals(fit.main, type = "all", keep.data = TRUE)
    GS <- as.character(unique(ckdL[rowSums(is.na(ckdL))>0,"id"]))
    test <- as.character(unique(res[rowSums(is.na(res))>0,"id"]))
    expect_equal(sort(GS),sort(test))
})

######################################################################
### test-previous-bug.R ends here
