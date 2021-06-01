### test-previous-bug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 23 2020 (12:33) 
## Version: 
## Last-Updated: Jun  1 2021 (10:16) 
##           By: Brice Ozenne
##     Update #: 15
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
data(gastricbypassL)
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

######################################################################
### test-previous-bug.R ends here
