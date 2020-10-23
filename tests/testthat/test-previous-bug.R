### test-previous-bug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 23 2020 (12:33) 
## Version: 
## Last-Updated: okt 23 2020 (13:14) 
##           By: Brice Ozenne
##     Update #: 6
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
    library(repeated)
    library(data.table)
}

context("Previous bug")

## * from: Julie Lyng Forman <jufo@sund.ku.dk> date: Fri, 23 Oct 2020 10:05:40 +0000
data(gastricbypassL)
test_that("AUC - BuyseTest vs pROC",{
    g.summaries <- summarize(weight~time, data=gastricbypassL, na.rm=TRUE)
})
######################################################################
### test-previous-bug.R ends here
