### test-auto-mlmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:20) 
## Version: 
## Last-Updated: jul 17 2025 (17:06) 
##           By: Brice Ozenne
##     Update #: 58
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
    library(mice)

    library(LMMstar)
}

context("Check mlmm ")
LMMstar.options(method.numDeriv = "Richardson", precompute.moments = TRUE)


##----------------------------------------------------------------------
### test-auto-mlmm.R ends here
