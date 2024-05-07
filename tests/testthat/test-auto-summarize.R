### test-auto-summarize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 10 2024 (15:10) 
## Version: 
## Last-Updated: maj  7 2024 (12:23) 
##           By: Brice Ozenne
##     Update #: 8
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
    library(LMMstar)
}

context("Check summarize function")

## * extract columns
set.seed(10)
dL <- sampleRem(1e2, n.time = 2, format = "long")

## ** full data
test_that("summarize full data, extract columns",{
    summarize(Y ~ 1, data = dL)

    aa <- summarize(Y ~ 1, data = dL)[,c("mean")]
    bb <- summarize(Y ~ 1, data = dL, columns = c("mean"))
    expect_equal(aa, bb$mean, tol = 1e-6)

    summarize(Y ~ 1, data = dL)[,c("mean","sd")]
    summarize(Y ~ 1, data = dL, columns = c("mean","sd"))
    summarize(Y ~ 1, data = dL, columns = c("mean","sd"))[,2:3]

    summarize(Y ~ visit, data = dL)

    summarize(Y ~ visit, data = dL)[,c("mean")]
    summarize(Y ~ visit, data = dL, columns = c("mean"))

    summarize(Y ~ visit, data = dL)[,c("mean","sd")]
    summarize(Y ~ visit, data = dL, columns = c("mean","sd"))
})


## ** missing data
dL.NA <- dL
dL.NA[1:5,"Y"] <- NA

test_that("summarize full data, extract columns",{
    summarize(Y ~ 1, data = dL.NA, na.rm = TRUE)

    aa <- summarize(Y ~ 1, data = dL.NA, na.rm = TRUE)[,c("mean")]
    bb <- summarize(Y ~ 1, data = dL.NA, columns = c("mean"), na.rm = TRUE)
    expect_equal(aa, bb$mean, tol = 1e-6)

    summarize(Y ~ 1, data = dL.NA, na.rm = TRUE)[,c("mean","sd")]
    summarize(Y ~ 1, data = dL.NA, columns = c("mean","sd"), na.rm = TRUE)

    summarize(Y ~ visit, data = dL.NA, na.rm = TRUE)

    summarize(Y ~ visit, data = dL.NA, na.rm = TRUE)[,c("mean")]
    summarize(Y ~ visit, data = dL.NA, columns = c("mean"), na.rm = TRUE)

    summarize(Y ~ visit, data = dL.NA, na.rm = TRUE)[,c("mean","sd")]
    summarize(Y ~ visit, data = dL.NA, columns = c("mean","sd"), na.rm = TRUE)
})
##----------------------------------------------------------------------
### test-auto-summarize.R ends here
