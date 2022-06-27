### test-auto-mlmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:20) 
## Version: 
## Last-Updated: jun 27 2022 (11:59) 
##           By: Brice Ozenne
##     Update #: 54
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

## * Multiple imputation
set.seed(10)
n <- 100
X <- rnorm(n)
Y <- rnorm(n) + 0.25*X
df <- data.frame(Y=Y, X=X)
df[1:5,"X"] <- NA


test_that("mlmm: pool",{
    dfA <- mice(df, m = 10, printFlag = FALSE)
    GS <- summary(pool(with(dfA, lm(Y~X))))

    e.mlmm <- mlmm(Y~X, by = ".imp", data = complete(dfA, action = "long"),
                   effects = c("X=0"))
    test <- model.tables(e.mlmm, method = "pool")

    expect_equal(as.double(GS[GS$term=="X",c("estimate","std.error","df","p.value")]),
                 as.double(test[c("estimate","se","df","p.value")]), tol = 1e-4)
})


##----------------------------------------------------------------------
### test-auto-mlmm.R ends here
