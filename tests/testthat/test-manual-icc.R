### test-manual-icc.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 20 2022 (12:12) 
## Version: 
## Last-Updated: apr 20 2022 (13:35) 
##           By: Brice Ozenne
##     Update #: 4
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
    library(psych)   
    library(lme4) 

    library(LMMstar)
}

context("Compare ICC estimation")
LMMstar.options(optimizer = "FS", method.numDeriv = "simple", precompute.moments = TRUE)


if(FALSE){

    library(mvtnorm)
    set.seed(10)
    Y <- rmvnorm(100, mean = 0:1, sigma = 0.5 + diag(0.5,2,2))
    dfL <- melt(Y)

    ## ** correlation coefficient
    e0.lmm <- lmm(value ~ Var2,
                  repetition =~ Var2|Var1,
                  type.information = "observed",
                  data = dfL, structure = "UN")
    test <- model.tables(e0.lmm, effects = "correlation")[,c("estimate","lower","upper")]
    GS <- unlist(cor.test(Y[,1],Y[,2])[c("estimate","conf.int")])
    as.double(GS-test)

    
    ## ** ICC
    e.icc <- ICC(Y)

    e1.lmm <- lmm(value ~ 1,
                  repetition =~ Var2|Var1,
                  data = dfL, structure = "CS", df = FALSE)
    GS <- e.icc$results[e.icc$results$type=="ICC1",c("ICC","lower bound","upper bound")]
    test <- model.tables(e1.lmm, effects = "correlation")[,c("estimate","lower","upper")]
    as.double(GS-test)

    e3.lmm <- lmm(value ~ Var2,
                  repetition =~ Var2|Var1,
                  data = dfL, structure = "CS", df = FALSE)
    test <- model.tables(e3.lmm, effects = "correlation")[,c("estimate","lower","upper")]
    GS <- e.icc$results[e.icc$results$type=="ICC3",c("ICC","lower bound","upper bound")]
    as.double(GS-test)

    Z <- scale(Y)
    e.lm <- lm(Z[,2] ~ Z[,1], data = dfL)
    test2 <- c(coef(e.lm)[2],confint(e.lm)[2,])
    as.double(GS-test2)
    
    

    ## assess type 1 error
    ls.sim <- lapply(1:1000, function(iX){
        Y <- rmvnorm(100, mean = 0:1, sigma = diag(1,2,2))
        dfL <- melt(Y)
        e.icc <- ICC(Y)

        e.lmm <- lmm(value ~ Var2,
                     repetition =~ Var2|Var1,
                     data = dfL, structure = "CS", df = FALSE)
        c(e.icc$results$p[3], model.tables(e.lmm, effects="correlation")$p.value)
    })
    colMeans(do.call(rbind,ls.sim)<=0.05)

}


##----------------------------------------------------------------------
### test-manual-icc.R ends here
