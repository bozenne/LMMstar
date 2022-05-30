### test-auto-cor.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 20 2022 (12:12) 
## Version: 
## Last-Updated: May 30 2022 (23:19) 
##           By: Brice Ozenne
##     Update #: 24
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
    library(mvtnorm)

    library(LMMstar)
}

context("Compare ICC estimation")
LMMstar.options(optimizer = "gls", method.numDeriv = "simple", precompute.moments = TRUE)

## * Correlation

## simulate data
set.seed(10)
Y <- mvtnorm::rmvnorm(100, mean = 0:1, sigma = 0.5 + diag(0.5,2,2))
dfL <- reshape2::melt(Y)

test_that("estimate correlation via lmm", {
    e0.lmm <- lmm(value ~ Var2,
                  repetition =~ Var2|Var1,
                  type.information = "observed",
                  data = dfL, structure = "UN")
    test <- model.tables(e0.lmm, effects = "correlation")[,c("estimate","lower","upper")]
    GS <- unlist(cor.test(Y[,1],Y[,2])[c("estimate","conf.int")])

    test2 <- partialCor(c(V1,V2)~1, data.frame(V1 = Y[,1], V2= Y[,2]))
    
    ## same point estimate
    expect_equal(as.double(test[,"estimate"]),as.double(GS["estimate.cor"]), tol = 1e-5)
    expect_equal(as.double(test2[,"estimate"]),as.double(GS["estimate.cor"]), tol = 1e-5)
 
    ## but different ci
    as.double(GS-test)
    expect_equal(as.double(test$lower),as.double(test2$lower), tol = 1e-5)
    expect_equal(as.double(test$upper),as.double(test2$upper), tol = 1e-5)
})

## * Partial correlation
data("Orthodont", package = "nlme")

test_that("estimate partial correlation via lmm", {

    ## univariate linear model
    e.lm <- lmm(distance ~ age, data = Orthodont)
    expect_equal(cor(Orthodont$age,Orthodont$distance), confint(e.lm, column = "partial.R")["age","partial.R"], tol = 1e-3)

    e.lm2 <- lmm(distance ~ age+Sex, data = Orthodont)
    GS <- lava::partialcor(c(age,distance)~Sex, data = Orthodont)
    expect_equal(GS[,"cor"], confint(e.lm2, column = "partial.R")["age","partial.R"], tol = 1e-3)

    ## mixed model
    e.lmm <- lmm(distance ~ age*Sex, repetition = ~1|Subject, structure = "CS", data = Orthodont)
    e.aovlmm <- summary(anova(e.lmm), columns = "partial.R", print = FALSE)
    ## e.lmer <- lme4::lmer(distance ~ age*Sex + (1|Subject), data = Orthodont)

    ## library(r2glmm); setNames(r2beta(e.lmer, method = "kr")[2:4,"Rsq"],r2beta(e.lmer, method = "kr")[2:4,"Effect"])
    GS <- c("age" = 0.57834264, "age:Sex" = 0.07388639, "Sex" = 0.00431524)
    GS - e.aovlmm[[1]][names(GS),"partial.R2"] ## some difference in age effect
    
    expect_equal(e.aovlmm[[1]][names(GS),"partial.R2"], e.aovlmm[[2]][names(GS),"partial.R"]^2, tol = 1e-6)
})

## * ICC
test_that("ICC", {
    ## library(psych)
    ## e.icc <- psych::ICC(Y)

    e1.lmm <- lmm(value ~ 1,
                  repetition =~ Var2|Var1,
                  data = dfL, structure = "CS", df = FALSE)
    test1 <- model.tables(e1.lmm, effects = "correlation")[,c("estimate","lower","upper")]
    ## GS1 <- e.icc$results[e.icc$results$type=="ICC1",c("ICC","lower bound","upper bound")]
    GS1 <- data.frame("ICC" = c(0.25796171), 
                     "lower bound" = c(0.09799708), 
                     "upper bound" = c(0.40506582))
    expect_equal(as.double(GS1["ICC"]),as.double(test1["rho","estimate"]), tol = 1e-6)

    e3.lmm <- lmm(value ~ Var2,
                  repetition =~ Var2|Var1,
                  data = dfL, structure = "CS", df = FALSE)
    test3 <- model.tables(e3.lmm, effects = "correlation")[,c("estimate","lower","upper")]
    ## GS3 <- e.icc$results[e.icc$results$type=="ICC3",c("ICC","lower bound","upper bound")]
    GS3 <- data.frame("ICC" = c(0.63110272), 
                     "lower bound" = c(0.52058073), 
                     "upper bound" = c(0.72082364))
    expect_equal(as.double(GS3["ICC"]),as.double(test3["rho","estimate"]), tol = 1e-6)


    ## ## assess type 1 error
    ## ls.sim <- lapply(1:1000, function(iX){
    ##     Y <- rmvnorm(100, mean = 0:1, sigma = diag(1,2,2))
    ##     dfL <- melt(Y)
    ##     e.icc <- ICC(Y)

    ##     e.lmm <- lmm(value ~ Var2,
    ##                  repetition =~ Var2|Var1,
    ##                  data = dfL, structure = "CS", df = FALSE)
    ##     c(e.icc$results$p[3], model.tables(e.lmm, effects="correlation")$p.value)
    ## })
    ## colMeans(do.call(rbind,ls.sim)<=0.05)

})


##----------------------------------------------------------------------
### test-auto-cor.R ends here
