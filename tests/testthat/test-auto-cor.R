### test-auto-cor.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 20 2022 (12:12) 
## Version: 
## Last-Updated: aug 30 2022 (11:38) 
##           By: Brice Ozenne
##     Update #: 59
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

test_that("estimate partial correlation via lmm (independence)", {

    ## univariate linear model
    e.lm <- lmm(distance ~ age, data = Orthodont)
    expect_equal(cor(Orthodont$age,Orthodont$distance), confint(e.lm, column = "partial.r")["age","partial.r"], tol = 1e-3)

    e.lm2 <- lmm(distance ~ Sex+age, data = Orthodont)
    GS <- lava::partialcor(c(age,distance)~Sex, data = Orthodont)
    expect_equal(GS[,"cor"], confint(e.lm2, column = "partial.r")["age","partial.r"], tol = 1e-3)

    ## mixed model
    e.lmm <- lmm(distance ~ Sex*age, repetition = ~1|Subject, structure = "CS", data = Orthodont)
    e.aovlmm <- anova(e.lmm)

    ## e.lmer <- lme4::lmer(distance ~ Sex*age + (1|Subject), data = Orthodont)
    ## library(r2glmm); setNames(r2beta(e.lmer, method = "kr")[2:4,"Rsq"],r2beta(e.lmer, method = "kr")[2:4,"Effect"])

    GS <- c("age" = 0.57834264, "Sex:age" = 0.07388639, "Sex" = 0.00431524)
    GS - e.aovlmm$multivariate[match(names(GS),e.aovlmm$multivariate),"partial.r2"] ## some difference in age effect
    
    expect_equal(e.aovlmm$multivariate$partial.r2, e.aovlmm$univariate$partial.r^2, tol = 1e-6)


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

## * Partial correlation with repeated measuremnts
set.seed(10)

n.time <- 3
n.id <- 100
sd.id <- 1.5
Sigma <- matrix(c(1,0.8,0.8,1),2,2)
##      [,1] [,2]
## [1,]  1.0  0.8
## [2,]  0.8  1.0
df.W <- data.frame(id = unlist(lapply(1:n.id, rep, n.time)),
                   time = rep(1:n.time,n.id),
                   rmvnorm(n.time*n.id, mean = c(3,3), sigma = Sigma)
                   )
## df.W$time2 <- as.factor(df.W$time)
df.W$X2 <- df.W$X2 + rnorm(n.id, sd = sd.id)[df.W$id]
df.W$id <- as.factor(df.W$id)
df.L <- reshape2::melt(df.W, id.vars = c("id","time")) 
df.L$time2 <- as.factor(as.numeric(as.factor(paste(df.L$variable,df.L$time,sep="."))))

Sigma.GS <- as.matrix(bdiag(Sigma,Sigma,Sigma))[c(1,3,5,2,4,6),c(1,3,5,2,4,6)]
Sigma.GS[4:6,4:6] <- Sigma.GS[4:6,4:6] + sd.id^2
cov2cor(Sigma.GS)
##           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
## [1,] 1.0000000 0.0000000 0.0000000 0.4437602 0.0000000 0.0000000
## [2,] 0.0000000 1.0000000 0.0000000 0.0000000 0.4437602 0.0000000
## [3,] 0.0000000 0.0000000 1.0000000 0.0000000 0.0000000 0.4437602
## [4,] 0.4437602 0.0000000 0.0000000 1.0000000 0.6923077 0.6923077
## [5,] 0.0000000 0.4437602 0.0000000 0.6923077 1.0000000 0.6923077
## [6,] 0.0000000 0.0000000 0.4437602 0.6923077 0.6923077 1.0000000
## cor(dcast(df.L, id~time2)[,-1])
##              1           2           3           4           5          6
## 1  1.000000000 -0.05237788 -0.02684335  0.44774229 0.004389268 0.06020988
## 2 -0.052377883  1.00000000  0.06204366 -0.01202896 0.504596466 0.03584304
## 3 -0.026843354  0.06204366  1.00000000 -0.05860345 0.019500348 0.47095753
## 4  0.447742292 -0.01202896 -0.05860345  1.00000000 0.652024959 0.70171094
## 5  0.004389268  0.50459647  0.01950035  0.65202496 1.000000000 0.68454889
## 6  0.060209882  0.03584304  0.47095753  0.70171094 0.684548892 1.00000000

## ggplot(df.W, aes(x = X1, y = X2, group = id, color = id)) + geom_point() + guides(color = "none") + geom_smooth(method="lm", se = FALSE)


test_that("estimate partial correlation via lmm (cluster)", {

    ## eWrong.lmm <- lmm(value ~ variable, repetition = ~time+variable|id, data = df.L,
    ##                   structure = CS(~variable, heterogeneous = TRUE), control = list(optimizer = "FS"))
    ## eOK.lmm <- lmm(value ~ variable, repetition = ~time2|id, data = df.L,
    ##                   structure = "UN", control = list(optimizer = "FS", trace = 2))

    test.hetero <- partialCor(c(X1,X2)~1, data = df.W, repetition = ~time|id, heterogeneous = TRUE)
    test.homo <- partialCor(c(X1,X2)~1, data = df.W, repetition = ~time|id, heterogeneous = 0.5)
    
    eTopHetero.lmm <- lmm(value ~ variable, repetition = ~time+variable|id, data = df.L,
                        structure = TOEPLITZ(heterogeneous = TRUE),
                        control = list(optimizer = "FS"))

    eTopHomo.lmm <- lmm(value ~ variable, repetition = ~time+variable|id, data = df.L,
                        structure = TOEPLITZ(heterogeneous = 0.5),
                        control = list(optimizer = "FS"))

    expect_equal(as.double(model.tables(eTopHetero.lmm, effects = "correlation")["rho(1.X1,1.X2)",]),
                 as.double(test.hetero), tol = 1e-6)
    expect_equal(c(4.332644e-01, 4.888257e-02, 1.978815e+01, 3.288053e-01, 5.272503e-01, 1.494162e-07),
                 as.double(test.hetero), tol = 1e-6)
    expect_equal(as.double(model.tables(eTopHomo.lmm, effects = "correlation")["rho(1.X1,1.X2)",]),
                 as.double(test.homo), tol = 1e-6)
    expect_equal(c(4.334200e-01, 4.875447e-02, 2.002694e+01, 3.293324e-01, 5.271000e-01, 1.314180e-07),
                 as.double(test.homo), tol = 1e-6)
    
})


##----------------------------------------------------------------------
### test-auto-cor.R ends here
