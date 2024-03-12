### test-auto-cor.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 20 2022 (12:12) 
## Version: 
## Last-Updated: mar 12 2024 (18:17) 
##           By: Brice Ozenne
##     Update #: 72
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
LMMstar.options(optimizer = "FS", method.numDeriv = "simple", precompute.moments = TRUE)

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

    test2 <- partialCor(c(V1,V2)~1, data = data.frame(V1 = Y[,1], V2= Y[,2]))
    
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
    expect_equal(cor(Orthodont$age,Orthodont$distance), partialCor(e.lm, se = FALSE)[,"estimate"], tol = 1e-3)

    e.lm2 <- lmm(distance ~ Sex+age, data = Orthodont)
    GS <- lava::partialcor(c(age,distance)~Sex, data = Orthodont)
    expect_equal(GS[,"cor"], partialCor(e.lm2, se = FALSE)["age","estimate"], tol = 1e-3)

    ## mixed model
    e.lmm <- lmm(distance ~ Sex*age, repetition = ~1|Subject, structure = "CS", data = Orthodont)
    e.R2lmm <- suppressWarnings(partialCor(e.lmm, se = FALSE, R2 = TRUE))

    ## e.lmer <- lme4::lmer(distance ~ Sex*age + (1|Subject), data = Orthodont)
    ## library(r2glmm); setNames(r2beta(e.lmer, method = "kr")[2:4,"Rsq"],r2beta(e.lmer, method = "kr")[2:4,"Effect"])

    GS <- c("age" = 0.57834264, "Sex:age" = 0.07388639, "Sex" = 0.00431524)
    GS - attr(e.R2lmm, "R2")[names(GS)] ## some difference in age effect
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

    test.UN <- partialCor(c(X1,X2)~1, data = df.W, repetition = ~time|id, structure = "UN")
    ## eTopUN.lmm <- lmm(value ~ variable, repetition = ~variable+time|id, data = df.L,
    ##                   structure = TOEPLITZ(type = "UN"),
    ##                   control = list(optimizer = "FS"))
    ## confint(eTopUN.lmm, effects = "correlation")["rho(1,2,dt=0)",]
    GS.UN <- data.frame("estimate" = c(0.47476928), "se" = c(0.0478157), "df" = c(12.4303314), "lower" = c(0.36463913), "upper" = c(0.57179996), "p.value" = c(1.87e-06))
    expect_equivalent(test.UN["marginal",], GS.UN, tol = 1e-5)

    test.PEARSON <- partialCor(c(X1,X2)~1, data = df.W, repetition = ~time|id, structure = "PEARSON")
    ## eTopPEARSON.lmm <- lmm(value ~ variable, repetition = ~variable+time|id, data = df.L,
    ##                        structure = TOEPLITZ(list(~variable,~variable+time), type = "UN"),
    ##                        control = list(optimizer = "FS"))
    ## confint(eTopPEARSON.lmm, effects = "correlation")["rho(1,2,dt=0)",]
    GS.PEARSON <- data.frame("estimate" = c(0.47638071), "se" = c(0.04816023), "df" = c(12.52209613), "lower" = c(0.36546699), "upper" = c(0.57395625), "p.value" = c(1.88e-06))
    expect_equivalent(test.PEARSON["marginal",], GS.PEARSON, tol = 1e-5)

    test.HLAG <- partialCor(c(X1,X2)~1, data = df.W, repetition = ~time|id, structure = "HLAG")
    ## eTopHLAG.lmm <- lmm(value ~ variable, repetition = ~variable+time|id, data = df.L,
    ##                     structure = TOEPLITZ(type = "LAG"),
    ##                     control = list(optimizer = "FS"))
    ## confint(eTopHLAG.lmm, effects = "correlation")["rho(1,2,dt=0)",]
    GS.HLAG <- data.frame("estimate" = c(0.47549331), "se" = c(0.04800254), "df" = c(13.38669723), "lower" = c(0.36577561), "upper" = c(0.572176), "p.value" = c(1.16e-06))
    expect_equivalent(test.HLAG["marginal",], GS.HLAG, tol = 1e-5)

    test.LAG <- partialCor(c(X1,X2)~1, data = df.W, repetition = ~time|id, structure = "LAG")
    ## eTopLAG.lmm <- lmm(value ~ variable, repetition = ~variable+time|id, data = df.L,
    ##                     structure = TOEPLITZ(list(~variable,~variable+time), type = "LAG"),
    ##                     control = list(optimizer = "FS"))
    ## confint(eTopLAG.lmm, effects = "correlation")["rho(1,2,dt=0)",]
    GS.LAG <- data.frame("estimate" = c(0.47388305),"se" = c(0.04810496),"df" = c(13.80173648),"lower" = c(0.36429768),"upper" = c(0.57052434),"p.value" = c(9.8e-07))
    expect_equivalent(test.LAG["marginal",], GS.LAG, tol = 1e-5)
    
    test.HCS <- partialCor(c(X1,X2)~1, data = df.W, repetition = ~time|id, structure = "HCS")
    ## eTopHCS.lmm <- lmm(value ~ variable, repetition = ~variable+time|id, data = df.L,
    ##                    structure = TOEPLITZ(list(~variable+time,~variable+time),type = "CS"),
    ##                    control = list(optimizer = "FS"))
    ## confint(eTopHCS.lmm, effects = "correlation")["rho(1,2,dt=0)",]
    GS.HCS <- data.frame("estimate" = c(0.47549676),"se" = c(0.04802152),"df" = c(13.5961607),"lower" = c(0.3659086),"upper" = c(0.57207872),"p.value" = c(1.04e-06))
    expect_equivalent(test.HCS["marginal",], GS.HCS, tol = 1e-5)

    test.CS <- partialCor(c(X1,X2)~1, data = df.W, repetition = ~time|id, structure = "CS")
    ## eTopCS.lmm <- lmm(value ~ variable, repetition = ~variable+time|id, data = df.L,
    ##                   structure = TOEPLITZ(type = "CS"),
    ##                   control = list(optimizer = "FS"))
    ## confint(eTopCS.lmm, effects = "correlation")["rho(1,2,dt=0)",]
    GS.CS <- data.frame("estimate" = c(0.4732798),"se" = c(0.04808291),"df" = c(14.13096783),"lower" = c(0.36401686),"upper" = c(0.56969305),"p.value" = c(8.3e-07))
    expect_equivalent(test.CS["marginal",], GS.CS, tol = 1e-5)
})


##----------------------------------------------------------------------
### test-auto-cor.R ends here
