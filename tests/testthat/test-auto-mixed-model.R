### test-auto-mixed-model.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 14 2021 (16:46) 
## Version: 
## Last-Updated: jul 31 2023 (18:06) 
##           By: Brice Ozenne
##     Update #: 199
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
    library(numDeriv)
    library(lava)
    library(multcomp)
    library(emmeans)
    library(lme4)
    library(lmerTest)

    library(LMMstar)
}

context("Check lmm on mixed model parametrized with random effects")
LMMstar.options(optimizer = "FS", method.numDeriv = "simple", precompute.moments = TRUE, # "Richardson"
                columns.confint = c("estimate","se","df","lower","upper","p.value"))

## * Random intercept model
data(Orthodont,package="nlme")

test_that("Random intercept model",{
    
    ## ** fit
    eRI.lmer <- lmer(distance ~ age + (1|Subject), data=Orthodont)
    eRI.lmm <- lmm(distance ~ age + (1|Subject), data=Orthodont,
                   control = list(init = "lmer"))
    eRI2.lmm <- lmm(distance ~ age + (1|Subject), data=Orthodont)
    
    xx <- capture.output(summary(eRI.lmm))

    ## ** iteration
    expect_equal(eRI.lmm$opt$n.iter,0)
    expect_equal(eRI2.lmm$opt$n.iter,4)

    ## ** likelihood
    expect_equal(as.double(logLik(eRI.lmer)), as.double(logLik(eRI.lmm)), tol = 1e-6)
    expect_equal(as.double(logLik(eRI.lmer)), as.double(logLik(eRI2.lmm)), tol = 1e-6)

    ## ** random effects (conditional mean)
    u.GS <- as.data.frame(ranef(eRI.lmer))    
    u.test <- ranef(eRI.lmm, effects = "mean", format = "long")
    expect_equal(as.double(u.GS$condval[match(u.test$level,u.GS$grp)]), as.double(u.test$estimate), tol = 1e-6)

    ## ** random effects (conditional variance)
    tau.GS <- as.data.frame(VarCorr(eRI.lmer))
    tau.test <- ranef(eRI.lmm, effects = "variance", format = "long")
    expect_equal(as.double(tau.GS[1,"vcov"]), as.double(tau.test[tau.test$type=="variance","estimate"]), tol = 1e-6)

})

## * Stratified random intercept model
Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")

test_that("Stratified random intercept model",{

    ## ** fit
    eRI.mlmm <- mlmm(distance ~ 1 + (1|Subject), data = Orthodont, by = "nsex", trace = FALSE)
    eSRI.lmm <- lmm(distance ~ nsex + (1|Subject), repetition = nsex~1|Subject, data = Orthodont)
    
    ## ** iteration
    expect_equal(eSRI.lmm$opt$n.iter,7)

    ## ** likelihood
    expect_equal(as.double(sum(unlist(logLik(eRI.mlmm)))), as.double(logLik(eSRI.lmm)), tol = 1e-6)

    ## ** random effects (conditional mean)
    u.GS <- do.call(rbind,ranef(eRI.mlmm))
    u.test <- ranef(eSRI.lmm, effects = "mean", format = "long")
    expect_equal(as.double(u.GS[match(u.test$level,u.GS$level),"estimate"]), as.double(u.test$estimate), tol = 1e-6)

    ## ** random effects (conditional variance)
    tau.GS <- do.call(rbind,ranef(eRI.mlmm, effects = "variance"))
    tau.test <- ranef(eSRI.lmm, effects = "variance", format = "long")
    expect_equal(as.double(tau.GS$estimate), as.double(tau.test$estimate), tol = 1e-5)
})

## * Crossed random intercept model (2 terms)
data(Penicillin, package = "lme4")
Penicillin$id <- 1

test_that("Crossed random intercept model (2 terms)",{

    ## ** fit
    eCRI2.lmer <- lmer(diameter ~ (1|plate) + (1|sample), data = Penicillin)
    eCRI2.lmm0 <- lmm(diameter ~ 1, repetition = ~1|id,
                      structure = CS(list(~1,~plate+sample), type = "ho", group = 1:2),
                      data = Penicillin, df = FALSE)
    eCRI2.lmm <- lmm(diameter ~ (1|plate) + (1|sample), data = Penicillin, df = FALSE,
                     control = list(init = "lmer"))

    ## ** iteration
    expect_equal(eCRI2.lmm0$opt$n.iter,7)
    expect_equal(eCRI2.lmm$opt$n.iter,2)

    ## ** likelihood
    expect_equal(as.double(logLik(eCRI2.lmer)), as.double(logLik(eCRI2.lmm0)), tol = 1e-6)
    expect_equal(as.double(logLik(eCRI2.lmer)), as.double(logLik(eCRI2.lmm)), tol = 1e-6)

    ## ** random effects (conditional mean)
    GS <- do.call(rbind,ranef(eCRI2.lmer))
    test <- ranef(eCRI2.lmm)
    expect_equal(as.double(GS[,1]), as.double(test$estimate), tol = 1e-6)

    ## ** random effects (conditional variance)
    GS <- as.data.frame(VarCorr(eCRI2.lmer))
    test <- ranef(eCRI2.lmm, effects = "variance", simplify = FALSE, format = "wide")
    expect_equal(as.double(test$variance[-1]), as.double(GS$vcov), tol = 1e-2)
})

## * Crossed random intercept model (3 terms)
set.seed(10)
Sigma.CRI3 <- matrix(0, nrow = 6, ncol = 6)
Sigma.CRI3[1:3,1:3] <- Sigma.CRI3[4:6,4:6] <- 0.8^2
diag(Sigma.CRI3[1:3,4:6]) <- diag(Sigma.CRI3[4:6,1:3]) <- 0.5^2
Sigma.CRI3[1,6] <- Sigma.CRI3[2,4] <- Sigma.CRI3[3,5] <- Sigma.CRI3[6,1] <- Sigma.CRI3[4,2] <- Sigma.CRI3[5,3] <- 0.3^2
diag(Sigma.CRI3) <- 1

n <- 25
dfL.CRI3 <- reshape2::melt(data.frame(sample = 1:n,
                                      mvtnorm::rmvnorm(n, mean = 1:6, sigma = Sigma.CRI3)
                                      ),
                           id.vars = "sample")
dfL.CRI3$patient <- paste(dfL.CRI3$sample,dfL.CRI3$variable %in% paste0("X",4:6)+1, sep = ".")
dfL.CRI3$day <- paste(dfL.CRI3$sample, sapply(as.character(dfL.CRI3$variable), switch, "X1" = 1, "X2" = 2, "X3" = 3, "X4" = 1, "X5" = 2, "X6" = 3),sep=".")
dfL.CRI3$batch <- paste(dfL.CRI3$sample, sapply(as.character(dfL.CRI3$variable), switch, "X1" = 1, "X2" = 2, "X3" = 3, "X4" = 2, "X5" = 3, "X6" = 1),sep=".")
dfL.CRI3 <- dfL.CRI3[order(dfL.CRI3$sample),]

## head(dfL.CRI3,10)

test_that("Crossed random intercept model (3 terms)",{

    ## ** fit
    eCRI3.lmer <- lmer(value ~ 0 + variable + (1|batch) + (1|day) + (1|patient), data = dfL.CRI3)
    eCRI3.lmm0 <- lmm(value ~ 0 + variable + (1|batch) + (1|day) + (1|patient), data = dfL.CRI3, df = FALSE)
    eCRI3.lmm <- lmm(value ~ 0 + variable + (1|batch) + (1|day) + (1|patient), data = dfL.CRI3, df = FALSE, control = list(init = "lmer"))

    ## ** iteration
    expect_equal(eCRI3.lmm0$opt$n.iter,5)
    expect_equal(eCRI3.lmm$opt$n.iter,0)

    ## ** likelihood
    expect_equal(as.double(logLik(eCRI3.lmer)), as.double(logLik(eCRI3.lmm0)), tol = 1e-6)
    expect_equal(as.double(logLik(eCRI3.lmer)), as.double(logLik(eCRI3.lmm)), tol = 1e-6)

    ## ** random effects (conditional mean)
    GS <- ranef(eCRI3.lmer)
    test <- ranef(eCRI3.lmm)
    expect_equal(as.double(test[test$variable=="patient","estimate"]), as.double(GS$patient[test[test$variable=="patient","level"],1]), tol = 1e-6)
    expect_equal(as.double(test[test$variable=="day","estimate"]), as.double(GS$day[test[test$variable=="day","level"],1]), tol = 1e-6)
    expect_equal(as.double(test[test$variable=="batch","estimate"]), as.double(GS$batch[test[test$variable=="batch","level"],1]), tol = 1e-6)

    ## ** random effects (conditional variance)
    GS <- as.data.frame(VarCorr(eCRI3.lmer))
    test <- ranef(eCRI3.lmm, effects = "variance", format = "wide", simplify = FALSE)
    expect_equal(as.double(test[-1,"variance"]), as.double(GS[,"vcov"]), tol = 1e-3)
})

## * Nested random intercept model (2 levels)
## https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified
## df.school <- read.table("http://bayes.acs.unt.edu:8083/BayesContent/class/Jon/R_SC/Module9/lmm.data.txt",
##                  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
## eNRI2.lmer <- lmer(extro ~ open + agree + social + (1 | school/class), data = dt)
## dt.red <- as.data.table(dt)[,.SD[1:10,.(id,extro=round(extro,2))],by = c("class","school")]
## setkeyv(dt.red, c("school","class"))
## dt.red$id <- 1:NROW(dt.red)

df.red <- data.frame("class" = rep(unlist(lapply(c("a","b","c","d"),rep,10)), 6), 
                     "school" = unlist(lapply(c("I","II","III","IV","V","VI"),rep,40)), 
                     "id" = 1:240, 
                     "extro" = c(41.10, 43.24, 42.15, 43.72, 40.59, 41.38, 35.43, 39.46, 41.06, 38.96, 47.24, 45.77, 44.32, 47.34, 45.06, 47.01, 46.07, 44.41, 45.30, 45.20, 47.65, 47.76, 48.59, 47.46, 48.63, 48.31, 49.28, 47.81, 49.01, 48.00, 50.97, 49.28, 49.66, 50.32, 49.33, 49.87, 50.87, 51.43, 50.95, 50.87, 52.02, 51.67, 52.62, 51.96, 51.88, 52.56, 51.54, 52.08, 52.76, 52.69, 53.69, 54.10, 54.07, 53.70, 53.03, 53.88, 53.80, 53.50, 53.42, 53.17, 55.17, 54.58, 55.06, 54.40, 54.86, 54.45, 55.27, 54.60, 55.33, 54.81, 55.63, 55.65, 55.96, 56.01, 55.58, 56.01, 56.08, 55.47, 56.07, 55.56, 56.68, 57.30, 56.84, 57.24, 56.57, 56.94, 56.81, 57.11, 57.04, 56.35, 57.75, 57.84, 57.66, 58.03, 58.22, 57.50, 57.54, 57.42, 58.11, 58.08, 59.03, 58.41, 59.04, 59.09, 58.94, 58.39, 58.30, 58.80, 58.95, 58.35, 60.15, 59.45, 59.43, 59.98, 59.40, 60.05, 60.00, 59.90, 59.93, 59.80, 60.65, 60.84, 60.66, 60.31, 60.98, 60.75, 60.60, 60.42, 60.74, 60.17, 62.29, 61.84, 61.57, 61.17, 62.12, 61.59, 61.42, 61.36, 61.67, 62.12, 62.97, 63.15, 62.82, 62.84, 63.29, 62.68, 62.45, 63.21, 62.80, 62.75, 63.69, 64.25, 64.18, 64.14, 63.72, 63.42, 64.18, 63.68, 63.48, 64.04, 64.62, 65.13, 65.40, 65.23, 64.56, 64.63, 64.81, 64.31, 64.76, 64.44, 65.92, 66.49, 66.28, 65.68, 66.03, 66.46, 65.62, 66.01, 65.61, 66.10, 66.57, 66.98, 66.55, 66.96, 66.64, 67.24, 67.23, 67.34, 66.94, 67.30, 67.94, 68.01, 68.73, 68.61, 68.89, 68.91, 68.84, 68.06, 68.13, 68.93, 69.48, 70.12, 70.17, 70.67, 70.60, 69.62, 70.16, 70.65, 70.34, 69.87, 71.79, 71.36, 72.45, 71.07, 70.87, 70.94, 71.25, 72.48, 72.04, 71.67, 74.35, 74.70, 73.14, 74.44, 73.98, 75.03, 74.21, 76.51, 75.94, 73.42, 79.74, 78.15, 78.19, 83.34, 80.11, 79.73, 80.64, 77.07, 78.02, 82.25))

test_that("Nested random intercept model (2 levels)",{

    ## ** fit
    eNRI2.lmer <- lmer(extro ~ (1|school/class), data = df.red)
    eNRI2.lmm0 <- lmm(extro ~ (1|school/class), data = df.red, df = FALSE)
    eNRI2.lmm <- lmm(extro ~ (1|school/class), data = df.red, df = FALSE, control = list(init = "lmer"))

        ## ** iteration
    expect_equal(eNRI2.lmm0$opt$n.iter,9)
    expect_equal(eNRI2.lmm$opt$n.iter,3)

    ## ** likelihood
    expect_equal(as.double(logLik(eNRI2.lmer)), as.double(logLik(eNRI2.lmm0)), tol = 1e-6)
    expect_equal(as.double(logLik(eNRI2.lmer)), as.double(logLik(eNRI2.lmm)), tol = 1e-6)

    ## ** random effects (conditional mean)
    GS <- ranef(eNRI2.lmer)
    test <- ranef(eNRI2.lmm)
    test$variable2 <- sapply(test$variable,paste, collapse=":")
    expect_equal(as.double(test[test$variable2=="school","estimate"]), as.double(GS$school[,1]), tol = 1e-4)
    expect_equal(as.double(test[test$variable2=="school:class","estimate"]), as.double(GS$class[,1]), tol = 1e-4)

    ## ** random effects (conditional variance)
    GS <- as.data.frame(VarCorr(eNRI2.lmer))
    test <- ranef(eNRI2.lmm, effects = "variance", format = "wide", simplify = FALSE)
    expect_equal(as.double(test[-1,"variance"]), as.double(GS[c(2,1,3),"vcov"]), tol = 1e-3)
})


## * Nested random intercept model (3 levels)
Sigma.NRI3 <- matrix(0, nrow = 8, ncol = 8)
Sigma.NRI3[5:8,1:4] <- Sigma.NRI3[1:4,5:8] <- 0.25^2
Sigma.NRI3[1:2,3:4] <- Sigma.NRI3[3:4,1:2] <- Sigma.NRI3[5:6,7:8] <- Sigma.NRI3[7:8,5:6] <- 0.5^2
Sigma.NRI3[1:2,1:2] <- Sigma.NRI3[3:4,3:4] <- Sigma.NRI3[5:6,5:6] <- Sigma.NRI3[7:8,7:8] <- 0.8^2
diag(Sigma.NRI3) <- 1

## c(sqrt(0.25^2), sqrt(0.5^2-0.25^2), sqrt(0.8^2-0.5^2-0.25^2))
n <- 1000
dfL.NRI3 <- reshape2::melt(data.frame(patient = 1:n,
                                      mvtnorm::rmvnorm(n, mean = 1:8, sigma = Sigma.NRI3)
                                      ),
                           id.vars = "patient")
dfL.NRI3$day <- dfL.NRI3$variable %in% paste0("X",5:8)+1
dfL.NRI3$session <- paste(dfL.NRI3$day,dfL.NRI3$variable %in% c(paste0("X",c(3:4,7:8)))+1,sep=".")
dfL.NRI3 <- dfL.NRI3[order(dfL.NRI3$patient),]

## head(dfL.NRI3,10)

test_that("Nested random intercept model (2 levels)",{

    ## ** fit
    eNRI3.lmer <- lmer(value ~ session + (1|patient/day/session), data = dfL.NRI3)
    eNRI3.lmm0 <- lmm(value ~ session + (1|patient/day/session), data = dfL.NRI3, df = FALSE)
    eNRI3.lmm <- lmm(value ~ session + (1|patient/day/session), data = dfL.NRI3, df = FALSE, control = list(init = "lmer"))

        ## ** iteration
    expect_equal(eNRI3.lmm0$opt$n.iter,4)
    expect_equal(eNRI3.lmm$opt$n.iter,2)

    ## ** likelihood
    expect_equal(as.double(logLik(eNRI3.lmer)), as.double(logLik(eNRI3.lmm0)), tol = 1e-6)
    expect_equal(as.double(logLik(eNRI3.lmer)), as.double(logLik(eNRI3.lmm)), tol = 1e-6)

    ## ** random effects (conditional mean)
    ## slow!!
    ## GS <- ranef(eNRI3.lmer)
    ## test <- ranef(eNRI3.lmm)
    ## test$variable2 <- sapply(test$variable,paste, collapse=":")
    ## expect_equal(as.double(test[test$variable2=="patient","estimate"]), as.double(GS$patient[,1]), tol = 1e-4)
    ## expect_equal(as.double(test[test$variable2=="patient:day","estimate"]), as.double(GS$day[,1]), tol = 1e-4)
    ## expect_equal(as.double(test[test$variable2=="patient:day:session","estimate"]), as.double(GS$session[,1]), tol = 1e-4)

    ## ** random effects (conditional variance)
    GS <- as.data.frame(VarCorr(eNRI3.lmer))
    test <- ranef(eNRI3.lmm, effects = "variance", format = "wide", simplify = FALSE)
    expect_equal(as.double(test[-1,"variance"]), as.double(GS[c(3,2,1,4),"vcov"]), tol = 1e-3)
})


##----------------------------------------------------------------------
### test-auto-mixed-model.R ends here
