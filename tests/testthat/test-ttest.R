### test-ttest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:20) 
## Version: 
## Last-Updated: Jun  1 2021 (12:16) 
##           By: Brice Ozenne
##     Update #: 10
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
    library(lava)
    library(multcomp)
    library(nlme)

    library(LMMstar)
}

context("Check lmm on examples of t-tests")
LMMstar.options(method.numDeriv = "Richardson")

## * single t-test
## ** simulate data
n <- 15
p <- 3
X.name <- paste0("X",1:p)
link.lvm <- paste0("Y~",X.name)
formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))

m <- lvm(formula.lvm)
distribution(m,~Id) <- Sequence.lvm(a = 1, b = n)
distribution(m,~Gender.num) <- binomial.lvm()
distribution(m,~X3) <- binomial.lvm(size=2)
transform(m,Gene~X3) <- function(x){factor(x,levels=0:2,labels=c("LL","LA","AA"))}
transform(m,Gender~Gender.num) <- function(x){factor(x,levels=0:1,labels=c("M","F"))}
transform(m,id~Gender.num) <- function(x){paste0("id",1:NROW(x))}
latent(m) <- ~Gender.num+X3
set.seed(10)
d <- lava::sim(m,n, latent = FALSE)
d$time <- "t1"

## ** test
test_that("single t-test",{
    e.lmm <- lmm(Y ~ Gender-1, structure = UN(~Gender), data = d, trace = 0)
    e2.lmm <- lmm(Y ~ Gender-1, repetition = ~Gender|id, structure = "UN", data = d, trace = 0, df = FALSE)
    ## e3.lmm <- lmm(Y ~ Gender, repetition = ~Gender|id, structure = "UN", data = d, trace = 0)
    e.gls <- gls(Y ~ Gender, weights = varIdent(form=~1|Gender), data = d)
    e.ttest <- t.test(Y~Gender, data = d)
    
    expect_equal(logLik(e.lmm),logLik(e2.lmm))
    expect_equal(unname(e.ttest$estimate), unname(coef(e.lmm, effects = "mean")))
    expect_equal(e.ttest$p.value, anova(e.lmm, effects = c("GenderM-GenderF=0"))$all$p.value, tol = 1e-5)
    ## summary(e.gls)$tTable
})
## confint(anova(e.lmm))

## * multiple t-test
## ** simulate data
m <- lvm(c(Y1,Y2,Y3,Y4) ~ age + gender)
categorical(m, labels = c("male","female")) <- ~gender
transform(m, id~gender) <- function(x){1:NROW(x)}
distribution(m, ~age) <- gaussian.lvm(mean = 50, sd = 10)

set.seed(10)
dW <- lava::sim(m, 30)

## move to the long format
name.varying <- paste0("Y",1:4)
dL <- reshape(dW, direction  = "long",
              idvar = c("id","age","gender"),
              varying = name.varying,
              v.names = "Y",
              timevar = "visit")
rownames(dL) <- NULL
dL$visit <- factor(dL$visit,
                   levels = 1:length(name.varying),
                   labels = name.varying)
 

## ** test
test_that("multiple t-test",{
    e.lmm <- lmm(Y ~ visit + gender:visit - 1, structure = UN(gender~visit|id), data = dL, trace = 0)
    e2.lmm <- lmm(Y ~ visit + gender:visit - 1, repetition = gender~visit|id, structure = "UN", data = dL, trace = 0, df = FALSE)
    e.gls <- gls(Y ~ visit + gender:visit, correlation = corSymm(form=~1|id), weights = varIdent(form=~1|visit), data = dL)
    ls.ttest <- list(t.test(Y1~gender, data = dW),
                     t.test(Y2~gender, data = dW),
                     t.test(Y3~gender, data = dW),
                     t.test(Y4~gender, data = dW))
    e.lh <- anova(e.lmm, effects = c("gender:Y1"="visitY1:male - visitY1:female = 0",
                                     "gender:Y2"="visitY2:male - visitY2:female = 0",
                                     "gender:Y3"="visitY3:male - visitY3:female = 0",
                                     "gender:Y4"="visitY4:male - visitY4:female = 0"), ci = TRUE)

    expect_equal(confint(e.lh)$all[[1]]$se, sapply(ls.ttest,"[[","stderr"), tol = 1e-5)
    expect_equal(as.double(unlist(confint(e.lh, method = "none")$all[[1]][,c("lower","upper")])),
                 as.double(do.call(rbind,lapply(ls.ttest,"[[","conf.int"))), tol = 1e-3)
})

##----------------------------------------------------------------------
### test-ttest.R ends here
