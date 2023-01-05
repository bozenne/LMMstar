## chunk 2
library(LMMstar)

## chunk 3
utils::packageVersion("LMMstar")

## * Illustrative dataset

## chunk 4
data(gastricbypassL, package = "LMMstar")
head(gastricbypassL)

## chunk 5
gastricbypassL$time <- factor(gastricbypassL$time,
                              levels = c("3monthsBefore", "1weekBefore",
                                         "1weekAfter", "3monthsAfter" ),
                              labels = c("B3m","B1w","A1w","A3m"))
gastricbypassL$visit <- as.numeric(gastricbypassL$time) ## convert to numeric
gastricbypassL$baseline <- gastricbypassL$visit<=2

## chunk 6
gastricbypassL$glucagon <- as.double(scale(gastricbypassL$glucagonAUC))+5

## chunk 7
gastricbypassL$group <- as.numeric(gastricbypassL$id)%%2

## chunk 8
data(gastricbypassW, package = "LMMstar")
head(gastricbypassW)

## chunk 9
gastricbypassW$group <- as.numeric(gastricbypassW$id)%%2

## chunk 10
dfL <- gastricbypassL[!is.na(gastricbypassL$glucagonAUC),]

## * Descriptive statistics
## ** Summary statistics

## chunk 11
sss <- summarize(weight+glucagon ~ time, data = gastricbypassL, na.rm = TRUE)
print(sss, digits = 3)

## chunk 12
sss <- summarize(weight ~ time|id, data = gastricbypassL, na.rm = TRUE)
print(sss, digits = 3)

## chunk 13
partialCor(weight1 + weight4 ~ 1, data = gastricbypassW)

## chunk 14
partialCor(list(weight1 ~ glucagonAUC1, weight4 ~ glucagonAUC4),
           data = gastricbypassW)

## chunk 15
partialCor(weight + glucagonAUC ~ 1, by = "group", data = gastricbypassL)

## chunk 16
partialCor(weight + glucagonAUC ~ 1, by = "group", effects = "Dunnett",
           data = gastricbypassL)

## ** Missing data patterns

## chunk 17
mp <- summarizeNA(gastricbypassL)
mp

## chunk 18
plot(mp)

## * Linear mixed model
## ** Classical covariance patterns

## chunk 20
eId.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id, 
               structure = "ID", data = dfL)
eId.lmm
cat(" covariance structure: \n");sigma(eId.lmm)

## chunk 21
eInd.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id, 
               structure = "IND", data = dfL)
eInd.lmm
cat(" covariance structure: \n");sigma(eInd.lmm)

## chunk 22
eCS.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id,
               structure = "CS", data = dfL)
eCS.lmm
cat(" covariance structure: \n");sigma(eCS.lmm)

## chunk 23
eTOE.lmm <- lmm(weight ~ time*group, repetition = ~time|id,
                structure = "TOEPLITZ", data = dfL)
eTOE.lmm
cat(" correlation structure: \n");cov2cor(sigma(eTOE.lmm))

## *unstructured* covariance matrix:

## chunk 24
eUN.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id,
               structure = "UN", data = dfL)
eUN.lmm
cat(" covariance structure: \n");sigma(eUN.lmm)

## chunk 25
eSCS.lmm <- lmm(weight ~ time*group,
                repetition = ~time|id, structure = CS(group~1),
                data = dfL)
eSCS.lmm

## chunk 26
eSUN.lmm <- lmm(weight ~ time*group + glucagon,
                repetition = ~time|id, structure = UN(~group),
                data = dfL)
eSUN.lmm

## chunk 27
sigma(eSCS.lmm)

## chunk 28
sigma(eSUN.lmm)

## chunk 29
eBCS.lmm <- lmm(weight ~ time*group,repetition = ~time|id,
                structure = CS(~baseline, heterogeneous = FALSE), data = dfL)
eBCS.lmm
cat(" covariance structure: \n");sigma(eBCS.lmm)

## chunk 31
eBUN.lmm <- lmm(weight ~ time*group, repetition = ~time|id,
                structure = CS(~baseline, heterogeneous = TRUE), data = dfL)
eBUN.lmm
cat(" covariance structure: \n");sigma(eBUN.lmm)

## ** User-specific covariance patterns

## chunk 32
rho.2block <- function(p,time,...){
  n.time <- length(time)
  rho <- matrix(1, nrow = n.time, ncol = n.time)
  rho[1,2] <- rho[2,1] <- rho[4,5] <- rho[5,4] <- p["rho1"]
  rho[1,3] <- rho[3,1] <- rho[4,6] <- rho[6,4] <- p["rho2"]
  rho[2,3] <- rho[3,2] <- rho[5,6] <- rho[6,5] <- p["rho3"]
  rho[4:6,1:3] <- rho[1:3,4:6] <- p["rho4"]
  return(rho)
}
Rho <- rho.2block(p = c(rho1=0.25,rho2=0.5,rho3=0.4,rho4=0.1),
                  time = 1:6)
Rho

## chunk 33
set.seed(11)
n <- 1000
Y <- rmvnorm(n, mean = rep(0,6), sigma = Rho)
dfL2 <- reshape2::melt(cbind(id = 1:n, as.data.frame(Y)), id.vars = "id")
dfL2$time  <- dfL2$variable
dfL2 <- dfL2[order(dfL2$id),]
dfL2[1:8,]

## chunk 34
myStruct <- CUSTOM(~variable,
                   FCT.sigma = function(p,time,X){rep(p,length(time))}, ## function f
                   init.sigma = c("sigma"=1),
                   FCT.rho = rho.2block, ## function g
                   init.rho = c("rho1"=0.25,"rho2"=0.25,"rho3"=0.25,"rho4"=0.25))

## chunk 35
e.lmmCUSTOM <- lmm(value~time,
                   repetition=~time|id,
                   structure = myStruct,
                   data=dfL2,
                   df = FALSE) ## df = FALSE to save computation time
logLik(e.lmmCUSTOM)

## chunk 36
cov2cor(sigma(e.lmmCUSTOM))

## chunk 37
myCS <- CUSTOM(~1,
       FCT.sigma = function(p,time,X){rep(p,length(time))},
       init.sigma = c("sigma"=1),
       FCT.rho = function(p,time,X){matrix(p,length(time),length(time))+diag(1-p,length(time),length(time))},
       init.rho = c("rho"=0.5))

## chunk 38
logLik(lmm(value~time,
           repetition = ~time|id,
           structure = myCS, 
           data = dfL2, df = FALSE
           ))

## chunk 39
logLik(lmm(value~time,
           repetition = ~time|id,
           structure = "CS", 
           data = dfL2, df = FALSE))

## ** Model output

## chunk 40
summary(eUN.lmm)

## chunk 41
summary(eUN.lmm, hide.mean = TRUE)

## chunk 42
oo <- capture.output(summary(eUN.lmm, hide.fit = TRUE, hide.data = TRUE, hide.cor = TRUE, hide.var = TRUE, hide.sd = TRUE))
cat(sapply(oo[-(1:2)],paste0,"\n"))

## ** Extract estimated coefficients

## chunk 43
coef(eUN.lmm)

## chunk 44
coef(eUN.lmm, effects = "variance")

## chunk 45
coef(eUN.lmm, effects = "variance", transform.k = "sd")

## chunk 46
dummy.coef(eUN.lmm)

## ** Extract estimated coefficient and associated uncertainty

## chunk 47
model.tables(eUN.lmm)

## chunk 48
model.tables(eUN.lmm, effect = "all")

## chunk 49
model.tables(eUN.lmm, columns = c("estimate","p.value"))

## chunk 50
model.tables(eUN.lmm, columns = add("statistic"))

## ** Extract estimated residual variance-covariance structure

## chunk 51
Sigma <- sigma(eUN.lmm)
Sigma

## chunk 52
cov2cor(Sigma)

## chunk 53
sigma(eUN.lmm, cluster = 5)

## chunk 54
newdata <- data.frame(id = "X", time = c("B3m","B1w","A1w","A3m"))
sigma(eUN.lmm, cluster = newdata)

## ** Random effects

## chunk 55
head(coef(eCS.lmm, effects = "ranef"))

## chunk 56
head(coef(eBCS.lmm, effects = "ranef"))

## ** Sum of squares

## chunk 58
sigma2 <- coef(eCS.lmm, effect = "variance")^2
tau <- coef(eCS.lmm, effect = "correlation")*sigma2
omega <- unname(sigma2 - tau)

## chunk 59
df.res <- df.residual(eCS.lmm)
SSE <- df.res * omega
c(df.res = df.res, SSE = SSE)

## chunk 60
eBeta.lmm <- coef(eCS.lmm)
eVcov.lmm <- vcov(eCS.lmm, type.information = "expected")

## chunk 61
attr(model.matrix(eCS.lmm),"assign")

## chunk 62
SSRstar.time <- eBeta.lmm[2:4] %*% solve(eVcov.lmm[2:4,2:4]) %*% eBeta.lmm[2:4] 
SSRstar.glucagon <- eBeta.lmm[5] %*% solve(eVcov.lmm[5,5]) %*% eBeta.lmm[5] 

## chunk 63
SSR.time <- as.double(SSRstar.time * omega)
SSR.glucagon <- as.double(SSRstar.glucagon * omega)
c(time = SSR.time, glucagon = SSR.glucagon)

## ** Proportion of explained variance and partial correlation

## chunk 65
c(SSR.time/ (SSR.time + SSE),
  SSR.glucagon/ (SSR.glucagon + SSE))

## chunk 66
eCS.R2 <- partialCor(eCS.lmm, R2 = TRUE)
summary(eCS.R2)

## chunk 67
aCS.aov <- anova(eCS.lmm)$multivariate
setNames(aCS.aov$statistic/(aCS.aov$statistic+aCS.aov$df.denom), aCS.aov$test)

## ** Model diagnostic

## chunk 69
plot(eUN.lmm, type = "scatterplot")

## chunk 70
plot(eUN.lmm, type = "scatterplot2")

## chunk 71
plot(eUN.lmm, type = "correlation", type.residual = "response")
plot(eUN.lmm, type = "correlation", type.residual = "normalized")

## chunk 73
plot(eUN.lmm, type = "qqplot", engine.qqplot = "qqtest")
## Note: the qqtest package to be installed to use the argument engine.plot = "qqtest" 

## chunk 74
eUN.diagW <- residuals(eUN.lmm, type = "normalized", format = "wide")
colnames(eUN.diagW) <- gsub("normalized.","",colnames(eUN.diagW))
head(eUN.diagW)

## chunk 75
eUN.diagL <- residuals(eUN.lmm, type = "normalized", format = "long")
head(eUN.diagL)

## ** Model fit

## chunk 76
library(ggplot2) ## left panel
plot(eUN.lmm, type = "fit", color = "id", ci.alpha = NA, size.text = 20)

## chunk 77
library(emmeans) ## right panel
emmip(eUN.lmm, ~time) + theme(text = element_text(size=20))

## chunk 78
plot(eUN.lmm, type = "fit", at = data.frame(glucagon = 10), color = "glucagon")
## result not shown

## chunk 79
gg <- plot(eUN.lmm, type = "fit", obs.alpha = 0.2, ci = FALSE,plot = FALSE)$plot
gg <- gg + facet_wrap(~id, labeller = label_both)
gg <- gg + theme(axis.text.x=element_text(angle = 90, hjust = 0))
gg

## ** Partial residuals

## chunk 81
gg1 <- plot(eUN.lmm, type = "partial", var = "glucagon", plot = FALSE)$plot
gg2 <- plot(eUN.lmm, type = "partial", var = c("(Intercept)","glucagon"), plot = FALSE)$plot
ggarrange(gg1,gg2)

## chunk 82
df.pres <- residuals(eUN.lmm, type = "partial", var = "glucagon", keep.data = TRUE)
head(df.pres)

## chunk 83
m.pres <- dfL$weight - model.matrix(~time,dfL) %*% coef(eUN.lmm)[1:4]
range(df.pres$r.partial - m.pres, na.rm = TRUE)

## chunk 84
eIID.lm <- lm(glucagon ~ time + weight, data = dfL)
pRes.lm <- residuals(eIID.lm, type = "partial")[,"weight"]

## chunk 85
eIID.lmm <- lmm(glucagon ~ time + weight, data = dfL)
pRes.lmm <- residuals(eIID.lmm, type = "partial-center", var = "weight")
range(pRes.lm-na.omit(pRes.lmm))

## ** Statistical inference (linear)

## chunk 87
anova(eUN.lmm)

## chunk 88
anova(eUN.lmm, effects = c("timeA1w-timeB1w=0"))

## chunk 89
e.anova <- anova(eUN.lmm, effects = c("timeA1w-timeB1w=0","timeA3m-timeB1w=0"))
summary(e.anova)

## chunk 90
library(multcomp)
summary(anova(eUN.lmm, effects = mcp(time = "Tukey")))

## chunk 91
try(
  anova(eUN.lmm,
        effects = c("log(k).B1w=0","log(k).A1w=0","log(k).A3m=0"))
)

## chunk 92
name.coef <- rownames(confint(eUN.lmm, effects = "all"))
name.varcoef <- grep("^k",name.coef, value = TRUE)
C <- matrix(0, nrow = 3, ncol = length(name.coef), dimnames = list(name.varcoef, name.coef))
diag(C[name.varcoef,name.varcoef]) <- 1
C[,1:9]

## chunk 93
anova(eUN.lmm, effects = C)

## chunk 94
Manova <- rbind(anova(eInd.lmm, effects = "glucagon = 0"),
                anova(eCS.lmm, effects = "glucagon = 0"),
                anova(eUN.lmm, effects = "glucagon = 0"),
                name = c("Ind","CS","UN"))
summary(Manova) 

## ** Statistical inference (non-linear)

## chunk 95
gastricbypassW <- reshape(dfL[,c("id","time","weight","group")],
                          direction = "wide",
                          timevar = "time", idvar = c("id","group"))
e.ANCOVA <- lm(weight.A1w ~ weight.B1w + group, data = gastricbypassW)
summary(e.ANCOVA)$coef

## chunk 96
dfL23 <- dfL[dfL$visit %in% 2:3,]
dfL23$time <- droplevels(dfL23$time)
e.lmmANCOVA <- lmm(weight ~ time+time:group, repetition = ~time|id,
                   data = dfL23)

## chunk 97
lava::estimate(e.lmmANCOVA, f = function(p){
  c(Y1 = as.double(p["rho(B1w,A1w)"]*p["k.A1w"]),
    X1 = as.double(p["timeA1w:group"]-p["rho(B1w,A1w)"]*p["k.A1w"]*p["timeB1w:group"]))
})

## ** Baseline adjustment

## chunk 98
gastricbypassL$treat <- baselineAdjustment(gastricbypassL, variable = "group",
                                repetition = ~time|id, constrain = c("B3m","B1w"),
                                new.level = "none")
table(treat = gastricbypassL$treat, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 99
gastricbypassL$treat2 <- baselineAdjustment(gastricbypassL, variable = "group",
                                            repetition = ~time|id, constrain = c("B3m","B1w"))
table(treat = gastricbypassL$treat2, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 100
gastricbypassL$timeXtreat <- baselineAdjustment(gastricbypassL, variable = "group",
                                                repetition = ~time|id, constrain = c("B3m","B1w"),
                                                collapse.time = ".")

table(treat = gastricbypassL$timeXtreat, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 101
eC.lmm <- lmm(weight ~ timeXtreat, data = gastricbypassL,
              repetition = ~time|id, structure = "UN")
coef(eC.lmm) ## change from baseline

## chunk 102
eC2.lmm <- lmm(weight ~ 0 + timeXtreat, data = gastricbypassL,
              repetition = ~time|id, structure = "UN")
coef(eC2.lmm) ## absolute value

## chunk 103
colnames(model.matrix(weight ~ treat*time, data = gastricbypassL))

## chunk 104
eC3.lmm <- lmm(weight ~ treat2*time, data = gastricbypassL,
               repetition = ~time|id, structure = "UN")

## chunk 105
model.tables(eC3.lmm)

## chunk 106
autoplot(eC3.lmm, color = "group", ci = FALSE, size.text = 20, obs.alpha = 0.1) 

## ** Marginal means

## chunk 107
gastricbypassL$group2 <- as.numeric(gastricbypassL$id) %% 3 == 0
e.group <- lmm(glucagon ~ time*group2, data = gastricbypassL,
               repetition = ~time|id, structure = "UN")

## chunk 108
emmeans(e.group, specs=~time)

## chunk 109
df.pred <- cbind(gastricbypassL, predict(e.group, newdata = gastricbypassL))
summarize(formula = estimate~time, data = df.pred)

## chunk 110
table(group = dfL$group2, time = dfL$time)

## chunk 111
mu.group1 <-  as.double(coef(e.group)["(Intercept)"])
mu.group2 <-  as.double(coef(e.group)["(Intercept)"] + coef(e.group)["group2TRUE"])
p.group1 <- 14/20          ; p.group2 <- 6/20
c(emmeans = (mu.group1+mu.group2)/2, predict = mu.group1 * p.group1 + mu.group2 * p.group2)

## chunk 112
emmeans.group <- emmeans(e.group, specs = ~group2|time)
emmeans.group

## chunk 113
epairs.group <- pairs(emmeans.group, reverse = TRUE)
epairs.group

## chunk 114
summary(epairs.group, by = NULL, adjust = "mvt", infer = TRUE)

## chunk 115
summary(pairs(emmeans(eC3.lmm , specs = ~treat2|time), reverse = TRUE), by = NULL)

## ** Predictions

## chunk 116
news <- dfL[dfL$id==1,]
news$glucagon <- 0
predict(eUN.lmm, newdata = news)

## chunk 117
X.12 <- model.matrix(formula(eUN.lmm), news)
X.12

## chunk 118
X.12 %*% coef(eUN.lmm)

## chunk 119
newd <- rbind(
  data.frame(id = 1, time = "B3m", weight = coef(eUN.lmm)["(Intercept)"], glucagon = 0),
  data.frame(id = 1, time = "B1w", weight = NA, glucagon = 0),
  data.frame(id = 2, time = "B3m", weight = 100, glucagon = 0),
  data.frame(id = 2, time = "B1w", weight = NA, glucagon = 0)
)
predict(eUN.lmm, newdata = newd, type = "dynamic", keep.newdata = TRUE)

## chunk 120
mu1 <- coef(eUN.lmm)[1]
mu2 <- sum(coef(eUN.lmm)[1:2])
Omega_11 <- sigma(eUN.lmm)["B3m","B3m"]
Omega_21 <- sigma(eUN.lmm)["B1w","B3m"]
as.double(mu2 + Omega_21 * (100 - mu1) / Omega_11)

## * Missing values and imputation
## ** Full information approach

## chunk 121
sss <- summarize(glucagon ~ time, data = gastricbypassL, na.rm = TRUE)
cbind(sss[,1:4], pc = paste0(100 * sss$missing / (sss$missing + sss$observed), "%"))

## chunk 122
vec.pattern <- tapply(as.numeric(is.na(gastricbypassL$glucagon)),
                      INDEX = gastricbypassL$id,
                      FUN = paste, collapse=".")
table(vec.pattern)

## chunk 123
eUN.lmmNA <- lmm(glucagon ~ time,
                 repetition = ~time|id, structure = "UN",
                 data = gastricbypassL)
summary(eUN.lmmNA)

## chunk 124
summary(eUN.lmmNA, hide.mean = TRUE, hide.sd = TRUE, hide.var = TRUE, hide.cor = TRUE, hide.fit = TRUE)

## ** Imputation 

## chunk 125
fitted(eUN.lmmNA, impute = TRUE)

## chunk 126
eData <- fitted(eUN.lmmNA, impute = TRUE, keep.newdata = TRUE)
eData$treat <- eData$treat2 <- eData$timeXtreat <- NULL
eData[eData$id %in% eData[eData$imputed,"id"],]

## chunk 127
ggplot(eData, aes(x=time,y=glucagon, group=id)) + geom_line() + geom_point(aes(color=imputed))

## chunk 129
set.seed(10)
fitted(eUN.lmmNA, impute = TRUE, se = "total")
fitted(eUN.lmmNA, impute = TRUE, se = "total")
fitted(eUN.lmmNA, impute = TRUE, se = "total")

## ** Multiple imputation

## chunk 130
data(gastricbypassW, package = "LMMstar")
colSums(is.na(gastricbypassW))

## chunk 131
library(mice)
set.seed(10)
gastricbypassW.mice <- mice(gastricbypassW, m = 5, printFlag = FALSE)
gastricbypassW.NNA <- complete(gastricbypassW.mice, action = "long")
table(gastricbypassW.NNA$.imp)

## chunk 132
e.mlmm <- mlmm(glucagonAUC3~glucagonAUC2+weight2, data=gastricbypassW.NNA,
               by = ".imp", effects = "weight2=0", trace = FALSE)
model.tables(e.mlmm)

## chunk 133
model.tables(e.mlmm, method = "pool.rubin")

## chunk 134
e.mice <- with(data=gastricbypassW.mice,exp=lm(glucagonAUC3~glucagonAUC2+weight2))
summary(pool(e.mice))

## chunk 135
plot(e.mlmm, method = c("pool.rubin","none"))

## * Data generation

## chunk 137
set.seed(10) ## ensure reproductibility
n.obs <- 100
n.times <- 4
mu <- rep(0,4)
gamma <- matrix(0, nrow = n.times, ncol = 10) ## add interaction
gamma[,6] <- c(0,1,1.5,1.5)
dW <- sampleRem(n.obs, n.times = n.times, mu = mu, gamma = gamma, format = "wide")
head(round(dW,3))

## chunk 138
set.seed(10) ## ensure reproductibility
dL <- sampleRem(n.obs, n.times = n.times, mu = mu, gamma = gamma, format = "long")
head(dL)

## * Modifying default options

## chunk 139
LMMstar.options("type.information")

## chunk 140
LMMstar.options(type.information = "expected")

## chunk 141
LMMstar.options(reinitialise = TRUE)

## * R session

## chunk 142
sessionInfo()

## * References
## * Likelihood in a linear mixed model
## ** Log-likelihood
## ** Score
## ** Hessian
## ** Degrees of freedom
## * Likelihood ratio test with the REML criterion

## chunk 143
## data(dfL, package = "LMMstar")
dfTest <- dfL
dfTest$glucagon2 <- dfTest$glucagon*2

## chunk 144
eML.lmmUN <- lmm(weight ~ time+glucagon, data = dfTest, repetition = ~time|id, method = "ML")
eML.lmmUN2 <- lmm(weight ~ time+glucagon2, data = dfTest, repetition = ~time|id, method = "ML")

## chunk 145
logLik(eML.lmmUN)
logLik(eML.lmmUN2)

## chunk 146
eREML.lmmUN <- lmm(weight ~ time + glucagon, data = dfTest, repetition = ~time|id, method = "REML")
eREML.lmmUN2 <- lmm(weight ~ time + glucagon2, data = dfTest, repetition = ~time|id, method = "REML")

## chunk 147
logLik(eREML.lmmUN)-logLik(eREML.lmmUN2)
log(2)

## chunk 148
set.seed(5) 
dfTest$ff <- rbinom(NROW(dfTest), size = 1, prob = 0.5)
logLik(lmm(weight ~ time+glucagon, data = dfTest, repetition = ~time|id, method = "REML"))
logLik(lmm(weight ~ time+glucagon*ff, data = dfTest, repetition = ~time|id, method = "REML"))

## chunk 149
logLik(lmm(weight ~ time + glucagon, data = dfTest, repetition = ~time|id, method = "ML"))
logLik(lmm(weight ~ time + glucagon*ff, data = dfTest, repetition = ~time|id, method = "ML"))

## * Sum of squares in a linear mixed model
## *Illustration for a univariate linear model:*

## chunk 150
df.aov <- dfL[!is.na(dfL$glucagon),]

## chunk 151
e.lm <- lm(weight ~ time + glucagon, data = df.aov)
car::Anova(e.lm, type = "II")

## chunk 153
e.lmm <- lmm(weight ~ time + glucagon, data = df.aov)

## chunk 154
SSEstar <- crossprod(residuals(e.lmm, type = "normalized"))
c(SSEstar = SSEstar, SSE = SSEstar * sigma(e.lmm))

## chunk 155
df.residual(e.lmm)

## chunk 156
eBeta.lmm <- coef(e.lmm)
eVcov.lmm <- vcov(e.lmm, type.information = "expected")

SSRstar.glucagon <- eBeta.lmm[5] %*% solve(eVcov.lmm[5,5]) %*% eBeta.lmm[5] 
SSRstar.time <- eBeta.lmm[2:4] %*% solve(eVcov.lmm[2:4,2:4]) %*% eBeta.lmm[2:4] 
c(SSR.glucagon = SSRstar.glucagon * sigma(e.lmm),
  SSR.time = SSRstar.time * sigma(e.lmm),
  F.glucagon = SSRstar.glucagon,
  F.time = SSRstar.time/3)

## chunk 157
R2.glucagon <- SSRstar.glucagon/(SSRstar.glucagon+SSEstar)
R2.glucagon

## chunk 158
sign(coef(e.lmm)["glucagon"])*sqrt(R2.glucagon)

## chunk 159
summary(partialCor(e.lmm, R2 = TRUE))

## * Equivalent with other R packages
## ** nlme package

## chunk 160
library(nlme)

## chunk 161
eCS.gls <- gls(weight ~ time + glucagon, correlation = corCompSymm(form=~time|id),
               data = dfL, na.action = na.omit)
eCS.lme <- lme(weight ~ time + glucagon, random = ~1|id,
               data = dfL, na.action = na.omit)
logLik(eCS.lme)
logLik(eCS.gls)
logLik(eCS.lmm)

## chunk 162
range(coef(eCS.lmm, effects = "ranef")-ranef(eCS.lme))

## chunk 163
eUN.gls <- gls(weight ~ time + glucagon,
               correlation = corSymm(form=~as.numeric(time)|id),
               weights = varIdent(form=~1|time),
               data = dfL, na.action = na.omit)
logLik(eUN.gls)
logLik(eUN.lmm)

## ** lme4 package

## chunk 164
library(lme4)
library(lmerTest)

## chunk 165
eCS.lmer <- lmer(weight ~ time + glucagon + (1|id),
                 data = dfL)
logLik(eCS.lmer)
logLik(eCS.lmm)

## chunk 166
range(coef(eCS.lmm, effects = "ranef")-ranef(eCS.lmer)$id)

## chunk 167
eBCS.lmer <- lmer(weight ~ time*group + (1|id/baseline),
                  data = dfL)
logLik(eBCS.lmer)
logLik(eBCS.lmm)

## chunk 168
eRanefBCS.lmm <- coef(eBCS.lmm, effects = "ranef")
eRanefBCS.lmer <- ranef(eBCS.lmer)
## id
range(eRanefBCS.lmm[,"id"]-eRanefBCS.lmer$id)
## baseline
range(c(eRanefBCS.lmm[,"baseline1"],eRanefBCS.lmm[,"baseline2"])-ranef(eBCS.lmer)$`baseline:id`)

## chunk 169
eUN.lmer <- lmer(weight ~ time + glucagon + (0 + time|id),
                 data = dfL, control = lmerControl(check.nobs.vs.nRE = "ignore"))
logLik(eUN.lmer)
logLik(eUN.lmm)

## chunk 170
anova(eUN.lmm)

## chunk 171
anova(eUN.lmer)

## chunk 172
eUN2.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id,
                structure = "UN", data = dfL, type.information = "expected")
suppressWarnings(anova(eUN2.lmm))

## ** mmrm package

## chunk 173
library(mmrm)
e.mmrm <- mmrm(
  formula = FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID),
  data = fev_data
)

## chunk 174
e.lmm <- lmm(
  formula = FEV1 ~ RACE + SEX + ARMCD * AVISIT,
  repetition = ~ AVISIT | USUBJID, structure = "UN",
  data = fev_data, type.information = "expected"
)

## chunk 175
logLik(e.mmrm) - logLik(e.lmm)
range(coef(e.mmrm) - coef(e.lmm))
range(vcov(e.mmrm) - vcov(e.lmm))

## ** effectsize package (\(R^2\) or \(\eta^2\))

## chunk 176
library(effectsize)
eta_squared(eCS.lmer)
cat("\n")

## chunk 177
eCS.Wald <- anova(eCS.lmm)$multivariate
eCS.Wald$df.num*eCS.Wald$statistic/(eCS.Wald$df.num*eCS.Wald$statistic+eCS.Wald$df.denom)

## chunk 178
eUN.Wald <- anova(eUN.lmm)$multivariate
eUN.Wald$df.num*eUN.Wald$statistic/(eUN.Wald$df.num*eUN.Wald$statistic+eUN.Wald$df.denom)

## chunk 179
eta_squared(eUN.lmer)
cat("\n")

