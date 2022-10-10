## chunk 2
library(LMMstar)

## chunk 3
data(gastricbypassL, package = "LMMstar")
head(gastricbypassL)

## chunk 4
gastricbypassL$time <- factor(gastricbypassL$time,
                              levels = c("3monthsBefore", "1weekBefore",
                                         "1weekAfter", "3monthsAfter" ),
                              labels = c("B3m","B1w","A1w","A3m"))
gastricbypassL$visit <- as.numeric(gastricbypassL$time) ## convert to numeric
gastricbypassL$baseline <- gastricbypassL$visit<=2

## chunk 5
gastricbypassL$glucagon <- as.double(scale(gastricbypassL$glucagonAUC))+5

## chunk 6
gastricbypassL$group <- as.numeric(gastricbypassL$id)%%2

## chunk 7
utils::packageVersion("LMMstar")

## * Descriptive statistics

## chunk 8
sss <- summarize(weight+glucagon ~ time, data = gastricbypassL, na.rm = TRUE)
print(sss, digits = 3)

## chunk 9
sss <- summarize(weight ~ time|id, data = gastricbypassL, na.rm = TRUE)
print(sss, digits = 3)

## chunk 10
data(gastricbypassW, package = "LMMstar")
partialCor(weight1 + weight3 ~ 1, data = gastricbypassW)

## chunk 11
partialCor(weight + glucagonAUC ~ group,
           data = gastricbypassL[gastricbypassL$time=="B3m",])

## * Linear mixed model
## ** Classical covariance patterns

## chunk 12
eId.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id, 
               structure = "ID", data = gastricbypassL)
eId.lmm
cat(" covariance structure: \n");sigma(eId.lmm)

## chunk 13
eInd.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id, 
               structure = "IND", data = gastricbypassL)
eInd.lmm
cat(" covariance structure: \n");sigma(eInd.lmm)

## chunk 14
eCS.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id,
               structure = "CS", data = gastricbypassL)
eCS.lmm
cat(" covariance structure: \n");sigma(eCS.lmm)

## chunk 15
eTOE.lmm <- lmm(weight ~ time*group, repetition = ~time|id,
                structure = "TOEPLITZ", data = gastricbypassL)
eTOE.lmm
cat(" correlation structure: \n");cov2cor(sigma(eTOE.lmm))

## *unstructured* covariance matrix:

## chunk 16
eUN.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id,
               structure = "UN", data = gastricbypassL)
eUN.lmm
cat(" covariance structure: \n");sigma(eUN.lmm)

## chunk 17
eSCS.lmm <- lmm(weight ~ time*group,
                repetition = ~time|id, structure = CS(group~1),
                data = gastricbypassL)
eSCS.lmm

## chunk 18
eSUN.lmm <- lmm(weight ~ time*group + glucagon,
                repetition = ~time|id, structure = UN(~group),
                data = gastricbypassL)
eSUN.lmm

## chunk 19
sigma(eSCS.lmm)

## chunk 20
sigma(eSUN.lmm)

## chunk 21
eBCS.lmm <- lmm(weight ~ time*group,repetition = ~time|id,
                structure = CS(~baseline, heterogeneous = FALSE), data = gastricbypassL)
eBCS.lmm
cat(" covariance structure: \n");sigma(eBCS.lmm)

## chunk 23
eBUN.lmm <- lmm(weight ~ time*group, repetition = ~time|id,
                structure = CS(~baseline, heterogeneous = TRUE), data = gastricbypassL)
eBUN.lmm
cat(" covariance structure: \n");sigma(eBUN.lmm)

## ** User-specific covariance patterns

## chunk 24
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

## chunk 25
set.seed(11)
n <- 1000
Y <- rmvnorm(n, mean = rep(0,6), sigma = Rho)
dfL <- reshape2::melt(cbind(id = 1:n, as.data.frame(Y)), id.vars = "id")
dfL$time  <- dfL$variable
dfL <- dfL[order(dfL$id),]
dfL[1:8,]

## chunk 26
myStruct <- CUSTOM(~variable,
                   FCT.sigma = function(p,time,X){rep(p,length(time))}, ## function f
                   init.sigma = c("sigma"=1),
                   FCT.rho = rho.2block, ## function g
                   init.rho = c("rho1"=0.25,"rho2"=0.25,"rho3"=0.25,"rho4"=0.25))

## chunk 27
e.lmmCUSTOM <- lmm(value~time,
                   repetition=~time|id,
                   structure = myStruct,
                   data=dfL,
                   df = FALSE) ## df = FALSE to save computation time
logLik(e.lmmCUSTOM)

## chunk 28
cov2cor(sigma(e.lmmCUSTOM))

## chunk 29
myCS <- CUSTOM(~1,
       FCT.sigma = function(p,time,X){rep(p,length(time))},
       init.sigma = c("sigma"=1),
       FCT.rho = function(p,time,X){matrix(p,length(time),length(time))+diag(1-p,length(time),length(time))},
       init.rho = c("rho"=0.5))

## chunk 30
logLik(lmm(value~time,
           repetition = ~time|id,
           structure = myCS, 
           data = dfL, df = FALSE
           ))

## chunk 31
logLik(lmm(value~time,
           repetition = ~time|id,
           structure = "CS", 
           data = dfL, df = FALSE))

## ** Model output

## chunk 32
summary(eUN.lmm)

## chunk 33
summary(eUN.lmm, hide.mean = TRUE)

## chunk 34
oo <- capture.output(summary(eUN.lmm, hide.fit = TRUE, hide.data = TRUE, hide.cor = TRUE, hide.var = TRUE, hide.sd = TRUE))
cat(sapply(oo[-(1:2)],paste0,"\n"))

## ** Extract estimated coefficients

## chunk 35
coef(eUN.lmm)

## chunk 36
coef(eUN.lmm, effects = "variance")

## chunk 37
coef(eUN.lmm, effects = "variance", transform.k = "sd")

## chunk 38
dummy.coef(eUN.lmm)

## ** Extract estimated coefficient and associated uncertainty

## chunk 39
model.tables(eUN.lmm)

## chunk 40
model.tables(eUN.lmm, effect = "all")

## chunk 41
model.tables(eUN.lmm, columns = c("estimate","p.value"))

## chunk 42
model.tables(eUN.lmm, columns = add("statistic"))

## ** Extract estimated residual variance-covariance structure

## chunk 43
Sigma <- sigma(eUN.lmm)
Sigma

## chunk 44
cov2cor(Sigma)

## chunk 45
sigma(eUN.lmm, cluster = 5)

## chunk 46
newdata <- data.frame(id = "X", time = c("B3m","B1w","A1w","A3m"))
sigma(eUN.lmm, cluster = newdata)

## ** Random effects

## chunk 47
head(coef(eCS.lmm, effects = "ranef"))

## chunk 48
head(coef(eBCS.lmm, effects = "ranef"))

## ** Sum of squares

## chunk 50
df.NNA <- gastricbypassL[order(gastricbypassL$id),]
df.NNA <- df.NNA[!is.na(df.NNA$glucagon),]

## chunk 51
eCS2.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id,
                data = df.NNA, structure = "CS")

## chunk 52
sigma2 <- coef(eCS2.lmm, effect = "variance")^2
tau <- coef(eCS2.lmm, effect = "correlation")*sigma2
omega <- unname(sigma2 - tau)

## chunk 53
df.res <- df.residual(eCS2.lmm)
c(df.res = df.res, SSE = df.res * omega)

## chunk 54
eBeta.lmm <- coef(eCS2.lmm)
eVcov.lmm <- vcov(eCS2.lmm, type.information = "expected")

## chunk 55
attr(model.matrix(eCS2.lmm),"assign")

## chunk 56
SSRstar.time <- eBeta.lmm[2:4] %*% solve(eVcov.lmm[2:4,2:4]) %*% eBeta.lmm[2:4] 
SSRstar.glucagon <- eBeta.lmm[5] %*% solve(eVcov.lmm[5,5]) %*% eBeta.lmm[5] 

## chunk 57
c(time = SSRstar.time * omega,
  glucagon = SSRstar.glucagon * omega)

## ** Proportion of explained variance and partial correlation

## chunk 59
summary(anova(eCS2.lmm), columns = add("partial.r"))

## ** Model diagnostic

## chunk 61
plot(eUN.lmm, type = "scatterplot")

## chunk 62
plot(eUN.lmm, type = "scatterplot2")

## chunk 63
plot(eUN.lmm, type = "correlation", type.residual = "response")
plot(eUN.lmm, type = "correlation", type.residual = "normalized")

## chunk 65
plot(eUN.lmm, type = "qqplot", engine.qqplot = "qqtest")
## Note: the qqtest package to be installed to use the argument engine.plot = "qqtest" 

## chunk 66
eUN.diagW <- residuals(eUN.lmm, type = "normalized", format = "wide")
colnames(eUN.diagW) <- gsub("normalized.","",colnames(eUN.diagW))
head(eUN.diagW)

## chunk 67
eUN.diagL <- residuals(eUN.lmm, type = "normalized", format = "long")
head(eUN.diagL)

## ** Model fit

## chunk 68
library(ggplot2) ## left panel
plot(eUN.lmm, type = "fit", color = "id", ci.alpha = NA, size.text = 20)

## chunk 69
library(emmeans) ## right panel
emmip(eUN.lmm, ~time) + theme(text = element_text(size=20))

## chunk 70
plot(eUN.lmm, type = "fit", at = data.frame(glucagon = 10), color = "glucagon")
## result not shown

## chunk 71
gg <- plot(eUN.lmm, type = "fit", obs.alpha = 0.2, ci = FALSE,plot = FALSE)$plot
gg <- gg + facet_wrap(~id, labeller = label_both)
gg <- gg + theme(axis.text.x=element_text(angle = 90, hjust = 0))
gg

## ** Partial residuals

## chunk 73
gg1 <- plot(eUN.lmm, type = "partial", var = "glucagon", plot = FALSE)$plot
gg2 <- plot(eUN.lmm, type = "partial", var = c("(Intercept)","glucagon"), plot = FALSE)$plot
ggarrange(gg1,gg2)

## chunk 74
df.pres <- residuals(eUN.lmm, type = "partial", var = "glucagon", keep.data = TRUE)
head(df.pres)

## chunk 75
m.pres <- gastricbypassL$weight - model.matrix(~time,gastricbypassL) %*% coef(eUN.lmm)[1:4]
range(df.pres$r.partial - m.pres, na.rm = TRUE)

## chunk 76
eIID.lm <- lm(glucagon ~ time + weight, data = gastricbypassL)
pRes.lm <- residuals(eIID.lm, type = "partial")[,"weight"]

## chunk 77
eIID.lmm <- lmm(glucagon ~ time + weight, data = gastricbypassL)
pRes.lmm <- residuals(eIID.lmm, type = "partial-center", var = "weight")
range(pRes.lm-na.omit(pRes.lmm))

## ** Statistical inference (linear)

## chunk 79
anova(eUN.lmm)

## chunk 80
anova(eUN.lmm, effects = c("timeA1w-timeB1w=0"))

## chunk 81
e.anova <- anova(eUN.lmm, effects = c("timeA1w-timeB1w=0","timeA3m-timeB1w=0"))
summary(e.anova)

## chunk 82
library(multcomp)
summary(anova(eUN.lmm, effects = mcp(time = "Tukey")))

## chunk 83
try(
  anova(eUN.lmm,
        effects = c("log(k).B1w=0","log(k).A1w=0","log(k).A3m=0"))
)

## chunk 84
name.coef <- rownames(confint(eUN.lmm, effects = "all"))
name.varcoef <- grep("^k",name.coef, value = TRUE)
C <- matrix(0, nrow = 3, ncol = length(name.coef), dimnames = list(name.varcoef, name.coef))
diag(C[name.varcoef,name.varcoef]) <- 1
C[,1:9]

## chunk 85
anova(eUN.lmm, effects = C)

## chunk 86
Manova <- rbind(anova(eInd.lmm, effects = "glucagon = 0"),
                anova(eCS.lmm, effects = "glucagon = 0"),
                anova(eUN.lmm, effects = "glucagon = 0"),
                name = c("Ind","CS","UN"))
summary(Manova) 

## ** Statistical inference (non-linear)

## chunk 87
gastricbypassW <- reshape(gastricbypassL[,c("id","time","weight","group")],
                          direction = "wide",
                          timevar = "time", idvar = c("id","group"))
e.ANCOVA <- lm(weight.A1w ~ weight.B1w + group, data = gastricbypassW)
summary(e.ANCOVA)$coef

## chunk 88
gastricbypassL23 <- gastricbypassL[gastricbypassL$visit %in% 2:3,]
gastricbypassL23$time <- droplevels(gastricbypassL23$time)
e.lmmANCOVA <- lmm(weight ~ time+time:group, repetition = ~time|id,
                   data = gastricbypassL23)

## chunk 89
lava::estimate(e.lmmANCOVA, f = function(p){
  c(Y1 = as.double(p["rho(B1w,A1w)"]*p["k.A1w"]),
    X1 = as.double(p["timeA1w:group"]-p["rho(B1w,A1w)"]*p["k.A1w"]*p["timeB1w:group"]))
})

## ** Baseline adjustment

## chunk 90
gastricbypassL$treat <- baselineAdjustment(gastricbypassL, variable = "group",
                                           repetition = ~time|id, constrain = c("B3m","B1w"),
                                           new.level = "none")
table(treat = gastricbypassL$treat, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 91
gastricbypassL$treat2 <- baselineAdjustment(gastricbypassL, variable = "group",
                                            repetition = ~time|id, constrain = c("B3m","B1w"))
table(treat = gastricbypassL$treat2, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 92
gastricbypassL$timeXtreat <- baselineAdjustment(gastricbypassL, variable = "group",
                                                repetition = ~time|id, constrain = c("B3m","B1w"),
                                                collapse.time = ".")

table(treat = gastricbypassL$timeXtreat, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 93
eC.lmm <- lmm(weight ~ timeXtreat, data = gastricbypassL,
              repetition = ~time|id, structure = "UN")
coef(eC.lmm) ## change from baseline

## chunk 94
eC2.lmm <- lmm(weight ~ 0 + timeXtreat, data = gastricbypassL,
              repetition = ~time|id, structure = "UN")
coef(eC2.lmm) ## absolute value

## chunk 95
colnames(model.matrix(weight ~ treat*time, data = gastricbypassL))

## chunk 96
eC3.lmm <- lmm(weight ~ treat2*time, data = gastricbypassL,
               repetition = ~time|id, structure = "UN")

## chunk 97
model.tables(eC3.lmm)

## chunk 98
autoplot(eC3.lmm, color = "group", ci = FALSE, size.text = 20, obs.alpha = 0.1) 

## ** Marginal means

## chunk 99
gastricbypassL$group2 <- as.numeric(gastricbypassL$id) %% 3 == 0
e.group <- lmm(glucagon ~ time*group2, data = gastricbypassL,
               repetition = ~time|id, structure = "UN")

## chunk 100
emmeans(e.group, specs=~time)

## chunk 101
df.pred <- cbind(gastricbypassL, predict(e.group, newdata = gastricbypassL))
summarize(formula = estimate~time, data = df.pred)

## chunk 102
table(group = gastricbypassL$group2, time = gastricbypassL$time)

## chunk 103
mu.group1 <-  as.double(coef(e.group)["(Intercept)"])
mu.group2 <-  as.double(coef(e.group)["(Intercept)"] + coef(e.group)["group2TRUE"])
p.group1 <- 14/20          ; p.group2 <- 6/20
c(emmeans = (mu.group1+mu.group2)/2, predict = mu.group1 * p.group1 + mu.group2 * p.group2)

## chunk 104
emmeans.group <- emmeans(e.group, specs = ~group2|time)
emmeans.group

## chunk 105
epairs.group <- pairs(emmeans.group, reverse = TRUE)
epairs.group

## chunk 106
summary(epairs.group, by = NULL, adjust = "mvt", infer = TRUE)

## chunk 107
summary(pairs(emmeans(eC3.lmm , specs = ~treat2|time), reverse = TRUE), by = NULL)

## ** Predictions

## chunk 108
news <- gastricbypassL[gastricbypassL$id==1,]
news$glucagon <- 0
predict(eUN.lmm, newdata = news)

## chunk 109
X.12 <- model.matrix(formula(eUN.lmm), news)
X.12

## chunk 110
X.12 %*% coef(eUN.lmm)

## chunk 111
newd <- rbind(
  data.frame(id = 1, time = "B3m", weight = coef(eUN.lmm)["(Intercept)"], glucagon = 0),
  data.frame(id = 1, time = "B1w", weight = NA, glucagon = 0),
  data.frame(id = 2, time = "B3m", weight = 100, glucagon = 0),
  data.frame(id = 2, time = "B1w", weight = NA, glucagon = 0)
)
predict(eUN.lmm, newdata = newd, type = "dynamic", keep.newdata = TRUE)

## chunk 112
mu1 <- coef(eUN.lmm)[1]
mu2 <- sum(coef(eUN.lmm)[1:2])
Omega_11 <- sigma(eUN.lmm)["B3m","B3m"]
Omega_21 <- sigma(eUN.lmm)["B1w","B3m"]
as.double(mu2 + Omega_21 * (100 - mu1) / Omega_11)

## * Missing values and imputation
## ** Full information approach

## chunk 113
sss <- summarize(glucagon ~ time, data = gastricbypassL, na.rm = TRUE)
cbind(sss[,1:4], pc = paste0(100 * sss$missing / (sss$missing + sss$observed), "%"))

## chunk 114
vec.pattern <- tapply(as.numeric(is.na(gastricbypassL$glucagon)),
                      INDEX = gastricbypassL$id,
                      FUN = paste, collapse=".")
table(vec.pattern)

## chunk 115
eUN.lmmNA <- lmm(glucagon ~ time,
                 repetition = ~time|id, structure = "UN",
                 data = gastricbypassL)
summary(eUN.lmmNA)

## chunk 116
summary(eUN.lmmNA, hide.mean = TRUE, hide.sd = TRUE, hide.var = TRUE, hide.cor = TRUE, hide.fit = TRUE)

## ** Imputation 

## chunk 117
fitted(eUN.lmmNA, impute = TRUE)

## chunk 118
eData <- fitted(eUN.lmmNA, impute = TRUE, keep.newdata = TRUE)
eData$treat <- eData$treat2 <- eData$timeXtreat <- NULL
eData[eData$id %in% eData[eData$imputed,"id"],]

## chunk 119
ggplot(eData, aes(x=time,y=glucagon, group=id)) + geom_line() + geom_point(aes(color=imputed))

## chunk 121
set.seed(10)
fitted(eUN.lmmNA, impute = TRUE, se = "total")
fitted(eUN.lmmNA, impute = TRUE, se = "total")
fitted(eUN.lmmNA, impute = TRUE, se = "total")

## ** Multiple imputation

## chunk 122
data(gastricbypassW)
colSums(is.na(gastricbypassW))

## chunk 123
library(mice)
gastricbypassW.mice <- mice(gastricbypassW, printFlag = FALSE)
gastricbypassW.NNA <- complete(gastricbypassW.mice, action = "long")
table(gastricbypassW.NNA$.imp)

## chunk 124
e.mlmm <- mlmm(glucagonAUC3~glucagonAUC2+weight2, data=gastricbypassW.NNA, by = ".imp", effects = "weight2=0")
model.tables(e.mlmm)

## chunk 125
model.tables(e.mlmm, method = "pool.rubin")

## chunk 126
e.mice <- with(data=gastricbypassW.mice,exp=lm(glucagonAUC3~glucagonAUC2+weight2))
summary(pool(e.mice))

## * Data generation

## chunk 127
set.seed(10) ## ensure reproductibility
n.obs <- 100
n.times <- 4
mu <- rep(0,4)
gamma <- matrix(0, nrow = n.times, ncol = 10) ## add interaction
gamma[,6] <- c(0,1,1.5,1.5)
dW <- sampleRem(n.obs, n.times = n.times, mu = mu, gamma = gamma, format = "wide")
head(round(dW,3))

## chunk 128
set.seed(10) ## ensure reproductibility
dL <- sampleRem(n.obs, n.times = n.times, mu = mu, gamma = gamma, format = "long")
head(dL)

## * Modifying default options

## chunk 129
LMMstar.options("type.information")

## chunk 130
LMMstar.options(type.information = "expected")

## chunk 131
LMMstar.options(reinitialise = TRUE)

## * R session

## chunk 132
sessionInfo()

## * References
## * Likelihood in a linear mixed model
## ** Log-likelihood
## ** Score
## ** Hessian
## ** Degrees of freedom
## * Likelihood ratio test with the REML criterion

## chunk 133
## data(gastricbypassL, package = "LMMstar")
dfTest <- gastricbypassL
dfTest$glucagon2 <- dfTest$glucagon*2

## chunk 134
eML.lmmUN <- lmm(weight ~ time+glucagon, data = dfTest, repetition = ~time|id, method = "ML")
eML.lmmUN2 <- lmm(weight ~ time+glucagon2, data = dfTest, repetition = ~time|id, method = "ML")

## chunk 135
logLik(eML.lmmUN)
logLik(eML.lmmUN2)

## chunk 136
eREML.lmmUN <- lmm(weight ~ time + glucagon, data = dfTest, repetition = ~time|id, method = "REML")
eREML.lmmUN2 <- lmm(weight ~ time + glucagon2, data = dfTest, repetition = ~time|id, method = "REML")

## chunk 137
logLik(eREML.lmmUN)-logLik(eREML.lmmUN2)
log(2)

## chunk 138
set.seed(15) 
dfTest$ff <- rbinom(NROW(dfTest), size = 1, prob = 0.5)
logLik(lmm(weight ~ time+glucagon, data = dfTest, repetition = ~time|id, method = "REML"))
logLik(lmm(weight ~ time+glucagon*ff, data = dfTest, repetition = ~time|id, method = "REML"))

## chunk 139
logLik(lmm(weight ~ time + glucagon, data = dfTest, repetition = ~time|id, method = "ML"))
logLik(lmm(weight ~ time + glucagon*ff, data = dfTest, repetition = ~time|id, method = "ML"))

## * Sum of squares in a linear mixed model
## *Illustration for a univariate linear model:*

## chunk 140
df.aov <- gastricbypassL[!is.na(gastricbypassL$glucagon),]

## chunk 141
e.lm <- lm(weight ~ time + glucagon, data = df.aov)
car::Anova(e.lm, type = "II")

## chunk 143
e.lmm <- lmm(weight ~ time + glucagon, data = df.aov)

## chunk 144
SSEstar <- crossprod(residuals(e.lmm, type = "normalized"))
c(SSEstar = SSEstar, SSE = SSEstar * sigma(e.lmm))

## chunk 145
df.residual(e.lmm)

## chunk 146
eBeta.lmm <- coef(e.lmm)
eVcov.lmm <- vcov(e.lmm, type.information = "expected")

SSRstar.glucagon <- eBeta.lmm[5] %*% solve(eVcov.lmm[5,5]) %*% eBeta.lmm[5] 
SSRstar.time <- eBeta.lmm[2:4] %*% solve(eVcov.lmm[2:4,2:4]) %*% eBeta.lmm[2:4] 
c(SSR.glucagon = SSRstar.glucagon * sigma(e.lmm),
  SSR.time = SSRstar.time * sigma(e.lmm),
  F.glucagon = SSRstar.glucagon,
  F.time = SSRstar.time/3)

## chunk 147
R2.glucagon <- SSRstar.glucagon/(SSRstar.glucagon+SSEstar)
R2.glucagon

## chunk 148
sign(coef(e.lmm)["glucagon"])*sqrt(R2.glucagon)

## * Equivalent with other R packages
## ** nlme package

## chunk 150
library(nlme)

## chunk 151
eCS.gls <- gls(weight ~ time + glucagon, correlation = corCompSymm(form=~time|id),
               data = gastricbypassL, na.action = na.omit)
eCS.lme <- lme(weight ~ time + glucagon, random = ~1|id,
               data = gastricbypassL, na.action = na.omit)
logLik(eCS.lme)
logLik(eCS.gls)
logLik(eCS.lmm)

## chunk 152
range(coef(eCS.lmm, effects = "ranef")-ranef(eCS.lme))

## chunk 153
eUN.gls <- gls(weight ~ time + glucagon,
               correlation = corSymm(form=~as.numeric(time)|id),
               weights = varIdent(form=~1|time),
               data = gastricbypassL, na.action = na.omit)
logLik(eUN.gls)
logLik(eUN.lmm)

## ** lme4 package

## chunk 154
library(lme4)
library(lmerTest)

## chunk 155
eCS.lmer <- lmer(weight ~ time + glucagon + (1|id),
                 data = gastricbypassL)
logLik(eCS.lmer)
logLik(eCS.lmm)

## chunk 156
range(coef(eCS.lmm, effects = "ranef")-ranef(eCS.lmer)$id)

## chunk 157
eBCS.lmer <- lmer(weight ~ time*group + (1|id/baseline),
                  data = gastricbypassL)
logLik(eBCS.lmer)
logLik(eBCS.lmm)

## chunk 158
eRanefBCS.lmm <- coef(eBCS.lmm, effects = "ranef")
eRanefBCS.lmer <- ranef(eBCS.lmer)
## id
range(eRanefBCS.lmm[,"id"]-eRanefBCS.lmer$id)
## baseline
range(c(eRanefBCS.lmm[,"baseline1"],eRanefBCS.lmm[,"baseline2"])-ranef(eBCS.lmer)$`baseline:id`)

## chunk 159
eUN.lmer <- lmer(weight ~ time + glucagon + (0 + time|id),
                 data = gastricbypassL, control = lmerControl(check.nobs.vs.nRE = "ignore"))
logLik(eUN.lmer)
logLik(eUN.lmm)

## chunk 160
anova(eUN.lmm)

## chunk 161
anova(eUN.lmer)

## chunk 162
eUN2.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id,
                structure = "UN", data = gastricbypassL, type.information = "expected")
suppressWarnings(anova(eUN2.lmm))

## ** effectsize package (\(R^2\) or \(\eta^2\))

## chunk 163
library(effectsize)
eta_squared(eCS.lmer)
cat("\n")

## chunk 164
print(anova(eCS.lmm), columns = add("partial.r"))

## chunk 165
print(anova(eUN.lmm), columns = add("partial.r"))

## chunk 166
eta_squared(eUN.lmer)
cat("\n")

