## chunk 2
library(LMMstar)
library(mvtnorm)
library(nlme)
library(ggpubr)
library(lava)

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

## * Linear mixed model (LMM)
## ** Classical covariance patterns

## chunk 20
eId.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id, 
               structure = "ID", data = dfL)
eId.lmm
cat(" modeled residual variance-covariance: \n");sigma(eId.lmm)

## chunk 21
eInd.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id, 
               structure = "IND", data = dfL)
eInd.lmm
cat(" modeled residual variance-covariance: \n");sigma(eInd.lmm)

## chunk 22
eCS.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id,
               structure = "CS", data = dfL)
eCS.lmm
cat(" modeled residual variance-covariance: \n");sigma(eCS.lmm)

## chunk 23
eTOE.lmm <- lmm(weight ~ time*group, repetition = ~time|id,
                structure = "TOEPLITZ", data = dfL)
eTOE.lmm
cat(" modeled residual correlation: \n");cov2cor(sigma(eTOE.lmm))

## *unstructured* covariance matrix:

## chunk 24
eUN.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id,
               structure = "UN", data = dfL)
eUN.lmm
cat(" modeled residual variance-covariance: \n");sigma(eUN.lmm)

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
                structure = CS(~baseline, type = "homogeneous"), data = dfL)
eBCS.lmm
cat(" modeled residual variance-covariance: \n");sigma(eBCS.lmm)

## chunk 31
eBUN.lmm <- lmm(weight ~ time*group, repetition = ~time|id,
                structure = CS(~baseline, type = "heterogeneous"), data = dfL)
eBUN.lmm
cat(" modeled residual variance-covariance: \n");sigma(eBUN.lmm)

## ** User-specific covariance patterns

## chunk 32
rho.2block <- function(p,n.time,X){
  rho <- matrix(1, nrow = n.time, ncol = n.time)
  rho[1,2] <- rho[2,1] <- rho[4,5] <- rho[5,4] <- p["rho1"]
  rho[1,3] <- rho[3,1] <- rho[4,6] <- rho[6,4] <- p["rho2"]
  rho[2,3] <- rho[3,2] <- rho[5,6] <- rho[6,5] <- p["rho3"]
  rho[4:6,1:3] <- rho[1:3,4:6] <- p["rho4"]
  return(rho)
}
Rho <- rho.2block(p = c(rho1=0.25,rho2=0.5,rho3=0.4,rho4=0.1),
                  n.time = 6)
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
                   FCT.sigma = function(p,n.time,X){rep(p,n.time)}, ## function f
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

## *Comparison to build-in structure*: consider the following model using

## chunk 37
system.time(
  e.lmmDEFAULT.CS <- lmm(value~time,
                         repetition = ~time|id,
                         structure = "CS", 
                         data = dfL2, df = FALSE)
)

## chunk 38
myCS <- CUSTOM(~1,
               FCT.sigma = function(p,n.time,X){rep(p,n.time)},
               init.sigma = c("sigma"=1), 
               FCT.rho = function(p,n.time,X){p+diag(1-p,n.time,n.time)},
               init.rho = c("rho"=0.5))

## chunk 39
system.time(
  e.lmmCUSTOM.CS <- lmm(value~time,
                        repetition = ~time|id,
                        structure = myCS, 
                        data = dfL2, df = FALSE
                        )
)

## chunk 40
logLik(e.lmmDEFAULT.CS)
logLik(e.lmmCUSTOM.CS)


## chunk 41
e.lmmDEFAULT.CS$opt$n.iter
e.lmmCUSTOM.CS$opt$n.iter

## chunk 42
myCS.wD <- CUSTOM(~1,
                  FCT.sigma = function(p,n.time,X){rep(p,n.time)},
                  dFCT.sigma = function(p,n.time,X){list(sigma = rep(1,n.time))},
                  d2FCT.sigma = function(p,n.time,X){list(sigma = rep(0,n.time))},
                  init.sigma = c("sigma"=1),
                  FCT.rho = function(p,n.time,X){p+diag(1-p,n.time,n.time)},
                  dFCT.rho = function(p,n.time,X){list(rho = 1-diag(1,n.time,n.time))},
                  d2FCT.rho = function(p,n.time,X){list(rho = matrix(0,n.time,n.time))},
                  init.rho = c("rho"=0.5))

system.time(
  e.lmmCUSTOMwD.CS <- lmm(value~time,
                          repetition = ~time|id,
                          structure = myCS.wD, 
                          data = dfL2, df = FALSE
                          )
)

## ** Estimation procedure
## *Initialiation*: by default the mean parameters are initialized using

## chunk 43
eCS.lmm.bis <- update(eCS.lmm, control = list(trace = 2))

## chunk 44
init.all <- coef(eCS.lmm, effects = "all")
eCS.lmm.bis <- update(eCS.lmm, control = list(init = init.all, trace = 2))

## chunk 45
init.mean <- coef(eCS.lmm, effects = "mean")
eCS.lmm.bis <- update(eCS.lmm, control = list(init = init.mean, trace = 2))

## chunk 46
init.vcov <- sigma(eCS.lmm)
eCS.lmm.bis <- update(eSCS.lmm, control = list(init = init.vcov, trace = 2))

## *Optimizer*: by default the optimizer is a Newton Raphson algorithm
## ** Model output

## chunk 47
summary(eUN.lmm)

## chunk 48
summary(eUN.lmm, hide.mean = TRUE)

## chunk 49
oo <- capture.output(summary(eUN.lmm, hide.fit = TRUE, hide.data = TRUE, hide.cor = TRUE, hide.var = TRUE, hide.sd = TRUE))
cat(sapply(oo[-(1:2)],paste0,"\n"))

## ** Extract estimated coefficients

## chunk 50
coef(eUN.lmm)

## chunk 51
coef(eUN.lmm, effects = "variance")

## chunk 52
coef(eUN.lmm, effects = "variance", transform.k = "sd")

## chunk 53
dummy.coef(eUN.lmm)

## ** Extract estimated coefficient and associated uncertainty

## chunk 54
model.tables(eUN.lmm)

## chunk 55
model.tables(eUN.lmm, effect = "all")

## chunk 56
model.tables(eUN.lmm, columns = c("estimate","p.value"))

## chunk 57
model.tables(eUN.lmm, columns = add("statistic"))

## ** Extract estimated residual variance-covariance structure

## chunk 58
Sigma <- sigma(eUN.lmm)
Sigma

## chunk 59
cov2cor(Sigma)

## chunk 60
sigma(eUN.lmm, cluster = 5)

## chunk 61
newdata <- data.frame(id = "X", time = c("B3m","B1w","A1w","A3m"))
sigma(eUN.lmm, cluster = newdata)

## ** Random effects

## chunk 62
eRI.lmm <- lmm(weight ~ time + glucagon + (1|id), data = dfL)
eRI.lmm

## chunk 63
eNRI.lmm <- lmm(weight ~ time*group + (1|id/baseline), data = dfL)
eNRI.lmm

## chunk 64
head(ranef(eRI.lmm, format = "wide"))

## chunk 65
head(ranef(eNRI.lmm, format = "wide"))

## chunk 67
ranef(eRI.lmm, effects = "variance", format = "wide", simplify = FALSE)

## chunk 68
ranef(eNRI.lmm, effects = "variance", format = "wide", simplify = FALSE)

## ** Sum of squares

## chunk 69
sigma2 <- coef(eCS.lmm, effect = "variance")^2
tau <- coef(eCS.lmm, effect = "correlation")*sigma2
delta <- unname(sigma2 - tau)

## chunk 70
df.res <- df.residual(eCS.lmm)
SSE <- df.res * delta
c(df.res = df.res, SSE = SSE)

## chunk 71
eBeta.lmm <- coef(eCS.lmm)
eVcov.lmm <- vcov(eCS.lmm, type.information = "expected")

## chunk 72
attr(model.matrix(eCS.lmm),"assign")

## chunk 73
SSRstar.time <- eBeta.lmm[2:4] %*% solve(eVcov.lmm[2:4,2:4]) %*% eBeta.lmm[2:4] 
SSRstar.glucagon <- eBeta.lmm[5] %*% solve(eVcov.lmm[5,5]) %*% eBeta.lmm[5] 

## chunk 74
SSR.time <- as.double(SSRstar.time * delta)
SSR.glucagon <- as.double(SSRstar.glucagon * delta)
c(time = SSR.time, glucagon = SSR.glucagon)

## ** Proportion of explained variance and partial correlation

## chunk 76
c(SSR.time/ (SSR.time + SSE),
  SSR.glucagon/ (SSR.glucagon + SSE))

## chunk 77
eCS.R2 <- partialCor(eCS.lmm, R2 = TRUE)
summary(eCS.R2)

## chunk 78
aCS.aov <- anova(eCS.lmm)$multivariate
setNames(with(aCS.aov, statistic*df.num/(statistic*df.num+df.denom)), aCS.aov$test)

## ** Model diagnostic

## chunk 80
eUN.diagW <- residuals(eUN.lmm, type = "normalized", format = "wide")
colnames(eUN.diagW) <- gsub("normalized.","",colnames(eUN.diagW))
head(eUN.diagW)

## chunk 81
eUN.diagL <- residuals(eUN.lmm, type = "normalized", format = "long", keep.data = TRUE)
head(eUN.diagL)

## chunk 82
plot(eUN.lmm, type = "scatterplot")

## chunk 83
plot(eUN.lmm, type = "scatterplot2")

## chunk 84
plot(eUN.lmm, type = "correlation", type.residual = "response")
plot(eUN.lmm, type = "correlation", type.residual = "normalized")

## chunk 86
plot(eUN.lmm, type = "qqplot", engine.qqplot = "qqtest")
## Note: the qqtest package to be installed to use the argument engine.plot = "qqtest" 

## chunk 87
plot(profile(eUN.lmm, effects = c("sigma","rho(B1w,A1w)")))

## ** Model fit

## chunk 88
library(ggplot2) ## left panel
plot(eUN.lmm, type = "fit", color = "id", ci.alpha = NA, size.text = 20)

## chunk 89
library(emmeans) ## right panel
emmip(eUN.lmm, ~time) + theme(text = element_text(size=20))

## chunk 90
## left panel
plot(eUN.lmm, type = "fit", at = data.frame(glucagon = 10), color = "glucagon") 

## chunk 91
## right panel
gg.spafit <- plot(eUN.lmm, type = "fit", obs.alpha = 0.25, ci = FALSE)$plot

## chunk 92
gg.traj <- gg.spafit + facet_wrap(~id, labeller = label_both)
gg.traj <- gg.traj + theme(axis.text.x=element_text(angle = 90, hjust = 0))
gg.traj

## ** Partial residuals

## chunk 94
gg1 <- plot(eUN.lmm, type = "partial", var = "glucagon")$plot
gg2 <- plot(eUN.lmm, type = "partial", var = c("(Intercept)","glucagon"))$plot
ggarrange(gg1,gg2)

## chunk 96
df.pres <- residuals(eUN.lmm, type = "partial", var = "glucagon", keep.data = TRUE)
head(df.pres)

## chunk 97
m.pres <- dfL$weight - model.matrix(~time,dfL) %*% coef(eUN.lmm)[1:4]
range(df.pres$r.partial - m.pres, na.rm = TRUE)

m.pfit <- model.matrix(~0+glucagon,dfL) %*% coef(eUN.lmm)["glucagon"]
range(df.pres$fitted - m.pfit, na.rm = TRUE)

## ** Statistical inference (linear)

## chunk 98
anova(eUN.lmm)

## chunk 99
anova(eUN.lmm, effects = c("timeA1w-timeB1w=0"))

## chunk 100
e.anova <- anova(eUN.lmm, effects = c("timeA1w-timeB1w=0","timeA3m-timeB1w=0"))
summary(e.anova)

## chunk 101
library(multcomp)
summary(anova(eUN.lmm, effects = mcp(time = "Tukey")))

## chunk 102
try(
  anova(eUN.lmm,
        effects = c("log(k).B1w=0","log(k).A1w=0","log(k).A3m=0"))
)

## chunk 103
name.coef <- rownames(confint(eUN.lmm, effects = "all"))
name.varcoef <- grep("^k",name.coef, value = TRUE)
C <- matrix(0, nrow = 3, ncol = length(name.coef), dimnames = list(name.varcoef, name.coef))
diag(C[name.varcoef,name.varcoef]) <- 1
C[,1:9]

## chunk 104
anova(eUN.lmm, effects = C)

## chunk 105
Manova <- rbind(anova(eInd.lmm, effects = "glucagon = 0", robust = FALSE),
                anova(eCS.lmm, effects = "glucagon = 0", robust = FALSE),
                anova(eUN.lmm, effects = "glucagon = 0", robust = FALSE),
                name = c("Ind","CS","UN"))
summary(Manova) 

## ** Statistical inference (non-linear)

## chunk 106
e.ANCOVA <- lm(weight4 ~ weight1 + group, data = gastricbypassW)
summary(e.ANCOVA)$coef

## chunk 107
dfL14 <- dfL[dfL$visit %in% c(1,4),]
dfL14$time <- droplevels(dfL14$time)
e.lmmANCOVA <- lmm(weight ~ time+time:group, repetition = ~time|id,
                   data = dfL14)

## chunk 108
lava::estimate(e.lmmANCOVA, f = function(p){
  c(Y1 = as.double(p["rho(B3m,A3m)"]*p["k.A3m"]),
    X1 = as.double(p["timeA3m:group"]-p["rho(B3m,A3m)"]*p["k.A3m"]*p["timeB3m:group"]))
})

## ** Baseline adjustment

## chunk 109
gastricbypassL$treat <- baselineAdjustment(gastricbypassL, variable = "group",
                                repetition = ~time|id, constrain = c("B3m","B1w"),
                                new.level = "none")
table(treat = gastricbypassL$treat, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 110
gastricbypassL$treat2 <- baselineAdjustment(gastricbypassL, variable = "group",
                                            repetition = ~time|id, constrain = c("B3m","B1w"))
table(treat = gastricbypassL$treat2, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 111
gastricbypassL$timeXtreat <- baselineAdjustment(gastricbypassL, variable = "group",
                                                repetition = ~time|id, constrain = c("B3m","B1w"),
                                                collapse.time = ".")

table(treat = gastricbypassL$timeXtreat, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 112
eC.lmm <- lmm(weight ~ timeXtreat, data = gastricbypassL,
              repetition = ~time|id, structure = "UN")
coef(eC.lmm) ## change from baseline

## chunk 113
eC2.lmm <- lmm(weight ~ 0 + timeXtreat, data = gastricbypassL,
              repetition = ~time|id, structure = "UN")
coef(eC2.lmm) ## absolute value

## chunk 114
colnames(model.matrix(weight ~ treat*time, data = gastricbypassL))

## chunk 115
eC3.lmm <- lmm(weight ~ treat2*time, data = gastricbypassL,
               repetition = ~time|id, structure = "UN")

## chunk 116
model.tables(eC3.lmm)

## chunk 117
plot(eC3.lmm, color = "group", ci = FALSE, size.text = 20, obs.alpha = 0.1)

## ** Marginal means

## chunk 118
dfL$group2 <- as.numeric(dfL$id) %% 3 == 0
e.group <- lmm(glucagon ~ time*group2, data = dfL,
               repetition = ~time|id, structure = "UN")

## chunk 119
emmeans(e.group, specs=~time)

## chunk 120
df.pred <- predict(e.group, newdata = dfL, keep.newdata = TRUE)
summarize(formula = estimate~time, data = df.pred)

## chunk 121
table(group = dfL$group2, time = dfL$time)

## chunk 122
mu.group1 <-  as.double(coef(e.group)["(Intercept)"])
mu.group2 <-  as.double(coef(e.group)["(Intercept)"] + coef(e.group)["group2TRUE"])
p.group1 <- 14/20          ; p.group2 <- 6/20
c(emmeans = (mu.group1+mu.group2)/2, predict = mu.group1 * p.group1 + mu.group2 * p.group2)

## chunk 123
emmeans.group <- emmeans(e.group, specs = ~group2|time)
emmeans.group

## chunk 124
epairs.group <- pairs(emmeans.group, reverse = TRUE)
epairs.group

## chunk 125
summary(epairs.group, by = NULL, adjust = "mvt", infer = TRUE)

## chunk 126
summary(pairs(emmeans(eC3.lmm , specs = ~treat2|time), reverse = TRUE), by = NULL)

## ** Predictions

## chunk 127
news <- dfL[dfL$id==1,]
news$glucagon <- 0
predict(eUN.lmm, newdata = news)

## chunk 128
X.12 <- model.matrix(formula(eUN.lmm), news)
X.12

## chunk 129
X.12 %*% coef(eUN.lmm)

## chunk 130
newd <- rbind(
  data.frame(id = 1, time = "B3m", weight = coef(eUN.lmm)["(Intercept)"], glucagon = 0),
  data.frame(id = 1, time = "B1w", weight = NA, glucagon = 0),
  data.frame(id = 2, time = "B3m", weight = 100, glucagon = 0),
  data.frame(id = 2, time = "B1w", weight = NA, glucagon = 0)
)
predict(eUN.lmm, newdata = newd, type = "dynamic", keep.newdata = TRUE)

## chunk 131
mu1 <- coef(eUN.lmm)[1]
mu2 <- sum(coef(eUN.lmm)[1:2])
Omega_11 <- sigma(eUN.lmm)["B3m","B3m"]
Omega_21 <- sigma(eUN.lmm)["B1w","B3m"]
as.double(mu2 + Omega_21 * (100 - mu1) / Omega_11)

## * Equivalence with other statistical methods
## ** T-test

## chunk 132
t.test(weight4 ~ group, data = gastricbypassW)

## chunk 133
e.ttest4 <- lmm(weight4 ~ group, structure = IND(~group), 
               data = gastricbypassW, trace = FALSE)
model.tables(e.ttest4)

## chunk 134
e.ttest1 <- lmm(weight1 ~ group, structure = IND(~group), 
                data = gastricbypassW, trace = FALSE)
e.ttest2 <- lmm(weight2 ~ group, structure = IND(~group), 
                data = gastricbypassW, trace = FALSE)
e.ttest3 <- lmm(weight3 ~ group, structure = IND(~group), 
                data = gastricbypassW, trace = FALSE)

## chunk 135
e.mttest <- rbind(anova(e.ttest1, effects = "group=0"),
                  anova(e.ttest2, effects = "group=0"),
                  anova(e.ttest3, effects = "group=0"),
                  anova(e.ttest4, effects = "group=0"))
model.tables(e.mttest, method = "bonferroni")

## chunk 136
e.mttest2 <- mlmm(weight ~ group, structure = IND(~group),
                  data = gastricbypassL, trace = FALSE,
                  effects = "group=0", by = "time", repetition = ~time|id)
model.tables(e.mttest2, method = "single-step2")

## chunk 137
mt.test(weight1+weight2+weight3+weight4~group, data = gastricbypassW)

## ** Linear regression on the change 

## chunk 138
gastricbypassW$changeG41 <- gastricbypassW$glucagonAUC4-gastricbypassW$glucagonAUC1
e.change41 <- lm(changeG41 ~ weight1, data = gastricbypassW)
summary(e.change41)$coef

## chunk 139
gastricbypassL41 <- gastricbypassL[gastricbypassL$visit %in% c(1,4),]
gastricbypassL41$time <- droplevels(gastricbypassL41$time)
gastricbypassL41$weight1 <- gastricbypassW$weight1[gastricbypassL41$id]

e.lmm41 <- lmm(glucagonAUC ~ time + time*weight1,
               repetition =~ time|id, structure = "UN",
               data = gastricbypassL41)
model.tables(e.lmm41)

## chunk 140
index.missing41 <- which(is.na(gastricbypassW$changeG41))
index.missing41

## ** Correlation between changes 

## chunk 141
gastricbypassW$changeG41 <- gastricbypassW$glucagonAUC4-gastricbypassW$glucagonAUC1
gastricbypassW$changeW41 <- gastricbypassW$weight4-gastricbypassW$weight1

## chunk 142
cor.test(gastricbypassW$changeW41, gastricbypassW$changeG41)

## chunk 143
e2.change41 <- lm(changeG41 ~ changeW41, data = gastricbypassW)
summary(e2.change41)$coef

## chunk 144
keep.col <- c("id","weight1","weight4","glucagonAUC1","glucagonAUC4")
gastricbypassL4 <- reshape(gastricbypassW[,keep.col], direction = "long",
                           idvar = "id", varying = 2:5, timevar = "type", v.names = "value")
gastricbypassL4$type <- factor(gastricbypassL4$type, labels = keep.col[-1])
gastricbypassL4 <- gastricbypassL4[order(gastricbypassL4$id),]
head(gastricbypassL4)

## chunk 145
e.lmm4 <- lmm(value ~ type,
              repetition = ~type|id, structure = "UN",
              data = gastricbypassL4)

## chunk 146
sigma.lmm4 <- sigma(e.lmm4)
sigma.lmm4

## chunk 147
Mcon <- cbind(c(-1,1,0,0),c(0,0,-1,1))
sigmeChange.lmm4 <- t(Mcon) %*% sigma.lmm4 %*% Mcon
dimnames(sigmeChange.lmm4) <- list(c("d.weight","d.glucagonAUC"),
                                   c("d.weight","d.glucagonAUC"))
sigmeChange.lmm4

## chunk 148
cov2cor(sigmeChange.lmm4)[1,2]
sigmeChange.lmm4[1,2]/sigmeChange.lmm4[1,1]

## chunk 149
estimate(e.lmm4, function(p){
  Sigma.change <- t(Mcon) %*% sigma(e.lmm4, p = p) %*% Mcon
  c(cor = cov2cor(Sigma.change)[1,2],
    beta = Sigma.change[1,2]/Sigma.change[1,1])
})

## * Missing values and imputation

## chunk 151
sss <- summarize(glucagon ~ time, data = gastricbypassL, na.rm = TRUE)
cbind(sss[,1:4], pc = paste0(100 * sss$missing / (sss$missing + sss$observed), "%"))

## chunk 152
summarizeNA(data = gastricbypassL, repetition = ~ time|id)

## chunk 153
## long format
gastricbypassL32 <- gastricbypassL[gastricbypassL$visit %in% c(3,2),]
gastricbypassL32$time <- droplevels(gastricbypassL32$time)
gastricbypassL32$weight1 <- gastricbypassW$weight1[gastricbypassL32$id]
## wide format
gastricbypassW$changeG32 <- gastricbypassW$glucagonAUC3-gastricbypassW$glucagonAUC2

## ** Full information approach

## chunk 154
e.lmm32 <- lmm(glucagonAUC ~ time + time*weight1,
               repetition =~ time|id, structure = "UN",
               data = gastricbypassL32)
model.tables(e.lmm32)

## chunk 155
e.change32 <- lm(changeG32 ~ weight1, data = gastricbypassW)
summary(e.change32)$coef

## chunk 156
coef(lm(changeG32 ~ weight1, data = gastricbypassW[-c(5,15),]))

## chunk 157
gastricbypassWA <- fitted(e.lmm32, impute = TRUE, format = "wide")
gastricbypassWA$change32 <- gastricbypassWA$glucagonAUC_A1w - gastricbypassWA$glucagonAUC_B1w
gastricbypassWA$weight1 <- gastricbypassW$weight1[match(gastricbypassW$id,gastricbypassWA$id)]
coef(lm(change32 ~ weight1, data = gastricbypassWA))

## ** Complete case approach

## chunk 158
e.lmmCC <- lmmCC(e.change32, repetition = changeG32 ~ glucagonAUC3-glucagonAUC2|id)
model.tables(e.lmmCC)

## chunk 159
summary(e.change32)$coef

## chunk 160
gastricbypassW$changeW32 <- gastricbypassW$weight3 - gastricbypassW$weight2

e2g.change32 <- lm(changeG32 ~ changeW32 + group, data = gastricbypassW)
summary(e2g.change32)$coef

## chunk 161
e2.lmmCC <-  lmmCC(e2g.change32, repetition = list(changeG32 ~ glucagonAUC3-glucagonAUC2|id,
                                                   changeW32 ~ weight3-weight2|id))
model.tables(e2.lmmCC)

## ** Imputation 

## chunk 162
eUN.lmmNA <- lmm(glucagon ~ time, repetition = ~time|id, data = gastricbypassL)
nobs(eUN.lmmNA)

## chunk 163
eData <- fitted(eUN.lmmNA, impute = TRUE, keep.newdata = TRUE)
eData$treat <- eData$treat2 <- eData$timeXtreat <- NULL
eData[eData$id %in% eData[eData$imputed,"id"],]

## chunk 164
ggplot(eData, aes(x=time,y=glucagon, group=id)) + geom_line() + geom_point(aes(color=imputed))

## chunk 166
set.seed(10)
index.na <- which(is.na(gastricbypassL$glucagonAUC))
fitted(eUN.lmmNA, impute = TRUE, se.impute = "total")[index.na]
fitted(eUN.lmmNA, impute = TRUE, se.impute = "total")[index.na]
fitted(eUN.lmmNA, impute = TRUE, se.impute = "total")[index.na]

## ** Multiple imputation

## chunk 167
data(gastricbypassW, package = "LMMstar")
colSums(is.na(gastricbypassW))

## chunk 168
library(mice)
set.seed(10)
gastricbypassW.mice <- mice(gastricbypassW, m = 5, printFlag = FALSE)
gastricbypassW.NNA <- complete(gastricbypassW.mice, action = "long")
table(gastricbypassW.NNA$.imp)

## chunk 169
e.mlmm <- mlmm(glucagonAUC3~glucagonAUC2+weight2, data=gastricbypassW.NNA,
               by = ".imp", effects = "weight2=0", trace = FALSE)
model.tables(e.mlmm)

## chunk 170
model.tables(e.mlmm, method = "pool.rubin")

## chunk 171
e.mice <- with(data=gastricbypassW.mice,exp=lm(glucagonAUC3~glucagonAUC2+weight2))
summary(pool(e.mice))

## chunk 172
plot(e.mlmm, method = c("pool.rubin","none"))

## * Data generation

## chunk 174
set.seed(10) ## ensure reproductibility
n.obs <- 100
n.times <- 4
mu <- rep(0,4)
gamma <- matrix(0, nrow = n.times, ncol = 10) ## add interaction
gamma[,6] <- c(0,1,1.5,1.5)
dW <- sampleRem(n.obs, n.times = n.times, mu = mu, gamma = gamma, format = "wide")
head(round(dW,3))

## chunk 175
set.seed(10) ## ensure reproductibility
dL <- sampleRem(n.obs, n.times = n.times, mu = mu, gamma = gamma, format = "long")
head(dL)

## * Modifying default options

## chunk 176
LMMstar.options("type.information")

## chunk 177
LMMstar.options(type.information = "expected")

## chunk 178
LMMstar.options(reinitialise = TRUE)

## * R session

## chunk 179
sessionInfo()

## * References
## * Likelihood in a linear mixed model
## ** Log-likelihood
## ** Score
## ** Hessian
## ** Degrees of freedom
## * Likelihood ratio test with the REML criterion

## chunk 180
## data(dfL, package = "LMMstar")
dfTest <- dfL
dfTest$glucagon2 <- dfTest$glucagon*2

## chunk 181
eML.lmmUN <- lmm(weight ~ time+glucagon, data = dfTest, repetition = ~time|id, method = "ML")
eML.lmmUN2 <- lmm(weight ~ time+glucagon2, data = dfTest, repetition = ~time|id, method = "ML")

## chunk 182
logLik(eML.lmmUN)
logLik(eML.lmmUN2)

## chunk 183
eREML.lmmUN <- lmm(weight ~ time + glucagon, data = dfTest, repetition = ~time|id, method = "REML")
eREML.lmmUN2 <- lmm(weight ~ time + glucagon2, data = dfTest, repetition = ~time|id, method = "REML")

## chunk 184
logLik(eREML.lmmUN)-logLik(eREML.lmmUN2)
log(2)

## chunk 185
set.seed(5) 
dfTest$ff <- rbinom(NROW(dfTest), size = 1, prob = 0.5)
logLik(lmm(weight ~ time+glucagon, data = dfTest, repetition = ~time|id, method = "REML"))
logLik(lmm(weight ~ time+glucagon*ff, data = dfTest, repetition = ~time|id, method = "REML"))

## chunk 186
logLik(lmm(weight ~ time + glucagon, data = dfTest, repetition = ~time|id, method = "ML"))
logLik(lmm(weight ~ time + glucagon*ff, data = dfTest, repetition = ~time|id, method = "ML"))

## * Sum of squares in a linear mixed model
## *Illustration for a univariate linear model:*

## chunk 187
df.aov <- dfL[!is.na(dfL$glucagon),]

## chunk 188
e.lm <- lm(weight ~ time + glucagon, data = df.aov)
car::Anova(e.lm, type = "II")

## chunk 190
e.lmm <- lmm(weight ~ time + glucagon, data = df.aov)

## chunk 191
SSEstar <- crossprod(residuals(e.lmm, type = "normalized"))
c(SSEstar = SSEstar, SSE = SSEstar * sigma(e.lmm))

## chunk 192
df.residual(e.lmm)

## chunk 193
eBeta.lmm <- coef(e.lmm)
eVcov.lmm <- vcov(e.lmm, type.information = "expected")

SSRstar.glucagon <- eBeta.lmm[5] %*% solve(eVcov.lmm[5,5]) %*% eBeta.lmm[5] 
SSRstar.time <- eBeta.lmm[2:4] %*% solve(eVcov.lmm[2:4,2:4]) %*% eBeta.lmm[2:4] 
c(SSR.glucagon = SSRstar.glucagon * sigma(e.lmm),
  SSR.time = SSRstar.time * sigma(e.lmm),
  F.glucagon = SSRstar.glucagon,
  F.time = SSRstar.time/3)

## chunk 194
R2.glucagon <- SSRstar.glucagon/(SSRstar.glucagon+SSEstar)
R2.glucagon

## chunk 195
sign(coef(e.lmm)["glucagon"])*sqrt(R2.glucagon)

## chunk 196
summary(partialCor(e.lmm, R2 = TRUE))

## * Equivalence with other R packages
## ** nlme package

## chunk 197
library(nlme)

## chunk 198
eCS.gls <- gls(weight ~ time + glucagon, correlation = corCompSymm(form=~time|id),
               data = dfL, na.action = na.omit)
eCS.lme <- lme(weight ~ time + glucagon, random = ~1|id,
               data = dfL, na.action = na.omit)
logLik(eCS.lme)
logLik(eCS.gls)
logLik(eCS.lmm)

## chunk 199
range(ranef(eCS.lmm)$estimate-ranef(eCS.lme))

## chunk 200
eUN.gls <- gls(weight ~ time + glucagon,
               correlation = corSymm(form=~as.numeric(time)|id),
               weights = varIdent(form=~1|time),
               data = dfL, na.action = na.omit)
logLik(eUN.gls)
logLik(eUN.lmm)

## ** lme4 package

## chunk 201
library(lme4)
library(lmerTest)

## chunk 202
eRI.lmer <- lmer(weight ~ time + glucagon + (1|id),
                 data = dfL)
logLik(eRI.lmer)
logLik(eRI.lmm)

## chunk 203
range(ranef(eRI.lmm)$estimate-ranef(eRI.lmer)$id)

## chunk 204
eNRI.lmer <- lmer(weight ~ time*group + (1|id/baseline),
                  data = dfL)
logLik(eNRI.lmer)
logLik(eNRI.lmm)

## chunk 205
eRanefNRI.lmm <- ranef(eNRI.lmm)
eRanefNRI.lmer <- ranef(eNRI.lmer)
## id
range(eRanefNRI.lmm[eRanefNRI.lmm$variable=="id","estimate"]-eRanefNRI.lmer$id)
## baseline
range(eRanefNRI.lmm[eRanefNRI.lmm$variable!="id","estimate"]-ranef(eNRI.lmer)$`baseline:id`)

## chunk 206
eUN.lmer <- lmer(weight ~ time + glucagon + (0 + time|id),
                 data = dfL, control = lmerControl(check.nobs.vs.nRE = "ignore"))
logLik(eUN.lmer)
logLik(eUN.lmm)

## chunk 207
anova(eUN.lmm)

## chunk 208
anova(eUN.lmer)

## chunk 209
eUN2.lmm <- lmm(weight ~ time + glucagon, repetition = ~time|id,
                structure = "UN", data = dfL, type.information = "expected")
suppressWarnings(anova(eUN2.lmm))

## chunk 210
data("Penicillin")
eCRI.lmer <- lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin)
logLik(eCRI.lmer)

## chunk 211
Penicillin$index <- paste(Penicillin$sample,Penicillin$plate,sep=".")
Penicillin$id <- 1

eCRI.lmm <- lmm(diameter ~ 1 + (1|plate) + (1|sample), data = Penicillin)
logLik(eCRI.lmm)

## chunk 212
range(ranef(eCRI.lmm)$estimate-rbind(ranef(eCRI.lmer)$plate,ranef(eCRI.lmer)$sample))

## ** mmrm package

## chunk 213
library(mmrm)
e.mmrm <- mmrm(
  formula = FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID),
  data = fev_data
)

## chunk 214
e.lmm <- lmm(
  formula = FEV1 ~ RACE + SEX + ARMCD * AVISIT,
  repetition = ~ AVISIT | USUBJID, structure = "UN",
  data = fev_data, type.information = "expected"
)

## chunk 215
logLik(e.mmrm) - logLik(e.lmm)
range(coef(e.mmrm) - coef(e.lmm))
range(vcov(e.mmrm) - vcov(e.lmm))

## ** effectsize package (\(R^2\) or \(\eta^2\))

## chunk 219
library(effectsize)
eta_squared(eCS.lmer)
cat("\n")

## chunk 220
eCS.Wald <- anova(eCS.lmm)$multivariate
eCS.Wald$df.num*eCS.Wald$statistic/(eCS.Wald$df.num*eCS.Wald$statistic+eCS.Wald$df.denom)

## chunk 221
eUN.Wald <- anova(eUN.lmm)$multivariate
eUN.Wald$df.num*eUN.Wald$statistic/(eUN.Wald$df.num*eUN.Wald$statistic+eUN.Wald$df.denom)

## chunk 222
eta_squared(eUN.lmer)
cat("\n")

## ** MuMIn package (\(R^2\))

## chunk 223
library(MuMIn)
r.squaredGLMM(eCS.lmer)
cat("\n")

## chunk 224
sigmaW <- sigma(eCS.lmm)[1,1]-sigma(eCS.lmm)[1,2]

## chunk 225
sigmaB <- sigma(eCS.lmm)[1,2]

## chunk 226
sigma2_XB <- var(fitted(eCS.lmm))

## chunk 227
c(R2m = sigma2_XB/(sigmaW + sigmaB + sigma2_XB),
  R2c = (sigma2_XB + sigmaB)/(sigmaW + sigmaB + sigma2_XB))

## chunk 233
eIID.lm <- lm(weight ~ time + glucagon, data = dfL)
pRes.lm <- residuals(eIID.lm, type = "partial")
head(pRes.lm)

## chunk 234
eIID.lmm <- lmm(weight ~ time + glucagon, data = dfL)
head(residuals(eIID.lmm, type = "partial", var = "glucagon"))

## chunk 235
m.pres2 <- dfL$weight - cbind(model.matrix(~time,dfL), mean(dfL$glucagon)) %*% coef(eIID.lmm)
range(pRes.lm[,"glucagon"] - m.pres2, na.rm = TRUE)

## chunk 236
eIID.lmm <- lmm(weight ~ time + glucagon, data = dfL)
pRes.lmm <- residuals(eIID.lmm, type = "partial-center", var = "glucagon")
range(pRes.lm[,"glucagon"]-pRes.lmm)

