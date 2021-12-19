## chunk 2
library(LMMstar)

## chunk 3
data(gastricbypassL, package = "LMMstar")
head(gastricbypassL)

## chunk 4
gastricbypassL$time <- factor(gastricbypassL$time,
                              levels = c("3 months before surgery", "1 week before surgery",
                                                            "1 week after surgery", "3 months after surgery" ),
                              labels = c("B3_months","B1_week","A1_week","A3_months"))

## chunk 5
gastricbypassL$glucagon <- as.double(scale(gastricbypassL$glucagon))

## chunk 6
utils::packageVersion("LMMstar")

## chunk 7
LMMstar.options(optimizer = "FS")

## * Descriptive statistics

## chunk 8
sss <- summarize(weight+glucagon ~ time, data = gastricbypassL, na.rm = TRUE)
print(sss, digits = 3)

## chunk 9
sss <- summarize(weight ~ time|id, data = gastricbypassL, na.rm = TRUE)
print(sss, digits = 3)

## * Linear mixed model
## ** Modeling tools

## chunk 10
eId.lmm <- lmm(weight ~ time + glucagon,
               repetition = ~time|id, structure = "ID",
               data = gastricbypassL)
eId.lmm
cat(" covariance structure: \n");getVarCov(eId.lmm)

## chunk 11
eInd.lmm <- lmm(weight ~ time + glucagon,
               repetition = ~time|id, structure = "IND",
               data = gastricbypassL)
eInd.lmm
cat(" covariance structure: \n");getVarCov(eInd.lmm)

## chunk 12
eCS.lmm <- lmm(weight ~ time + glucagon,
               repetition = ~time|id, structure = "CS",
               data = gastricbypassL)
eCS.lmm
cat(" covariance structure: \n");getVarCov(eCS.lmm)

## chunk 13
eUN.lmm <- lmm(weight ~ time + glucagon,
               repetition = ~time|id, structure = "UN",
               data = gastricbypassL)
eUN.lmm
cat(" covariance structure: \n");getVarCov(eUN.lmm)

## chunk 14
gastricbypassL$group <- as.numeric(gastricbypassL$id)%%2
eSUN.lmm <- lmm(weight ~ time*group,
                repetition = group~time|id, structure = "UN",
                data = gastricbypassL)
eSUN.lmm
cat(" covariance structure: \n");getVarCov(eSUN.lmm)

## ** Model output

## chunk 15
summary(eUN.lmm)

## chunk 16
summary(eUN.lmm, hide.mean = TRUE)

## chunk 17
oo <- capture.output(summary(eUN.lmm, hide.fit = TRUE, hide.data = TRUE, hide.cor = TRUE, hide.var = TRUE, hide.sd = TRUE))
cat(sapply(oo[-(1:2)],paste0,"\n"))

## ** Extract estimated coefficients

## chunk 18
coef(eUN.lmm)

## chunk 19
coef(eUN.lmm, effects = "variance")

## chunk 20
coef(eUN.lmm, effects = "variance", transform.k = "sd")

## chunk 21
dummy.coef(eUN.lmm)

## ** Extract estimated coefficient and associated uncertainty

## chunk 22
model.tables(eUN.lmm)

## chunk 23
model.tables(eUN.lmm, effect = "all") ## not shown

## ** Extract estimated residual variance-covariance structure

## chunk 24
getVarCov(eUN.lmm)

## chunk 25
getVarCov(eUN.lmm, individual = 5)

## chunk 26
newdata <- data.frame(id = "X", time = c("B3_months","B1_week","A1_week","A3_months"))
getVarCov(eUN.lmm, individual = newdata)

## ** Model diagnostic

## chunk 27
plot(eUN.lmm, type = "scatterplot")

## chunk 28
plot(eUN.lmm, type = "scatterplot2")

## chunk 29
plot(eUN.lmm, type = "correlation", type.residual = "response")
plot(eUN.lmm, type = "correlation", type.residual = "normalized")

## chunk 31
plot(eUN.lmm, type = "qqplot", engine.qqplot = "qqtest")
## Note: the qqtest package to be installed to use the argument engine.plot = "qqtest" 

## chunk 32
eUN.diagW <- residuals(eUN.lmm, type = "normalized", format = "wide")
colnames(eUN.diagW) <- gsub("normalized.","",colnames(eUN.diagW))
head(eUN.diagW)

## chunk 33
eUN.diagL <- residuals(eUN.lmm, type = "normalized", format = "long")
head(eUN.diagL)

## ** Model fit

## chunk 34
library(ggplot2) ## left panel
plot(eUN.lmm, type = "fit", color = "id", ci.alpha = NA, size.text = 20)

## chunk 35
library(emmeans) ## right panel
emmip(eUN.lmm, ~time) + theme(text = element_text(size=20))

## chunk 36
plot(eUN.lmm, type = "fit", at = data.frame(glucagon = 10), color = "glucagon")

## chunk 37
gg <- plot(eUN.lmm, type = "fit", obs.alpha = 0.2, ci = FALSE,plot = FALSE)$plot
gg <- gg + facet_wrap(~id, labeller = label_both)
gg <- gg + theme(axis.text.x=element_text(angle = 90, hjust = 0))
gg

## ** Statistical inference

## chunk 38
anova(eUN.lmm)

## chunk 39
anova(eUN.lmm, effects = c("timeA1_week-timeB1_week=0"), ci = TRUE)

## chunk 40
anova(eUN.lmm, effects = c("timeA1_week-timeB1_week=0",
                           "timeA3_months-timeB1_week=0"), ci = TRUE)

## chunk 41
library(multcomp)
anova(eUN.lmm, effects = mcp(time = "Tukey"), ci = TRUE)

## chunk 42
try(
  anova(eUN.lmm,
        effects = c("log(k).B1_week=0","log(k).A1_week=0","log(k).A3_months=0"))
)

## chunk 43
name.coef <- rownames(confint(eUN.lmm, effects = "all", backtransform = FALSE))
name.varcoef <- grep("log(k)",name.coef, value = TRUE, fixed = TRUE)
C <- matrix(0, nrow = 3, ncol = length(name.coef), dimnames = list(name.varcoef, name.coef))
diag(C[name.varcoef,name.varcoef]) <- 1
C

## chunk 44
anova(eUN.lmm, effects = C)

## ** Baseline adjustment

## chunk 45
gastricbypassL$group <- c("1","2")[as.numeric(gastricbypassL$id) %in% 15:20 + 1]

## chunk 46
gastricbypassL$treat <- baselineAdjustment(gastricbypassL, variable = "group",
                                           repetition = ~time|id, constrain = c("B3_months","B1_week"),
                                           new.level = "none")
table(treat = gastricbypassL$treat, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 47
colnames(model.matrix(weight ~ treat*time, data = gastricbypassL))

## chunk 48
eC.lmm <- lmm(weight ~ treat*time, data = gastricbypassL,
              repetition = ~time|id, structure = "UN")

## chunk 49
coef(eC.lmm, effects = "mean")

## chunk 50
autoplot(eC.lmm, color = "group", ci = FALSE, size.text = 20)

## chunk 51
gastricbypassL$treat2 <- baselineAdjustment(gastricbypassL, variable = "group",
                                            repetition = ~time|id, constrain = c("B3_months","B1_week"))
table(treat = gastricbypassL$treat2, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 52
eC2.lmm <- suppressWarnings(lmm(weight ~ treat2*time, data = gastricbypassL,
                                repetition = ~time|id, structure = "UN"))

## chunk 53
model.tables(eC2.lmm)

## chunk 54
eC3.lmm <- suppressWarnings(lmm(weight ~ 0+treat2:time, data = gastricbypassL,
                                repetition = ~time|id, structure = "UN"))
model.tables(eC3.lmm) ## equivalent to dummy.coef(eC2.lmm)

## chunk 55
eC4.lmm <- suppressWarnings(lmm(weight ~ treat2:time, data = gastricbypassL,
                                repetition = ~time|id, structure = "UN"))
model.tables(eC4.lmm)

## ** Marginal means

## chunk 56
e.group <- lmm(weight ~ time*group, data = gastricbypassL,
               repetition = ~time|id, structure = "UN")

## chunk 57
emmeans(e.group, specs=~time)

## chunk 58
df.pred <- cbind(gastricbypassL, predict(e.group, newdata = gastricbypassL))
summarize(formula = estimate~time, data = df.pred)

## chunk 59
table(group = gastricbypassL$group, time = gastricbypassL$time)

## chunk 60
mu.group1 <-  as.double(coef(e.group)["(Intercept)"])
mu.group2 <-  as.double(coef(e.group)["(Intercept)"] + coef(e.group)["group2"])
p.group1 <- 14/20          ; p.group2 <- 6/20
c(emmeans = (mu.group1+mu.group2)/2, predict = mu.group1 * p.group1 + mu.group2 * p.group2)

## chunk 61
emmeans.group <- emmeans(e.group, specs = ~group|time)
emmeans.group

## chunk 62
epairs.group <- pairs(emmeans.group, reverse = TRUE)
epairs.group

## chunk 63
summary(epairs.group, by = NULL, adjust = "mvt", infer = TRUE)

## chunk 64
summary(pairs(emmeans(eC2.lmm , specs = ~treat2|time), reverse = TRUE), by = NULL)

## ** Predictions

## chunk 65
news <- gastricbypassL[gastricbypassL$id==1,]
news$glucagon <- 0
predict(eUN.lmm, newdata = news)

## chunk 66
X.12 <- model.matrix(formula(eUN.lmm), news)
X.12

## chunk 67
X.12 %*% coef(eUN.lmm)

## chunk 68
newd <- rbind(
  data.frame(id = 1, time = "B3_months", weight = coef(eUN.lmm)["(Intercept)"], glucagon = 0),
  data.frame(id = 1, time = "B1_week", weight = NA, glucagon = 0),
  data.frame(id = 2, time = "B3_months", weight = 100, glucagon = 0),
  data.frame(id = 2, time = "B1_week", weight = NA, glucagon = 0)
)
predict(eUN.lmm, newdata = newd, type = "dynamic", keep.newdata = TRUE)

## chunk 69
mu1 <- coef(eUN.lmm)[1]
mu2 <- sum(coef(eUN.lmm)[1:2])
Omega_11 <- getVarCov(eUN.lmm)["B3_months","B3_months"]
Omega_21 <- getVarCov(eUN.lmm)["B1_week","B3_months"]
as.double(mu2 + Omega_21 * (100 - mu1) / Omega_11)

## ** Missing values and imputation

## chunk 70
sss <- summarize(glucagon ~ time, data = gastricbypassL, na.rm = TRUE)
cbind(sss[,1:4], pc = paste0(100 * sss$missing / (sss$missing + sss$observed), "%"))

## chunk 71
vec.pattern <- tapply(as.numeric(is.na(gastricbypassL$glucagon)),
                      INDEX = gastricbypassL$id,
                      FUN = paste, collapse=".")
table(vec.pattern)

## chunk 72
eUN.lmmNA <- lmm(glucagon ~ time,
                 repetition = ~time|id, structure = "UN",
                 data = gastricbypassL)
summary(eUN.lmmNA, hide.fit = TRUE,
        hide.cor = TRUE, hide.sd = TRUE, hide.mean = TRUE)

## chunk 73
fitted(eUN.lmmNA, impute = TRUE)

## chunk 74
eData <- fitted(eUN.lmmNA, impute = TRUE, keep.newdata = TRUE)
eData[eData$id %in% eData[eData$imputed,"id"],]

## chunk 75
ggplot(eData, aes(x=time,y=glucagon, group=id)) + geom_line() + geom_point(aes(color=imputed))

## chunk 77
set.seed(10)
fitted(eUN.lmmNA, impute = TRUE, se = "total")
fitted(eUN.lmmNA, impute = TRUE, se = "total")
fitted(eUN.lmmNA, impute = TRUE, se = "total")

## * Data generation

## chunk 78
set.seed(10) ## ensure reproductibility
n.obs <- 100
n.times <- 4
mu <- rep(0,4)
gamma <- matrix(0, nrow = n.times, ncol = 10) ## add interaction
gamma[,6] <- c(0,1,1.5,1.5)
dW <- sampleRem(n.obs, n.times = n.times, mu = mu, gamma = gamma, format = "wide")
head(round(dW,3))

## chunk 79
set.seed(10) ## ensure reproductibility
dL <- sampleRem(n.obs, n.times = n.times, mu = mu, gamma = gamma, format = "long")
head(dL)

## * Modifying default options

## chunk 80
LMMstar.options("type.information")

## chunk 81
LMMstar.options(type.information = "expected")

## chunk 82
LMMstar.options(reinitialise = TRUE)

## * R session

## chunk 83
sessionInfo()

## * References
## * Likelihood in a linear mixed model
## ** Log-likelihood
## ** Score
## ** Hessian
## ** Degrees of freedom
## * Likelihood ratio test with the REML criterion

## chunk 84
LMMstar.options(optimizer = "FS",
                param.optimizer = c(n.iter = 1000, tol.score = 1e-3, tol.param = 1e-5))

## chunk 85
## data(gastricbypassL, package = "LMMstar")
dfTest <- gastricbypassL
dfTest$glucagon2 <- dfTest$glucagon*2

## chunk 86
logLik(lmm(weight ~ glucagon, data = dfTest, structure = UN(~time|id), method = "ML"))
logLik(lmm(weight ~ glucagon2, data = dfTest, structure = UN(~time|id), method = "ML"))

## chunk 87
logLik(lmm(weight ~ glucagon, data = dfTest, structure = UN(~time|id), method = "REML"))
logLik(lmm(weight ~ glucagon2, data = dfTest, structure = UN(~time|id), method = "REML"))
log(2)

## chunk 88
set.seed(1)
dfTest$ff <- rbinom(NROW(dfTest), size = 1, prob = 0.5)
logLik(lmm(weight ~ glucagon, data = dfTest, structure = UN(~time|id), method = "REML"))
logLik(lmm(weight ~ glucagon*ff, data = dfTest, structure = UN(~time|id), method = "REML"))

## chunk 89
logLik(lmm(weight ~ glucagon, data = dfTest, structure = UN(~time|id), method = "ML"))
logLik(lmm(weight ~ glucagon*ff, data = dfTest, structure = UN(~time|id), method = "ML"))

