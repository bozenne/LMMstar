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

## chunk 8
LMMstar.options(optimizer = "FS")

## * Descriptive statistics

## chunk 9
sss <- summarize(weight+glucagon ~ time, data = gastricbypassL, na.rm = TRUE)
print(sss, digits = 3)

## chunk 10
sss <- summarize(weight ~ time|id, data = gastricbypassL, na.rm = TRUE)
print(sss, digits = 3)

## chunk 11
partialCor(list(weight~group, glucagonAUC~group),
           data = gastricbypassL[gastricbypassL$time=="B3m",])

## * Linear mixed model
## ** Classical covariance patterns

## chunk 12
eId.lmm <- lmm(weight ~ time + glucagon,
               repetition = ~time|id, structure = "ID",
               data = gastricbypassL)
eId.lmm
cat(" covariance structure: \n");sigma(eId.lmm)

## chunk 13
eInd.lmm <- lmm(weight ~ time + glucagon,
               repetition = ~time|id, structure = "IND",
               data = gastricbypassL)
eInd.lmm
cat(" covariance structure: \n");sigma(eInd.lmm)

## chunk 14
eCS.lmm <- lmm(weight ~ time + glucagon,
               repetition = ~time|id, structure = "CS",
               data = gastricbypassL)
eCS.lmm
cat(" covariance structure: \n");sigma(eCS.lmm)

## chunk 15
eSCS.lmm <- lmm(weight ~ time*group,
                repetition = ~time|id, structure = CS(group~1),
                data = gastricbypassL)
eSCS.lmm
cat(" covariance structure: \n");sigma(eSCS.lmm)

## chunk 16
eBCS.lmm <- lmm(weight ~ time*group,
                repetition = ~time|id, structure = CS(~baseline, heterogeneous = FALSE),
                data = gastricbypassL)
eBCS.lmm
cat(" covariance structure: \n");sigma(eBCS.lmm)

## chunk 18
eBUN.lmm <- lmm(weight ~ time*group,
                repetition = ~time|id, structure = CS(~baseline),
                data = gastricbypassL)
eBUN.lmm
cat(" covariance structure: \n");sigma(eBUN.lmm)

## chunk 19
eUN.lmm <- lmm(weight ~ time + glucagon,
               repetition = ~time|id, structure = "UN",
               data = gastricbypassL)
eUN.lmm
cat(" covariance structure: \n");sigma(eUN.lmm)

## chunk 20
eSUN.lmm <- lmm(weight ~ time*group + glucagon,
                repetition = ~time|id, structure = UN(~group),
                data = gastricbypassL)
eSUN.lmm
cat(" covariance structure: \n");sigma(eSUN.lmm)

## ** Model output

## chunk 21
summary(eUN.lmm)

## chunk 22
summary(eUN.lmm, hide.mean = TRUE)

## chunk 23
oo <- capture.output(summary(eUN.lmm, hide.fit = TRUE, hide.data = TRUE, hide.cor = TRUE, hide.var = TRUE, hide.sd = TRUE))
cat(sapply(oo[-(1:2)],paste0,"\n"))

## ** Extract estimated coefficients

## chunk 24
coef(eUN.lmm)

## chunk 25
coef(eUN.lmm, effects = "variance")

## chunk 26
coef(eUN.lmm, effects = "variance", transform.k = "sd")

## chunk 27
dummy.coef(eUN.lmm)

## ** Extract estimated coefficient and associated uncertainty

## chunk 28
model.tables(eUN.lmm)

## chunk 29
model.tables(eUN.lmm, effect = "all") ## not shown

## ** Extract estimated residual variance-covariance structure

## chunk 30
sigma(eUN.lmm)

## chunk 31
sigma(eUN.lmm, cluster = 5)

## chunk 32
newdata <- data.frame(id = "X", time = c("B3m","B1w","A1w","A3m"))
sigma(eUN.lmm, cluster = newdata)

## ** Model diagnostic

## chunk 33
plot(eUN.lmm, type = "scatterplot")

## chunk 34
plot(eUN.lmm, type = "scatterplot2")

## chunk 35
plot(eUN.lmm, type = "correlation", type.residual = "response")
plot(eUN.lmm, type = "correlation", type.residual = "normalized")

## chunk 37
plot(eUN.lmm, type = "qqplot", engine.qqplot = "qqtest")
## Note: the qqtest package to be installed to use the argument engine.plot = "qqtest" 

## chunk 38
eUN.diagW <- residuals(eUN.lmm, type = "normalized", format = "wide")
colnames(eUN.diagW) <- gsub("normalized.","",colnames(eUN.diagW))
head(eUN.diagW)

## chunk 39
eUN.diagL <- residuals(eUN.lmm, type = "normalized", format = "long")
head(eUN.diagL)

## ** Model fit

## chunk 40
library(ggplot2) ## left panel
plot(eUN.lmm, type = "fit", color = "id", ci.alpha = NA, size.text = 20)

## chunk 41
library(emmeans) ## right panel
emmip(eUN.lmm, ~time) + theme(text = element_text(size=20))

## chunk 42
plot(eUN.lmm, type = "fit", at = data.frame(glucagon = 10), color = "glucagon")

## chunk 43
  gg <- plot(eUN.lmm, type = "fit", obs.alpha = 0.2, ci = FALSE,plot = FALSE)$plot
  gg <- gg + facet_wrap(~id, labeller = label_both)
  gg <- gg + theme(axis.text.x=element_text(angle = 90, hjust = 0))
  gg
ggsave(gg + theme(text = element_text(size=20)), filename = "figures/fit-autoplot-indiv.pdf", width = 12)

## chunk 44
gg1 <- plot(eUN.lmm, type = "partial", var = "glucagon", plot = FALSE)$plot
gg2 <- plot(eUN.lmm, type = "partial", var = c("(Intercept)","glucagon"), plot = FALSE)$plot
ggarrange(gg1,gg2)

## chunk 45
df.pres <- residuals(eUN.lmm, type = "partial", var = "glucagon", keep.data = TRUE)
m.pres <- gastricbypassL$weight - model.matrix(~time,gastricbypassL) %*% coef(eUN.lmm)[1:4]
range(df.pres$r.partial - m.pres, na.rm = TRUE)

## ** Statistical inference (linear)

## chunk 47
anova(eUN.lmm)

## chunk 48
anova(eUN.lmm, effects = c("timeA1w-timeB1w=0"))

## chunk 49
e.anova <- anova(eUN.lmm, effects = c("timeA1w-timeB1w=0","timeA3m-timeB1w=0"))
summary(e.anova)

## chunk 50
library(multcomp)
summary(anova(eUN.lmm, effects = mcp(time = "Tukey")))

## chunk 51
try(
  anova(eUN.lmm,
        effects = c("log(k).B1w=0","log(k).A1w=0","log(k).A3m=0"))
)

## chunk 52
name.coef <- rownames(confint(eUN.lmm, effects = "all"))
name.varcoef <- grep("^k",name.coef, value = TRUE)
C <- matrix(0, nrow = 3, ncol = length(name.coef), dimnames = list(name.varcoef, name.coef))
diag(C[name.varcoef,name.varcoef]) <- 1
C

## chunk 53
anova(eUN.lmm, effects = C)

## chunk 54
Manova <- rbind(anova(eInd.lmm, effects = "glucagon = 0"),
                anova(eCS.lmm, effects = "glucagon = 0"),
                anova(eUN.lmm, effects = "glucagon = 0"),
                name = c("Ind","CS","UN"))
summary(Manova) 

## ** Statistical inference (non-linear)

## chunk 55
gastricbypassW <- reshape(gastricbypassL[,c("id","time","weight","group")],
                          direction = "wide",
                          timevar = "time", idvar = c("id","group"))
e.ANCOVA <- lm(weight.A1w ~ weight.B1w + group, data = gastricbypassW)
summary(e.ANCOVA)$coef

## chunk 56
e.lmmANCOVA <- lmm(weight ~ time+time:group, repetition = ~time|id,
                   data = gastricbypassL[gastricbypassL$visit %in% 2:3,])

## chunk 57
lava::estimate(e.lmmANCOVA, f = function(p){
  c(Y1 = as.double(p["rho(B1w,A1w)"]*p["k.A1w"]),
    X1 = as.double(p["timeA1w:group"]-p["rho(B1w,A1w)"]*p["k.A1w"]*p["timeB1w:group"]))
})

## ** Baseline adjustment

## chunk 58
gastricbypassL$treat <- baselineAdjustment(gastricbypassL, variable = "group",
                                           repetition = ~time|id, constrain = c("B3m","B1w"),
                                           new.level = "none")
table(treat = gastricbypassL$treat, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 59
gastricbypassL$treat2 <- baselineAdjustment(gastricbypassL, variable = "group",
                                            repetition = ~time|id, constrain = c("B3m","B1w"))
table(treat = gastricbypassL$treat2, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 60
gastricbypassL$timeXtreat <- baselineAdjustment(gastricbypassL, variable = "group",
                                                repetition = ~time|id, constrain = c("B3m","B1w"),
                                                collapse.time = ".")

table(treat = gastricbypassL$timeXtreat, time = gastricbypassL$time, group = gastricbypassL$group)

## chunk 61
eC.lmm <- lmm(weight ~ timeXtreat, data = gastricbypassL,
              repetition = ~time|id, structure = "UN")
coef(eC.lmm) ## change from baseline

## chunk 62
eC2.lmm <- lmm(weight ~ 0 + timeXtreat, data = gastricbypassL,
              repetition = ~time|id, structure = "UN")
coef(eC2.lmm) ## absolute value

## chunk 63
colnames(model.matrix(weight ~ treat*time, data = gastricbypassL))

## chunk 64
eC3.lmm <- lmm(weight ~ treat2*time, data = gastricbypassL,
               repetition = ~time|id, structure = "UN")

## chunk 65
model.tables(eC3.lmm)

## chunk 66
autoplot(eC3.lmm, color = "group", ci = FALSE, size.text = 20, obs.alpha = 0.1) 

## ** Marginal means

## chunk 67
e.group <- lmm(weight ~ time*group, data = gastricbypassL,
               repetition = ~time|id, structure = "UN")

## chunk 68
emmeans(e.group, specs=~time)

## chunk 69
df.pred <- cbind(gastricbypassL, predict(e.group, newdata = gastricbypassL))
summarize(formula = estimate~time, data = df.pred)

## chunk 70
table(group = gastricbypassL$group, time = gastricbypassL$time)

## chunk 71
mu.group1 <-  as.double(coef(e.group)["(Intercept)"])
mu.group2 <-  as.double(coef(e.group)["(Intercept)"] + coef(e.group)["group"])
p.group1 <- 14/20          ; p.group2 <- 6/20
c(emmeans = (mu.group1+mu.group2)/2, predict = mu.group1 * p.group1 + mu.group2 * p.group2)

## chunk 72
emmeans.group <- emmeans(e.group, specs = ~group|time)
emmeans.group

## chunk 73
epairs.group <- pairs(emmeans.group, reverse = TRUE)
epairs.group

## chunk 74
summary(epairs.group, by = NULL, adjust = "mvt", infer = TRUE)

## chunk 75
summary(pairs(emmeans(eC3.lmm , specs = ~treat2|time), reverse = TRUE), by = NULL)

## ** Predictions

## chunk 76
news <- gastricbypassL[gastricbypassL$id==1,]
news$glucagon <- 0
predict(eUN.lmm, newdata = news)

## chunk 77
X.12 <- model.matrix(formula(eUN.lmm), news)
X.12

## chunk 78
X.12 %*% coef(eUN.lmm)

## chunk 79
newd <- rbind(
  data.frame(id = 1, time = "B3m", weight = coef(eUN.lmm)["(Intercept)"], glucagon = 0),
  data.frame(id = 1, time = "B1w", weight = NA, glucagon = 0),
  data.frame(id = 2, time = "B3m", weight = 100, glucagon = 0),
  data.frame(id = 2, time = "B1w", weight = NA, glucagon = 0)
)
predict(eUN.lmm, newdata = newd, type = "dynamic", keep.newdata = TRUE)

## chunk 80
mu1 <- coef(eUN.lmm)[1]
mu2 <- sum(coef(eUN.lmm)[1:2])
Omega_11 <- sigma(eUN.lmm)["B3m","B3m"]
Omega_21 <- sigma(eUN.lmm)["B1w","B3m"]
as.double(mu2 + Omega_21 * (100 - mu1) / Omega_11)

## ** Missing values and imputation

## chunk 81
sss <- summarize(glucagon ~ time, data = gastricbypassL, na.rm = TRUE)
cbind(sss[,1:4], pc = paste0(100 * sss$missing / (sss$missing + sss$observed), "%"))

## chunk 82
vec.pattern <- tapply(as.numeric(is.na(gastricbypassL$glucagon)),
                      INDEX = gastricbypassL$id,
                      FUN = paste, collapse=".")
table(vec.pattern)

## chunk 83
eUN.lmmNA <- lmm(glucagon ~ time,
                 repetition = ~time|id, structure = "UN",
                 data = gastricbypassL)
summary(eUN.lmmNA, hide.fit = TRUE,
        hide.cor = TRUE, hide.sd = TRUE, hide.mean = TRUE)

## chunk 84
fitted(eUN.lmmNA, impute = TRUE)

## chunk 85
eData <- fitted(eUN.lmmNA, impute = TRUE, keep.newdata = TRUE)
eData$treat <- eData$treat2 <- eData$timeXtreat <- NULL
eData[eData$id %in% eData[eData$imputed,"id"],]

## chunk 86
ggplot(eData, aes(x=time,y=glucagon, group=id)) + geom_line() + geom_point(aes(color=imputed))

## chunk 88
set.seed(10)
fitted(eUN.lmmNA, impute = TRUE, se = "total")
fitted(eUN.lmmNA, impute = TRUE, se = "total")
fitted(eUN.lmmNA, impute = TRUE, se = "total")

## * User-specific covariance patterns

## chunk 89
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

## chunk 90
set.seed(11)
n <- 1000
Y <- rmvnorm(n, mean = rep(0,6), sigma = Rho)
dfL <- reshape2::melt(cbind(id = 1:n, as.data.frame(Y)), id.vars = "id")
dfL$time  <- dfL$variable
dfL <- dfL[order(dfL$id),]
dfL[1:8,]

## chunk 91
e.lmmCUSTOM <- lmm(value~time,
                   repetition=~time|id,
                   structure=CUSTOM(~variable,
                                    FCT.sigma = function(p,time,X){rep(p,length(time))}, ## function f
                                    init.sigma = c("sigma"=1),
                                    FCT.rho = rho.2block, ## function g
                                    init.rho = c("rho1"=0.25,"rho2"=0.25,"rho3"=0.25,"rho4"=0.25)),
                   data=dfL, control = list(optimizer = "FS"),
                   df = FALSE) ## df = FALSE to save computation time
logLik(e.lmmCUSTOM)

## chunk 92
cov2cor(sigma(e.lmmCUSTOM))

## chunk 93
logLik(lmm(value~time,
           repetition=~time|id,
           structure=CUSTOM(~1,
                            FCT.sigma = function(p,time,X){rep(p,length(time))},
                            init.sigma = c("sigma"=1),
                            FCT.rho = function(p,time,X){matrix(p,length(time),length(time))+diag(1-p,length(time),length(time))},
                            init.rho = c("rho"=0.5)), 
           data=dfL, control = list(optimizer = "FS"),
           df = FALSE))

## chunk 94
logLik(lmm(value~time,
           repetition=~time|id,
           structure="CS", 
           data=dfL, control = list(optimizer = "FS"),
           df = FALSE))

## * Data generation

## chunk 95
set.seed(10) ## ensure reproductibility
n.obs <- 100
n.times <- 4
mu <- rep(0,4)
gamma <- matrix(0, nrow = n.times, ncol = 10) ## add interaction
gamma[,6] <- c(0,1,1.5,1.5)
dW <- sampleRem(n.obs, n.times = n.times, mu = mu, gamma = gamma, format = "wide")
head(round(dW,3))

## chunk 96
set.seed(10) ## ensure reproductibility
dL <- sampleRem(n.obs, n.times = n.times, mu = mu, gamma = gamma, format = "long")
head(dL)

## * Modifying default options

## chunk 97
LMMstar.options("type.information")

## chunk 98
LMMstar.options(type.information = "expected")

## chunk 99
LMMstar.options(reinitialise = TRUE)

## * R session

## chunk 100
sessionInfo()

## * References
## * Likelihood in a linear mixed model
## ** Log-likelihood
## ** Score
## ** Hessian
## ** Degrees of freedom
## * Likelihood ratio test with the REML criterion

## chunk 101
LMMstar.options(optimizer = "FS",
                param.optimizer = c(n.iter = 1000, tol.score = 1e-3, tol.param = 1e-5))

## chunk 102
## data(gastricbypassL, package = "LMMstar")
dfTest <- gastricbypassL
dfTest$glucagon2 <- dfTest$glucagon*2

## chunk 103
logLik(lmm(weight ~ glucagon, data = dfTest, structure = UN(~time|id), method = "ML"))
logLik(lmm(weight ~ glucagon2, data = dfTest, structure = UN(~time|id), method = "ML"))

## chunk 104
logLik(lmm(weight ~ glucagon, data = dfTest, structure = UN(~time|id), method = "REML"))
logLik(lmm(weight ~ glucagon2, data = dfTest, structure = UN(~time|id), method = "REML"))
log(2)

## chunk 105
set.seed(1)
dfTest$ff <- rbinom(NROW(dfTest), size = 1, prob = 0.5)
logLik(lmm(weight ~ glucagon, data = dfTest, structure = UN(~time|id), method = "REML"))
logLik(lmm(weight ~ glucagon*ff, data = dfTest, structure = UN(~time|id), method = "REML"))

## chunk 106
logLik(lmm(weight ~ glucagon, data = dfTest, structure = UN(~time|id), method = "ML"))
logLik(lmm(weight ~ glucagon*ff, data = dfTest, structure = UN(~time|id), method = "ML"))

