## chunk 2
library(LMMstar)

## chunk 3
data(gastricbypassL, package = "LMMstar")
head(gastricbypassL)

## chunk 4
utils::packageVersion("LMMstar")

## chunk 5
LMMstar.options(optimizer = "FS")

## * Descriptive statistics

## chunk 6
sss <- summarize(weight+glucagon ~ time, data = gastricbypassL, na.rm = TRUE)
print(sss, digits = 3)

## * Linear mixed model
## ** Modeling tools

## chunk 7
eId.lmm <- lmm(weight ~ time + glucagon,
               repetition = ~time|id, structure = "ID",
               data = gastricbypassL)
eId.lmm
cat(" covariance structure: \n");getVarCov(eId.lmm)

## chunk 8
eInd.lmm <- lmm(weight ~ time + glucagon,
               repetition = ~time|id, structure = "IND",
               data = gastricbypassL)
eInd.lmm
cat(" covariance structure: \n");getVarCov(eInd.lmm)

## chunk 9
eCS.lmm <- lmm(weight ~ time + glucagon,
               repetition = ~time|id, structure = "CS",
               data = gastricbypassL)
eCS.lmm
cat(" covariance structure: \n");getVarCov(eCS.lmm)

## chunk 10
eUN.lmm <- lmm(weight ~ time + glucagon,
               repetition = ~time|id, structure = "UN",
               data = gastricbypassL)
eUN.lmm
cat(" covariance structure: \n");getVarCov(eUN.lmm)

## chunk 11
gastricbypassL$group <- as.numeric(gastricbypassL$id)%%2
eSUN.lmm <- lmm(weight ~ time*group,
                repetition = group~time|id, structure = "UN",
                data = gastricbypassL)
eSUN.lmm
cat(" covariance structure: \n");getVarCov(eSUN.lmm)

## ** Model output

## chunk 12
summary(eCS.lmm)

## ** Extract estimated coefficients
## ** Extract estimated residual variance-covariance structure
## ** Model diagnostic
## ** Model fit
## ** Statistical inference
## *** Model coefficients
## *** Linear combination of the model coefficients
## ** Baseline adjustment
## ** Marginal means
## ** Predictions

## chunk 13
news <- gastricbypassL[gastricbypassL$id==1,]
news$glucagon <- 0
predict(eCS.lmm, newdata = news)

## chunk 14
X.12 <- model.matrix(formula(eCS.lmm), news)
X.12

## chunk 15
X.12 %*% coef(eCS.lmm)

## chunk 16
newd <- rbind(
  data.frame(id = 1, time = "B3_months", weight = coef(eCS.lmm)["(Intercept)"], glucagon = 0),
  data.frame(id = 1, time = "B1_week", weight = NA, glucagon = 0),
  data.frame(id = 2, time = "B3_months", weight = 100, glucagon = 0),
  data.frame(id = 2, time = "B1_week", weight = NA, glucagon = 0)
)
predict(eCS.lmm, newdata = newd, type = "dynamic", keep.newdata = TRUE)

## chunk 17
mu1 <- coef(eCS.lmm)[1]
mu2 <- sum(coef(eCS.lmm)[1:2])
Omega_11 <- getVarCov(eCS.lmm)["B3_months","B3_months"]
Omega_21 <- getVarCov(eCS.lmm)["B1_week","B3_months"]
as.double(mu2 + Omega_21 * (100 - mu1) / Omega_11)

## * Data generation

## chunk 18
set.seed(10) ## ensure reproductibility
n.obs <- 100
n.times <- 4
mu <- rep(0,4)
gamma <- matrix(0, nrow = n.times, ncol = 10) ## add interaction
gamma[,6] <- c(0,1,1.5,1.5)
dW <- sampleRem(n.obs, n.times = n.times, mu = mu, gamma = gamma, format = "wide")
head(round(dW,3))

## chunk 19
set.seed(10) ## ensure reproductibility
dL <- sampleRem(n.obs, n.times = n.times, mu = mu, gamma = gamma, format = "long")
head(dL)

## * Modifying default options

## chunk 20
LMMstar.options("type.information")

## chunk 21
LMMstar.options(type.information = "expected")

## chunk 22
LMMstar.options(reinitialise = TRUE)

## * R session

## chunk 23
sessionInfo()

## * References
## * Likelihood in a linear mixed model
## ** Log-likelihood
## ** Score
## ** Hessian
## ** Degrees of freedom
## * Likelihood ratio test with the REML criterion

## chunk 24
LMMstar.options(optimizer = "FS",
                param.optimizer = c(n.iter = 1000, tol.score = 1e-3, tol.param = 1e-5))

