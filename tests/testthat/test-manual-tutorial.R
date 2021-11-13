### test-manual-tutorial.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 13 2021 (16:47) 
## Version: 
## Last-Updated: nov 13 2021 (18:30) 
##           By: Brice Ozenne
##     Update #: 5
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
    library(lattice) 
    library(psych)   
    library(emmeans) 

    library(LMMstar)
}

context("Check lmm on Julie tutorial")
LMMstar.options(optimizer = "gls", method.numDeriv = "simple", precompute.moments = TRUE)

test.tutorial <- FALSE

if(test.tutorial){

## ** section 3.1
data("gastricbypassW", package = "LMMstar")
wide <- gastricbypassW

## ** section 3.2
str(wide)
summary(wide)

## ** section 3.5
#### Transform data from wide to long format #####
long <- reshape(wide, 
                direction='long', 
                idvar='id', 
                  varying=list(
                    c('weight1','weight2','weight3','weight4'),
                    c('glucagonAUC1', 'glucagonAUC2', 'glucagonAUC3', 'glucagonAUC4')
                    ), 
                  v.names=c('weight','glucagonAUC'),
                timevar='visit')

# Make a categorical version of the time variable:
time.names <- c('-3 month','-1 week','+1 week','+3 month')
long$time <- factor(long$visit, labels=time.names)

## * section 5

long$time <- relevel(long$time, ref="-3 month")

fit.main <- lmm(weight~time,
                repetition=~visit|id,
                structure="UN",
                df=TRUE,
                data=long)

## ** Display of the parametrisation
library(ggplot2)

ggParam <- autoplot(fit.main, ci = FALSE, plot = FALSE, mean.size = c(4, 1.5), size.text = 20)$plot
ggParam <- ggParam + coord_cartesian(ylim = c(0,1.2*coef(fit.main)[1]), xlim = c(0.5,4))
## ggParam <- ggParam + geom_curve(aes(x = 1, y = 0, xend = 1, yend = 125),
##                       arrow = arrow(length = unit(0.03, "npc"), type="closed"),
##                       colour = "#EC7014", size = 1.2, angle = -90)
ggParam <- ggParam + geom_curve(aes(x = 1, y = 0, xend = 1, yend = 0.99*coef(fit.main)[1]),
                      arrow = arrow(length = unit(0.03, "npc"), type="closed"),
                      colour = "#EC7014", size = 1, angle = -90, curvature = -0.25)
ggParam <- ggParam + geom_curve(aes(x = 1, y = 1.01*coef(fit.main)[1], xend = 2, yend = 1.02*(coef(fit.main)[1]+coef(fit.main)[2])),
                      arrow = arrow(length = unit(0.03, "npc"), type="closed"),
                      colour = "blue", size = 1, angle = -90, curvature = -0.25)
ggParam <- ggParam + geom_curve(aes(x = 1, y = 1.01*coef(fit.main)[1], xend = 3, yend = 1.02*(coef(fit.main)[1]+coef(fit.main)[3])),
                      arrow = arrow(length = unit(0.03, "npc"), type="closed"),
                      colour = "purple", size = 1, angle = -90, curvature = -0.25)
ggParam <- ggParam + geom_curve(aes(x = 1, y = 1.01*coef(fit.main)[1], xend = 4, yend = 1.02*(coef(fit.main)[1]+coef(fit.main)[4])),
                      arrow = arrow(length = unit(0.03, "npc"), type="closed"),
                      colour = "darkgreen", size = 1, angle = -90, curvature = -0.25)
ggParam <- ggParam + geom_text(mapping = aes(x = 0.9, y = coef(fit.main)[1]/2, label = "beta[1]"), parse = TRUE,colour = "#EC7014", size = 8)
ggParam <- ggParam + geom_text(mapping = aes(x = 1.5, y = 0.9*coef(fit.main)[1], label = "beta[2]"), parse = TRUE,colour = "blue", size = 8)
ggParam <- ggParam + geom_text(mapping = aes(x = 2.3, y = 1*coef(fit.main)[1], label = "beta[3]"), parse = TRUE,colour = "purple", size = 8)
ggParam <- ggParam + geom_text(mapping = aes(x = 3.5, y = 1.05*coef(fit.main)[1], label = "beta[4]"), parse = TRUE,colour = "darkgreen", size = 8)
ggParam <- ggParam + scale_x_discrete(breaks=1:4,
                              labels=sort(unique(long$time)))
## ggParam <- ggParam + scale_y_continuous(breaks=as.double(c(0,50,100,sort(coef(fit.main)[1]+c(0,coef(fit.main)[-1])),150)),
##                               labels=list(bquote(0),bquote(50),bquote(100),bquote(mu[1]),bquote(mu[2]),bquote(mu[3]),bquote(mu[4]),bquote(150)))
ggParam <- ggParam + scale_y_continuous(breaks=as.double(c(0,50,100,coef(fit.main)[1]+c(0,coef(fit.main)[-1]),150)),
                              labels=c(0,50,100,
                                       expression(mu[1]),
                                       expression(mu[2]),
                                       expression(mu[3]),
                                       expression(mu[4]),
                                       150))
ggParam <- ggParam + theme(panel.grid.minor = element_blank())
ggParam

## ggsave(ggParam, filename = "figures/gg-explanation-table.png", width = 10)

## ** Dynamic predictions
newd <- rbind(data.frame(id = 1:100, time = "-3 month", visit = 1, weight = seq(100,175,length.out = 100)),
              data.frame(id = 1:100, time = "-1 week", visit = 2, weight = NA))


## sigma(lm(weight2 ~ weight1, data = gastricbypassW))
## head(dfFit)
## sqrt(prod(coef(fit.main, effects = "all")[c("sigma","k.2")])^2*(1-coef(fit.main, effects = "all")[c("rho(1,2)")]^2))

dfFit <- cbind(weight = newd[newd$visit==1,"weight"],
               res = predict(fit.main, newdata = newd, type = "dynamic", se = "res"),
               total = predict(fit.main, newdata = newd, type = "dynamic", se = "total")
               )
dfData <- reshape2::dcast(data = long[long$time %in% c("-3 month","-1 week"),], formula = id~time, value.var = "weight")
colnames(dfData) <- c("id","w1","w2")

ggDyn <- ggplot()
ggDyn <- ggDyn + geom_point(data = dfData, aes(x=w1,y=w2))
ggDyn <- ggDyn + geom_line(data = dfFit, aes(x=weight,y=res.estimate))
ggDyn <- ggDyn + geom_line(data = dfFit, aes(x=weight,y=res.upper, color = "residual variance"), linetype = 2, size = 1.1)
ggDyn <- ggDyn + geom_line(data = dfFit, aes(x=weight,y=res.lower, color = "residual variance"), linetype = 2, size = 1.1)
ggDyn <- ggDyn + geom_line(data = dfFit, aes(x=weight,y=total.upper, color = "residual variance and estimation uncertainty"), linetype = 3, size = 1.1)
ggDyn <- ggDyn + geom_line(data = dfFit, aes(x=weight,y=total.lower, color = "residual variance and estimation uncertainty"), linetype = 3, size = 1.1)
ggDyn <- ggDyn + labs(x = "weight 3 months before", y = "weight 1 week before", color = "95% prediction limit accounting for")
ggDyn <- ggDyn + theme(legend.position = "bottom", legend.direction = "vertical",
                 text = element_text(size=15),
                 axis.line = element_line(size = 1.25),
                 axis.ticks = element_line(size = 2),
                 axis.ticks.length=unit(.25, "cm"))
ggDyn
## ggsave(ggDyn, filename = "figures/gg-dynamic-prediction.png", width = 10)

## ** residual plot
df.res <- residuals(fit.main, format = "long", type = "all", keep.data = TRUE)
df.res$predicted <- df.res$weight  - df.res$r.response 
df.res[df.res$id==2,c("id","visit","weight","predicted","r.response","r.pearson","r.studentized","r.scaled")]

df.allres <- residuals(fit.main, type = "all", keep.data = TRUE)
gg.res1 <- ggplot(df.allres, aes(x=time,y = weight, group = id, color = id)) + geom_line() + geom_point()
gg.res2 <- ggplot(df.allres, aes(x=time,y = r.response, group = id, color = id)) + geom_line() + geom_point() + ylab("raw residuals")
gg.res3 <- ggplot(df.allres, aes(x=time,y = r.studentized, group = id, color = id)) + geom_line() + geom_point() + ylab("studentized residuals")
gg.res4 <- ggplot(df.allres, aes(x=time,y = r.scaled, group = id, color = id)) + geom_line() + geom_point() + ylab("scaled residuals")
library(ggpubr)
gg.res <- ggarrange(gg.res1,gg.res2,gg.res3,gg.res4, legend = "none")
## ggsave(gg.res, filename = "figures/gg-residuals.png", width = 10)


## * section 6
#### Linear mixed model analysis of a single group study ####

### Main analysis of changes in mean since baseline ####

## ** section 6.1
## Set reference point (intercept) for time factor:
long$time <- relevel(long$time, ref="-3 month")

fit.main <- lmm(weight~time,
                repetition=~visit|id,
                structure="UN",
                df=TRUE,
                data=long)

## Model summary, loads of information
summary(fit.main) 

## Display fitted values
autoplot(fit.main)

## Extract estimates and confidence intervals
confint(fit.main)

## Extract covariance matrix:
getVarCov(fit.main)

# F-test:
anova(fit.main, ci = TRUE)

## round(coef(fit.main, effects = "correlation"),3)
## round(coef(fit.main, effects = "variance", transform.k = "sd"),2)

## Predictions
# Make a dataset with covariate values for prediction:
pred <- long[,c('visit','time')]
# Reduce to one of each value:
pred <- unique(pred)
# Add predicted means to the dataframe:
pred <- cbind(pred, predict(fit.main, newdata=pred))

# Plot predicted means
xyplot(estimate~time, data=pred, type='b')
xyplot(se~time, data=pred, type='b')

# OPTIONAL: Picture including 95% CIs (emmeans-package)
emmip(fit.main, ~time, CIs=TRUE, xlab='Time', ylab='Mean weight')


## Residual diagnostics
par(mfrow=c(2,2))
plot(fitted(fit.main), residuals(fit.main, type='studentized'), main='Studentized residuals')
abline(h=0)
qqnorm(residuals(fit.main, type='studentized'))
abline(0,1)
# Note: Studentized residuals are slightly skewed. 
plot(fitted(fit.main), residuals(fit.main, type='scaled'), main='Scaled residuals')
abline(h=0)
qqnorm(residuals(fit.main, type='scaled'))
abline(0,1)
# Note: Scaled residuals look good.
residuals(fit.main, format = "long", type = "normalized", plot = "scatterplot")
residuals(fit.main, format = "long", type = "normalized", plot = "qqplot", engine.qqplot = "qqtest")
residuals(fit.main, format = "wide", type = "normalized", plot = "qqplot", engine.qqplot = "qqtest")

## ** section 6.5
##### Alternative syntax: Means over time ##########
# Fit the model without an intercept using -1 in the model formula:
fit.means <- lmm(weight~-1+time,
                 repetition=~visit|id,
                 structure="UN",
                 data=long,
                 df=TRUE)

# Extract estimated means and confidence intervals:
confint(fit.means)


## * Appendix D
#### Alternative analysis: Log-transformed weights ####

# Scatterplot matrices with or without effects:
plot(log2(weight.replicates))
pairs.panels(log2(weight.replicates), smooth=FALSE, ellipses=TRUE)

# Add log2-transformed weights to the long data:
long$log2weight <- log2(long$weight)

# Spaghettiplot
xyplot(log2weight~time, data=long, group=id, type='b')

# Fit the linear mixed model
fit.log <- lmm(log2weight~time,
               repetition=~visit|id,
               structure="UN",
               df=TRUE,
               data=long)

# Estimates and CIs on log2-scale:
confint(fit.log)
# Back-transformed estimates and CIs:
2^confint(fit.log)

# Save predicted population means on log2-scale and plot them
pred.log <- long[,c('visit','time')]
pred.log <- unique(pred.log)
pred.log <- cbind(pred.log, predict(fit.log, newdata=pred.log))

xyplot(estimate~time, type='b', data=pred.log) # mean on log2-scale
xyplot(2^estimate~time, type='b', data=pred.log) # back-transformed, i.e. geometric means or medians

# Residual diagnostics:
plot(fitted(fit.log), residuals(fit.log, type='studentized'))
abline(h=0)
qqnorm(residuals(fit.log, type='studentized'))
abline(0,1)
plot(fitted(fit.log), residuals(fit.log, type='scaled'))
abline(h=0)
qqnorm(residuals(fit.log, type='scaled'))
abline(0,1)
# Note: slightly better fit. 
}


##----------------------------------------------------------------------
### test-manual-tutorial.R ends here
