### test-manual-armd.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Dec 19 2021 (17:07) 
## Version: 
## Last-Updated: jul 15 2024 (10:24) 
##           By: Brice Ozenne
##     Update #: 45
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
    library(ggplot2) 

    library(LMMstar)
} 

context("Check lmm on Basic Statistic course")
LMMstar.options(optimizer = "FS", method.numDeriv = "simple", precompute.moments = TRUE, df = TRUE)
test.practical <- FALSE

## * Data
## Load data
data(armd.wide, package = "nlmeU")

## and move to long format
armd.long <- reshape(armd.wide, direction = "long", idvar = "subject",
                     timevar = "week.num", times = c(0,4,12,24,52), varying = paste0("visual",c(0,4,12,24,52)), 
                     v.names = "visual")
armd.long$week <- as.factor(armd.long$week.num)

## * summarize
test_that("summarize", {
    if(test.practical==FALSE){skip('Not run to save time in the check')}
    armd.sum <- summarize(visual ~ week * treat.f | subject, ## value ~ group
                          data = armd.long, ## dataset
                          na.rm = TRUE) ## additional argument
    armd.cor <- correlate(visual ~ week * treat.f | subject, ## value ~ group
                          data = armd.long, ## dataset
                          use = "pairwise.complete.obs") ## additional argument
    expect_equivalent(armd.cor[,,"visual.treat.f=Placebo"],
                      cor(armd.wide[armd.wide$treat.f=="Placebo",paste0("visual",c(0,4,12,24,52))], use = "pairwise"),
                      tol = 1e-4)

})

## * mixed model
test_that("lmm 2 times", {
    if(test.practical==FALSE){skip('Not run to save time in the check')}

    e052.lmm <- lmm(visual ~ treat.f*week,
                    repetition = ~week|subject, 
                    data = armd.long[armd.long$week %in% c("0","52"),])
    model.tables(e052.lmm)

    ## figure
    
    ## chunk 20
    expect_equal(levels(e052.lmm)$reference, c(treat.f = "Placebo", week = "0"))

    ## chunk 21
    c(placebo.0 = as.double(coef(e052.lmm)["(Intercept)"]),
      placebo.52 = sum(coef(e052.lmm)[c("(Intercept)","week52")]),
      active.0 = sum(coef(e052.lmm)[c("(Intercept)","treat.fActive")]),
      active.52 = sum(coef(e052.lmm)))

    ## chunk 22
    effects(e052.lmm, variable = "treat.f")

    armd.114 <- armd.long[armd.long$subject=="114" & armd.long$week %in% c("0","52"),]
    armd.114

    ## chunk 23
    sigma0 <- coef(e052.lmm, effects="variance")["sigma"]
    sigma52 <- sigma0*coef(e052.lmm, effects="variance")["k.52"]
    rho <- coef(e052.lmm, effects="correlation")["rho(0,52)"]
    alpha0 <- coef(e052.lmm)["(Intercept)"]
    alpha52 <- alpha0 + coef(e052.lmm)["week52"]

    ## chunk 24
    expect_equivalent(predict(e052.lmm, newdata = armd.114, type = "dynamic"),
                      c(NA,unname(alpha52 + rho * sigma52/sigma0 * (armd.114$visual[1]-alpha0))),
                      tol = 1e-4)

})

test_that("lmm 4 times", {
    if(test.practical==FALSE){skip('Not run to save time in the check')}
    armd.long$week.num <- as.numeric(as.character(armd.long$week))

    eLin.lmm <- lmm(visual ~ week + week.num:treat.f,
                    repetition = ~ week | subject, structure = "UN",
                    data = armd.long)
    model.tables(eLin.lmm)

    levels(eLin.lmm)$reference
    summary(effects(eLin.lmm, variable = "treat.f", type = "identity"))
    summary(effects(eLin.lmm, variable = "treat.f", type = "difference"))
    ## plot(eLin.lmm, ci = FALSE, obs.alpha = 0.1)

    eFlex.lmm <- lmm(visual ~ week*treat.f,
                     repetition = ~ week | subject, structure = "UN",
                     data = armd.long)
    model.tables(eFlex.lmm)
    summary(effects(eFlex.lmm, variable = "treat.f", type = "identity"))
    summary(effects(eFlex.lmm, variable = "treat.f", type = "difference"))
    ## plot(eFlex.lmm, ci = FALSE, obs.alpha = 0.1)

    armd.long.imp <- fitted(eFlex.lmm, type = "outcome", keep.data = TRUE)
    gg <- autoplot(eFlex.lmm, obs.alpha = 0.1, ci = FALSE)$plot
    gg <- gg + geom_point(data = armd.long.imp[armd.long.imp$subject %in% c("114","167"),,drop=FALSE], aes(x = week, y = visual, group = subject, color = treat.f, shape = impute), size = 4)
    gg <- gg + geom_line(data = armd.long.imp[armd.long.imp$subject %in% c("114","167"),,drop=FALSE], aes(x = week, y = visual, group = subject, color = treat.f))
    gg
    
    pLin <- fitted(eLin.lmm, newdata = "unique", keep.data = TRUE, se = TRUE)
    pFlex <- fitted(eFlex.lmm, newdata = "unique", keep.data = TRUE, se = TRUE)
    armdU <- rbind(cbind(model = "linear", pLin[,c("week","treat.f","visual","lower","upper")]),
                   cbind(model = "flexible", pFlex[,c("week","treat.f","visual","lower","upper")]))

    gg.comp <- ggplot(armdU, aes(x = week, y = visual, color = treat.f,
                                 group =  interaction(treat.f,model)))
    gg.comp <- gg.comp + geom_errorbar(aes(ymin = lower, ymax = upper),
                                       width = 0.2, alpha = 0.3)
    gg.comp <- gg.comp + geom_point(aes(shape = model), size = 3) 
    gg.comp <- gg.comp + geom_line(aes(linetype = model), size = 1.25)
    gg.comp + ylab("vision")

})

## * equivalence t-test
test_that("lmm - ttest", {
    if(test.practical==FALSE){skip('Not run to save time in the check')}

    ## 2 timepoints
    e052.lmm <- lmm(visual ~ treat.f*week,
                    repetition = ~week|subject, 
                    data = armd.long[armd.long$week %in% c("0","52"),])

    armd.long.imp <- fitted(e052.lmm, type = "outcome", keep.data = TRUE)

    GS <- t.test(armd.long.imp[armd.long.imp$treat.f=="Placebo" & armd.long.imp$week=="52","visual"]-armd.long.imp[armd.long.imp$treat.f=="Placebo" & armd.long.imp$week=="0","visual"],
                 armd.long.imp[armd.long.imp$treat.f=="Active" & armd.long.imp$week=="52","visual"]-armd.long.imp[armd.long.imp$treat.f=="Active" & armd.long.imp$week=="0","visual"])


    expect_equivalent(as.double(GS$estimate),
                      as.double(c(coef(e052.lmm)["week52"], coef(e052.lmm)["week52"] + coef(e052.lmm)["treat.fActive:week52"])),
                      tol = 1e-6)

    ## 4 timepoints
    e052.lmm2 <- lmm(visual ~ treat.f*week,
                    repetition = ~week|subject, 
                    data = armd.long)

    armd.long.imp <- fitted(e052.lmm2, type = "outcome", keep.data = TRUE)

    GS <- t.test(armd.long.imp[armd.long.imp$treat.f=="Placebo" & armd.long.imp$week=="52","visual"]-armd.long.imp[armd.long.imp$treat.f=="Placebo" & armd.long.imp$week=="0","visual"],
                 armd.long.imp[armd.long.imp$treat.f=="Active" & armd.long.imp$week=="52","visual"]-armd.long.imp[armd.long.imp$treat.f=="Active" & armd.long.imp$week=="0","visual"])

    expect_equivalent(as.double(GS$estimate),
                      as.double(c(coef(e052.lmm2)["week52"], coef(e052.lmm2)["week52"] + coef(e052.lmm2)["treat.fActive:week52"])),
                      tol = 1e-6)
})

## * predict function
test_that("lmm - predict", {
    if(test.practical==FALSE){skip('Not run to save time in the check')}

    subject.NNA <- armd.wide$subject[rowSums(is.na(armd.wide))==0]

    armd.longR <- armd.long[armd.long$subject %in% subject.NNA,]
    armd.wideR <- armd.wide[armd.wide$subject %in% subject.NNA,]

    e.UN <- lmm(visual~0+week:treat.f, data = armd.longR,
                repetition = ~week|subject, control = list(optimizer = "FS"))
    e.SUN <- lmm(visual~0+week:treat.f, data = armd.longR,
                 repetition = treat.f~week|subject, control = list(optimizer = "FS"))
    e.ANCOVA <- lm(visual52 ~ visual0 + treat.f, data = armd.wideR)

    ## counterfactual dataset
    armd.longRpl <- armd.longR[armd.longR$week %in% c(0,52),c("subject","week","treat.f","visual")]
    armd.longRpl$treat.f[] <- "Placebo"
    armd.longRpl$visual[armd.longRpl$week != 0] <- NA
    armd.longRa <- armd.longR[armd.longR$week %in% c(0,52),c("subject","week","treat.f","visual")]
    armd.longRa$treat.f[] <- "Active"
    armd.longRa$visual[armd.longRa$week != 0] <- NA

    ## ANCOVA with homogenous variance/correlation
    pred.longRpl <- predict(e.UN, newdata = armd.longRpl, type = "dynamic", se = FALSE, keep.data = TRUE)
    pred.longRa <- predict(e.UN, newdata = armd.longRa, type = "dynamic", se = FALSE, keep.data = TRUE)

    expect_true(all(abs(pred.longRa[pred.longRpl$week!=0,"estimate"]-pred.longRpl[pred.longRpl$week!=0,"estimate"]-coef(e.ANCOVA)["treat.fActive"])<1e-5))

    ## lava::estimate(e.UN, function(p){
    ##     pred.longRpl <- predict(e.UN, p = p, newdata = armd.longRpl, type = "dynamic", se = FALSE)
    ##     pred.longRa <- predict(e.UN, p = p, newdata = armd.longRa, type = "dynamic", se = FALSE)
    ##     pred.longRa[,1] - pred.longRpl[,1]
    ## }, method.numDeriv = "simple", average = TRUE)

    ## ANCOVA with heterogenous variance/correlation
    pred.longRpl <- predict(e.SUN, newdata = armd.longRpl, type = "dynamic", se = FALSE, keep.data = TRUE)
    pred.longRa <- predict(e.SUN, newdata = armd.longRa, type = "dynamic", se = FALSE, keep.data = TRUE)

    expect_equal(-4.3300289,mean(pred.longRa[pred.longRpl$week!=0,"estimate"]-pred.longRpl[pred.longRpl$week!=0,"estimate"]), tol = 1e-5)

    ## lava::estimate(e.SUN, function(p){
    ##     pred.longRpl <- predict(e.SUN, p = p, newdata = armd.longRpl, type = "dynamic", se = FALSE)
    ##     pred.longRa <- predict(e.SUN, p = p, newdata = armd.longRa, type = "dynamic", se = FALSE)
    ##     pred.longRa[,1] - pred.longRpl[,1]
    ## }, method.numDeriv = "simple", average = TRUE)


})

## * baseline constrain
test_that("lmm - baseline constrain", {
    if(test.practical==FALSE){skip('Not run to save time in the check')}

    armd.long$treat <- armd.long$treat.f
    armd.long$treat[armd.long$week == 0] <- "Placebo"
    
    e.lmm <- lmm(visual ~ week:treat,
                 repetition = ~week+treat|subject,
                 structure = UN,
                 control = list(optimizer = "FS"), data = armd.long)
    expect_equal(logLik(e.lmm), -4146.824, tol = 1e-5)
    ## plot(e.lmm)
    plot(e.lmm, color = "treat.f", time.var = 'week')
    plot(e.lmm, type = "correlation")
    sigma(e.lmm)
    
})

## * graphical display
test_that("lmm - autoplot", {
    if(test.practical==FALSE){skip('Not run to save time in the check')}

    armd.long$week.num <- as.numeric(as.character(armd.long$week))
    eLin.lmm <- lmm(visual ~ 0 + week + week.num:treat.f,
                    repetition = ~ week | subject,
                    structure = "UN",
                    data = armd.long)

    expect_true(!is.null(attr(eLin.lmm$design$mean, "variable")))
    autoplot(eLin.lmm) ## was returning an error because could not identify the mean variables necessary for the prediction
})
##----------------------------------------------------------------------
### test-manual-armd.R ends here
