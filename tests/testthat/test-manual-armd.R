### test-manual-armd.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Dec 19 2021 (17:07) 
## Version: 
## Last-Updated: Dec 19 2021 (17:38) 
##           By: Brice Ozenne
##     Update #: 12
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

context("Check lmm on Basic Statistic course")
LMMstar.options(optimizer = "FS", method.numDeriv = "simple", precompute.moments = TRUE)
test.practical <- FALSE

## * Data
## Load data
data(armd.wide, package = "nlmeU")

## and move to long format
library(reshape2)
armd.long <- melt(armd.wide,
                  measure.vars = paste0("visual",c(0,4,12,24,52)),
                  id.var = c("subject","lesion","treat.f","miss.pat"),
                  variable.name = "week",
                  value.name = "visual")

armd.long$week <- factor(armd.long$week, 
                         level = paste0("visual",c(0,4,12,24,52)),
                         labels = c(0,4,12,24,52))

## * summarize
test_that("summarize", {
    if(test.practical==FALSE){skip('Not run to save time in the check')}
    armd.sum <- summarize(visual ~ week * treat.f | subject, ## value ~ group
                          data = armd.long, ## dataset
                          na.rm = TRUE) ## additional argument
    expect_equivalent(attr(armd.sum,"correlation")$visual$Placebo,
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
    levels(e052.lmm)$reference

    ## chunk 21
    c(placebo.0 = as.double(coef(e052.lmm)["(Intercept)"]),
      placebo.52 = sum(coef(e052.lmm)[c("(Intercept)","week52")]),
      active.0 = sum(coef(e052.lmm)[c("(Intercept)","treat.fActive")]),
      active.52 = sum(coef(e052.lmm)))

    ## chunk 22
    dummy.coef(e052.lmm)

    armd.114 <- armd.long[armd.long$subject=="114" & armd.long$week %in% c("0","52"),]
    armd.114

    ## chunk 23
    sigma0 <- coef(e052.lmm, effects="variance")["sigma"]
    sigma52 <- sigma0*coef(e052.lmm, effects="variance")["k.52"]
    rho <- coef(e052.lmm, effects="correlation")["rho(0,52)"]
    alpha0 <- coef(e052.lmm)["(Intercept)"]
    alpha52 <- alpha0 + coef(e052.lmm)["week52"]

    ## chunk 24
    
    expect_equivalent(predict(e052.lmm, newdata = armd.114, type = "dynamic")$estimate,
                      unname(alpha52 + rho * sigma52/sigma0 * (armd.114$visual[1]-alpha0)),
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

    eFlex.lmm <- lmm(visual ~ week*treat.f,
                     repetition = ~ week | subject, structure = "UN",
                     data = armd.long)
    model.tables(eFlex.lmm)

    armd.long.imp <- fitted(eFlex.lmm, impute = TRUE, keep.newdata = TRUE)
    gg <- plot(eFlex.lmm, obs.alpha = 0.1, ci = FALSE, plot = FALSE)$plot
    gg <- gg + geom_point(data = armd.long.imp[armd.long.imp$subject %in% c("114","167"),,drop=FALSE], aes(x = week, y = visual, group = subject, color = treat.f, shape = imputed), size = 4)
    gg <- gg + geom_line(data = armd.long.imp[armd.long.imp$subject %in% c("114","167"),,drop=FALSE], aes(x = week, y = visual, group = subject, color = treat.f))
    gg
    
    armdU <- unique(armd.long[,c("week","week.num","treat.f")])
    armdU

    pLin <- predict(eLin.lmm, newdata = armdU, keep.newdata = TRUE)
    pFlex <- predict(eFlex.lmm, newdata = armdU, keep.newdata = TRUE)
    armdU <- rbind(cbind(pLin, model = "linear"),
                   cbind(pFlex, model = "flexible"))

    gg.comp <- ggplot(armdU, aes(x = week, y = estimate, color = treat.f,
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

    armd.long.imp <- fitted(e052.lmm, impute = TRUE, keep.newdata = TRUE)

    GS <- t.test(armd.long.imp[armd.long.imp$treat.f=="Placebo" & armd.long.imp$week=="52","visual"]-armd.long.imp[armd.long.imp$treat.f=="Placebo" & armd.long.imp$week=="0","visual"],
                 armd.long.imp[armd.long.imp$treat.f=="Active" & armd.long.imp$week=="52","visual"]-armd.long.imp[armd.long.imp$treat.f=="Active" & armd.long.imp$week=="0","visual"])


    expect_equivalent(as.double(GS$estimate),
                      as.double(c(coef(e052.lmm)["week52"], coef(e052.lmm)["week52"] + coef(e052.lmm)["treat.fActive:week52"])),
                      tol = 1e-6)

    ## 4 timepoints
    e052.lmm2 <- lmm(visual ~ treat.f*week,
                    repetition = ~week|subject, 
                    data = armd.long)

    armd.long.imp <- fitted(e052.lmm2, impute = TRUE, keep.newdata = TRUE)

    GS <- t.test(armd.long.imp[armd.long.imp$treat.f=="Placebo" & armd.long.imp$week=="52","visual"]-armd.long.imp[armd.long.imp$treat.f=="Placebo" & armd.long.imp$week=="0","visual"],
                 armd.long.imp[armd.long.imp$treat.f=="Active" & armd.long.imp$week=="52","visual"]-armd.long.imp[armd.long.imp$treat.f=="Active" & armd.long.imp$week=="0","visual"])

    expect_equivalent(as.double(GS$estimate),
                      as.double(c(coef(e052.lmm2)["week52"], coef(e052.lmm2)["week52"] + coef(e052.lmm2)["treat.fActive:week52"])),
                      tol = 1e-6)
})
##----------------------------------------------------------------------
### test-manual-armd.R ends here
