### test-auto-predict.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 22 2024 (10:15) 
## Version: 
## Last-Updated: maj  7 2024 (14:53) 
##           By: Brice Ozenne
##     Update #: 32
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

if(FALSE){
    ## library(pracma) ## trapz
    library(nlmeU)
    library(reshape2)
    library(testthat)

    library(LMMstar)
}

context("Check predict and fitted functions")
LMMstar.options(optimizer = "FS", precompute.moments = TRUE,
                columns.confint = c("estimate","se","df","lower","upper","p.value"))

## * load and reshape data
data("armd.wide", package = "nlmeU")
armd.long <- reshape2::melt(armd.wide,
                            measure.vars = paste0("visual",c(0,4,12,24,52)),
                            id.var = c("subject","lesion","treat.f","miss.pat"),
                            variable.name = "week",
                            value.name = "visual")

armd.long$week <- factor(armd.long$week, 
                         level = paste0("visual",c(0,4,12,24,52)),
                         labels = c(0,4,12,24,52))
armd.long$lesion[is.na(armd.long$lesion)] <- 1
armd.long$week2 <- armd.long$week

armd.longNNA <- armd.long[!is.na(armd.long$lesion) & !is.na(armd.long$visual),]
## sum(is.na(armd.longNNA))
## 0

nd <- armd.long[armd.long$subject %in% c(2,10),]
nd.X <- model.matrix( ~ lesion + week + week:treat.f, data = nd)
nd2 <- armd.long[armd.long$subject %in% 1:3,]
nd2[nd2$subject == 1, "visual"] <- c(59,55,45,40,38)
nd2[nd2$subject == 2, "visual"] <- NA
nd2[nd2$subject == 3, "visual"] <- c(40,40,37,NA,NA)

ndNNA <- armd.longNNA[armd.longNNA$subject %in% c(1,2,10),]
ndNNA.X <- model.matrix( ~ lesion + week + week:treat.f, data = ndNNA)

## * predict - lm
test_that("predict (lmm)", {
    e.lm <- lm(visual ~ lesion + week + week:treat.f,
               data = armd.long)
    e.lm2 <- suppressWarnings(lmm(visual ~ lesion + week + week:treat.f,
                                  data = armd.long))
    
    expect_equivalent(predict(e.lm, newdata = armd.long, type = "terms"),
                      predict(e.lm2, newdata = armd.long, type = "terms"),
                      tol = 1e-6)

    GS <- predict(e.lm, newdata = armd.long, type = "response", se.fit = TRUE)
    test <- predict(e.lm2, newdata = armd.long, se = TRUE)
    expect_equivalent(GS$fit, test$estimate, tol = 1e-6)
    expect_equivalent(GS$se.fit, test$se, tol = 1e-6)
    expect_equivalent(GS$df, unique(round(test$df)))
})


## * predict/fitted - lmm
e.lmm <- lmm(visual ~ lesion + week + week:treat.f,
             repetition = ~ week | subject, structure = "UN",
             data = armd.long)
e.lmmNNA <- lmm(visual ~ lesion + week + week:treat.f,
             repetition = ~ week | subject, structure = "UN",
             data = armd.longNNA)

test_that("predict/fitted (lmm)", {

    ##-- static
    testthat::expect_visible(predict(e.lmm, newdata = nd)) ## check it outputs something
    expect_equal(predict(e.lmm, newdata = nd, type = "static", format = "long", se = FALSE),
                 as.double(nd.X %*% coef(e.lmm)), tol = 1e-6)
    expect_equivalent(predict(e.lmm, newdata = nd, type = "static", format = "long", se = TRUE)[,c("estimate","se")],
                      data.frame(estimate = nd.X %*% coef(e.lmm),
                                 se = sqrt(diag(nd.X %*% vcov(e.lmm) %*% t(nd.X)))),
                      tol = 1e-6)
    expect_visible(predict(e.lmm, newdata = nd, type = "static", format = "long", se = TRUE, keep.data = TRUE))
    expect_visible(predict(e.lmm, newdata = nd, type = "static", format = "wide"))
    GS <- data.frame("subject" = c("2", "10"), 
                     "lesion" = c(1, 1), 
                     "treat.f" = as.factor(c("Active", "Placebo")), 
                     "estimate_0" = c(57.47046, 58.19500), 
                     "estimate_4" = c(53.99483, 56.91638), 
                     "estimate_12" = c(51.62498, 55.84074), 
                     "estimate_24" = c(48.38935, 52.16911), 
                     "estimate_52" = c(41.27514, 46.88251))
    expect_equivalent(predict(e.lmm, newdata = nd, type = "static", format = "wide", keep.data = TRUE), GS, tol = 1e-6)

    ##-- time ordering
    test1 <- predict(e.lmm, newdata = nd[NROW(nd):1,], se = c(TRUE,TRUE))
    test2 <- predict(e.lmm, newdata = nd, se = c(TRUE,TRUE))
    expect_equivalent(test1, test2[NROW(nd):1,], tol = 1e-6)
    ## sqrt(diag(attr(predict(e.lmm, newdata = nd, se = "total2", simplify = FALSE),"vcov")))

    ##-- static NA
    expect_equivalent(predict(e.lmmNNA, newdata = nd, type = "static", format = "wide", keep.data = TRUE), GS, tol = 1e-6)
    GSNNA <- GS
    GSNNA[2,"estimate_52"] <- NA
    expect_equivalent(predict(e.lmmNNA, newdata = ndNNA[ndNNA$subject %in% GSNNA$subject,,drop=FALSE], type = "static", format = "wide", keep.data = TRUE),
                      GSNNA, tol = 1e-6)

    ##-- unique static
    Ufit <- fitted(e.lmm, newdata = "unique")
    
    ##-- dynamic
    predL <- predict(e.lmm, newdata = nd2, type = "dynamic", keep.data = TRUE)
    expect_equal(is.na(nd2$visual),!is.na(predL$estimate))
    expect_equal(c(NA, 1.5957344, NA, NA, 1.66930219, NA, NA, 1.78849944, NA, NA, 1.92606744, 1.07271763, NA, 2.00302196, 1.50857931),predL$se, tol = 1e-6)
    expect_equal(c(NA, 243.9075, NA, NA, 257.7478, NA, NA, 272.432, NA, NA, 275.4818, 174.074, NA, 263.3198, 198.7918), predL$df, tol = 1e-1)
    ## similar to not using a transformation
    ## predL2 <- predict(e.lmm, newdata = nd2, type = "dynamic", transform.sigma = "none", transform.k = "none", transform.rho = "none",
    ##                   keep.data = TRUE)
    ## expect_equal(predL2$estimate,predL$estimate, tol = 1e-6)
    ## expect_equal(predL2$se,predL$se, tol = 1e-3)
    ## expect_equal(predL2$df,predL$df, tol = 1e-1)

    predLse2 <- predict(e.lmm, newdata = nd2, type = "dynamic", se = c(TRUE,TRUE))
    predW <- predict(e.lmm, newdata = nd2, type = "dynamic", format = "wide", keep.data = TRUE)
    expect_equal(colnames(predW), c("subject", "lesion", "treat.f", "estimate_0", "estimate_4", "estimate_12", "estimate_24", "estimate_52"))
    fitL <- fitted(e.lmm, newdata = nd2, type = "outcome", format = "long", keep.data = TRUE, export.vcov = TRUE)
    ## GS <- estimate(e.lmm, function(p){ ## p <- NULL
    ##     predict(e.lmm, p = p, newdata = nd2, type = "dynamic")
    ## })
    expect_equal(fitL$se, c(NA, 1.5957344, NA, NA, 1.66930219, NA, NA, 1.78849944, NA, NA, 1.92606744, 1.07271669, NA, 2.00302196, 1.50857614), tol = 1e-6)
    expect_equal(fitL$df, c(NA, 243.784, NA, NA, 257.264, NA, NA, 271.656, NA, NA, 275.089, 137.352, NA, 262.783, 155.409), tol = 1e-1)

    fitW <- fitted(e.lmm, newdata = nd2, type = "outcome", format = "wide")
    
    ##-- impute
    imputeW <- fitted(e.lmm, newdata = nd2, type = "impute", format = "wide")
    imputeL <- fitted(e.lmm, newdata = nd2, type = "impute", format = "long")
   
    ##-- change
    changeL <- fitted(e.lmm, newdata = nd2, type = "change", format = "long", keep.data = TRUE)
    expect_equal(fitL[fitL$subject==3,"se"], changeL[changeL$subject==3,"se"])
    expect_equal(fitL[fitL$subject==3,"df"], changeL[changeL$subject==3,"df"])

    M.change <- matrix(c(-1,rep(0,4)), byrow = TRUE, nrow = 5, ncol = 5) + diag(1, 5, 5)
    GS.est <- M.change %*% fitL[fitL$subject==2,"visual"]
    GS.se <- sqrt(diag(M.change %*% attr(fitL,"vcov")[fitL$subject==2,fitL$subject==2] %*% t(M.change)))

    expect_equal(changeL[changeL$subject==2,"visual"], GS.est[,1], tol = 1e-6)
    expect_equal(changeL[changeL$subject==2,"se"], c(NA,GS.se[-1]), tol = 1e-6)
   
    changeW <- fitted(e.lmm, newdata = nd2, type = "change", format = "wide")
    expect_equal(changeW$visual_0, rep(0,3), tol = 1e-6)
    expect_equal(changeW$visual_4, fitW$visual_4-fitW$visual_0, tol = 1e-6)
    expect_equal(changeW$visual_12, fitW$visual_12-fitW$visual_0, tol = 1e-6)
   
    ## auc
    aucL <- fitted(e.lmm, newdata = nd2, type = "auc", format = "long", keep.data = TRUE)
    aucbL <- fitted(e.lmm, newdata = nd2, type = "auc-b", format = "long", keep.data = TRUE)
})



##----------------------------------------------------------------------
### test-auto-predict.R ends here
