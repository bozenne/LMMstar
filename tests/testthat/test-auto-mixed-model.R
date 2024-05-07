### test-auto-mixed-model.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 14 2021 (16:46) 
## Version: 
## Last-Updated: maj  7 2024 (12:20) 
##           By: Brice Ozenne
##     Update #: 217
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
    library(numDeriv)
    library(lava)
    library(multcomp)
    library(lme4) ## removed due to error in appveyer and Github action (change of Matrix package version)
    library(lmerTest)

    library(LMMstar)
}

context("Check lmm on mixed model parametrized with random effects")
LMMstar.options(optimizer = "FS", method.numDeriv = "simple", precompute.moments = TRUE, # "Richardson"
                columns.confint = c("estimate","se","df","lower","upper","p.value"))

## * Random intercept model
data(Orthodont,package="nlme")

test_that("Random intercept model",{
    
    ## ** fit
    ## eRI.lmer <- lmer(distance ~ age + (1|Subject), data=Orthodont)
    ## eRI.lmm <- lmm(distance ~ age + (1|Subject), data=Orthodont,
    ##                control = list(init = "lmer"))
    eRI2.lmm <- lmm(distance ~ age + (1|Subject), data=Orthodont)
    
    xx <- capture.output(summary(eRI2.lmm))

    ## ** iteration
    ## expect_equal(eRI.lmm$opt$n.iter,0)
    expect_equal(eRI2.lmm$opt$n.iter,4)

    ## ** likelihood
    ## expect_equal(as.double(logLik(eRI.lmer)), as.double(logLik(eRI.lmm)), tol = 1e-6)
    ## expect_equal(as.double(logLik(eRI.lmer)), as.double(logLik(eRI2.lmm)), tol = 1e-6)
    ## expect_equal(-223.5012578, as.double(logLik(eRI.lmm)), tol = 1e-6)
    expect_equal(-223.5012578, as.double(logLik(eRI2.lmm)), tol = 1e-6)

    ## ** random effects (conditional mean)
    ## u.GS <- as.data.frame(ranef(eRI.lmer))    
    u.GS <- data.frame("grpvar" = c("Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject", "Subject"), 
                       "term" = c("(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)", "(Intercept)"), 
                       "grp" = c("M16", "M05", "M02", "M11", "M07", "M08", "M03", "M12", "M13", "M14", "M09", "M15", "M06", "M04", "M01", "M10", "F10", "F09", "F06", "F01", "F05", "F07", "F02", "F08", "F03", "F04", "F11"), 
                       "condval" = c(-0.9179756, -0.9179756, -0.5815230, -0.3572213, -0.2450704, -0.1329195,  0.2035330,  0.2035330,  0.2035330,  0.7642874,  0.9885891,  1.6614942,  2.1100977,  2.3343994,  3.3437571,  4.9138692, -4.9554066, -2.6002385, -2.6002385, -2.3759368, -1.2544282, -0.9179756, -0.9179756, -0.5815230, -0.2450704,  0.7642874,  2.1100977), 
                       "condsd" = c(0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092, 0.6780092))
    u.test <- ranef(eRI2.lmm, effects = "mean", format = "long", simplify = FALSE)
    expect_equal(as.double(u.GS$condval[match(u.test$Subject,u.GS$grp)]), as.double(u.test$estimate), tol = 1e-6)

    ## ** random effects (conditional variance)
    ## tau.GS <- as.data.frame(VarCorr(eRI.lmer))
    tau.GS <- data.frame("grp" = c("Subject", "Residual"), 
                         "var1" = c("(Intercept)", "NA"), 
                         "var2" = c("NA", "NA"), 
                         "vcov" = c(4.472056, 2.049456), 
                         "sdcor" = c(2.114724, 1.431592))
    tau.test <- ranef(eRI2.lmm, effects = "variance", format = "long")
    expect_equal(as.double(tau.GS[,"vcov"]), as.double(tau.test[-1]), tol = 1e-6)

})

## * Stratified random intercept model
Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")

test_that("Stratified random intercept model",{

    ## ** fit
    eRI.mlmm <- mlmm(distance ~ 1 + (1|Subject), data = Orthodont, by = "nsex", trace = FALSE)
    eSRI.lmm <- lmm(distance ~ nsex + (1|Subject), repetition = nsex~1|Subject, data = Orthodont)
    
    ## ** iteration
    expect_equal(eSRI.lmm$opt$n.iter,7)

    ## ** likelihood
    expect_equal(as.double(sum(unlist(logLik(eRI.mlmm)))), as.double(logLik(eSRI.lmm)), tol = 1e-6)

    ## ** random effects (conditional mean)
    u.GS <- do.call(rbind,ranef(eRI.mlmm, simplify = FALSE))
    u.test <- ranef(eSRI.lmm, effects = "mean", format = "long", simplify = FALSE)
    expect_equal(as.double(u.GS[match(u.test$Subject,u.GS$Subject),"estimate"]), as.double(u.test$estimate), tol = 1e-6)

    ## ** random effects (conditional variance)
    tau.GS <- as.data.frame(do.call(rbind,ranef(eRI.mlmm, effects = "variance")))
    tau.test <- ranef(eSRI.lmm, effects = "variance", format = "wide")
    expect_equal(as.double(tau.GS$total), as.double(tau.test[tau.test$type=="total",c("absolute.0","absolute.1")]), tol = 1e-5)
    expect_equal(as.double(tau.GS$Subject), as.double(tau.test[tau.test$type=="Subject",c("absolute.0","absolute.1")]), tol = 1e-5)
    expect_equal(as.double(tau.GS$residual), as.double(tau.test[tau.test$type=="residual",c("absolute.0","absolute.1")]), tol = 1e-5)
})

## * Crossed random intercept model (2 terms)
data(Penicillin, package = "lme4")
Penicillin$id <- 1

test_that("Crossed random intercept model (2 terms)",{

    ## ** fit
    ## eCRI2.lmer <- lmer(diameter ~ (1|plate) + (1|sample), data = Penicillin)
    eCRI2.lmm0 <- lmm(diameter ~ 1, repetition = ~1|id,
                      structure = CS(list(~1,~plate+sample), type = "ho", group = 1:2),
                      data = Penicillin, df = FALSE)
    eCRI2.lmm <- lmm(diameter ~ (1|plate) + (1|sample), data = Penicillin, df = FALSE)
    ## eCRI2.lmm <- lmm(diameter ~ (1|plate) + (1|sample), data = Penicillin, df = FALSE,
    ##                  control = list(init = "lmer"))

    ## ** iteration
    expect_equal(eCRI2.lmm0$opt$n.iter,7)
    expect_equal(eCRI2.lmm$opt$n.iter,7)
    ## expect_equal(eCRI2.lmm$opt$n.iter,2)

    ## ** likelihood
    ## expect_equal(as.double(logLik(eCRI2.lmer)), as.double(logLik(eCRI2.lmm0)), tol = 1e-6)
    ## expect_equal(as.double(logLik(eCRI2.lmer)), as.double(logLik(eCRI2.lmm)), tol = 1e-6)
    expect_equal(c(-165.4302945), as.double(logLik(eCRI2.lmm0)), tol = 1e-6)
    expect_equal(c(-165.4302945), as.double(logLik(eCRI2.lmm)), tol = 1e-6)

    ## ** random effects (conditional mean)
    ## GS <- do.call(rbind,ranef(eCRI2.lmer))
    GS <- data.frame("(Intercept)" = c(0.80454691, 0.80454691, 0.18167188, 0.33739064, 0.02595313, -0.44120314, -1.37551568, 0.80454691, -0.75264065, -0.75264065, 0.96026566, 0.49310939, 1.42742193, 0.49310939, 0.96026566, 0.02595313, -0.28548439, -0.28548439, -1.37551568, 0.96026566, -0.90835941, -0.28548439, -0.5969219, -1.21979692, 2.1870584, -1.01047635, 1.93789985, -0.09689499, -0.01384214, -3.00374477))
    test <- ranef(eCRI2.lmm)
    expect_equal(as.double(GS[,1]), as.double(test$estimate), tol = 1e-6)

    ## ** random effects (conditional variance)
    ## GS <- as.data.frame(VarCorr(eCRI2.lmer))
    GS <- data.frame("grp" = c("plate", "sample", "Residual"), 
                     "var1" = c("(Intercept)", "(Intercept)", "NA"), 
                     "var2" = c("NA", "NA", "NA"), 
                     "vcov" = c(0.7169051, 3.7311318, 0.3024150), 
                     "sdcor" = c(0.8467025, 1.9316138, 0.5499227))
    test <- ranef(eCRI2.lmm, effects = "variance", simplify = FALSE, format = "wide")
    expect_equal(as.double(test$absolute.1[-1]), as.double(GS$vcov), tol = 1e-2)
})

## * Crossed random intercept model (3 terms)
set.seed(10)
Sigma.CRI3 <- matrix(0, nrow = 6, ncol = 6)
Sigma.CRI3[1:3,1:3] <- Sigma.CRI3[4:6,4:6] <- 0.8^2
diag(Sigma.CRI3[1:3,4:6]) <- diag(Sigma.CRI3[4:6,1:3]) <- 0.5^2
Sigma.CRI3[1,6] <- Sigma.CRI3[2,4] <- Sigma.CRI3[3,5] <- Sigma.CRI3[6,1] <- Sigma.CRI3[4,2] <- Sigma.CRI3[5,3] <- 0.3^2
diag(Sigma.CRI3) <- 1

n <- 25
dfL.CRI3 <- reshape2::melt(data.frame(sample = 1:n,
                                      mvtnorm::rmvnorm(n, mean = 1:6, sigma = Sigma.CRI3)
                                      ),
                           id.vars = "sample")
dfL.CRI3$patient <- paste(dfL.CRI3$sample,dfL.CRI3$variable %in% paste0("X",4:6)+1, sep = ".")
dfL.CRI3$day <- paste(dfL.CRI3$sample, sapply(as.character(dfL.CRI3$variable), switch, "X1" = 1, "X2" = 2, "X3" = 3, "X4" = 1, "X5" = 2, "X6" = 3),sep=".")
dfL.CRI3$batch <- paste(dfL.CRI3$sample, sapply(as.character(dfL.CRI3$variable), switch, "X1" = 1, "X2" = 2, "X3" = 3, "X4" = 2, "X5" = 3, "X6" = 1),sep=".")
dfL.CRI3 <- dfL.CRI3[order(dfL.CRI3$sample),]

## head(dfL.CRI3,10)

test_that("Crossed random intercept model (3 terms)",{

    ## ** fit
    ## eCRI3.lmer <- lmer(value ~ 0 + variable + (1|batch) + (1|day) + (1|patient), data = dfL.CRI3)
    eCRI3.lmm0 <- lmm(value ~ 0 + variable + (1|batch) + (1|day) + (1|patient), data = dfL.CRI3, df = FALSE)
    ## eCRI3.lmm <- lmm(value ~ 0 + variable + (1|batch) + (1|day) + (1|patient), data = dfL.CRI3, df = FALSE, control = list(init = "lmer"))

    ## ** iteration
    expect_equal(eCRI3.lmm0$opt$n.iter,5)
    ## expect_equal(eCRI3.lmm$opt$n.iter,0)

    ## ** likelihood
    ## expect_equal(as.double(logLik(eCRI3.lmer)), as.double(logLik(eCRI3.lmm0)), tol = 1e-6)
    ## expect_equal(as.double(logLik(eCRI3.lmer)), as.double(logLik(eCRI3.lmm)), tol = 1e-6)
    expect_equal(-163.0799, as.double(logLik(eCRI3.lmm0)), tol = 1e-6)
    ## expect_equal(-163.0799, as.double(logLik(eCRI3.lmm)), tol = 1e-6)

    ## ** random effects (conditional mean)
    ## GS <- ranef(eCRI3.lmer)
    GS <- list(patient = structure(list('(Intercept)'=c(-0.460885985821923,0.0446045409822634,0.391967982341095,0.0734263331932129,-0.914365324774579,1.17483234137768,0.597005548078684,-1.02504470878643,-0.996043272958032,-0.312008691353995,0.683867800186977,0.524841299947234,0.771825425270793,0.415425088200219,-0.0201590737708554,0.146971021343681,-0.0357809361943112,0.392505426139319,-0.202773088310432,-0.18897669946572,0.606152559155231,0.146060919471704,-1.20653563240842,0.726657636205626,1.1057680116466,1.28998820299617,0.0924540922056612,-0.659462553646222,-0.306790820996071,0.383589959287912,-1.1713880737607,0.122551701798961,0.207624507789385,-0.754553663646655,1.08278389788988,0.660880854015443,0.887389771181412,-0.482277909200269,0.611891574836963,-2.21358629586407,-0.802244896106532,-0.531879005155838,-0.189441797511981,-1.16194222515625,-0.540066671190878,-0.137030002202902,-0.205033299493541,-0.0399750206488556,0.012777702715505,1.40440145016628)),
                                   class = "data.frame",
                                   row.names = c("1.1","1.2","10.1","10.2","11.1","11.2","12.1","12.2","13.1","13.2","14.1","14.2","15.1","15.2","16.1","16.2","17.1","17.2","18.1","18.2","19.1","19.2","2.1","2.2","20.1","20.2","21.1","21.2","22.1","22.2","23.1","23.2","24.1","24.2","25.1","25.2","3.1","3.2","4.1","4.2","5.1","5.2","6.1","6.2","7.1","7.2","8.1","8.2","9.1","9.2")),
               day = structure(list('(Intercept)'=c(-0.000962135271624641,0.0437364730438857,-0.124207986028689,0.29756998256484,-0.0954684600270453,-0.111060335231542,-0.175867662960196,0.0537033624187523,0.173117278809866,0.0839803319420781,-0.0497419893366753,-0.117972054220817,0.308306743301917,-0.381194735032493,-0.182995245159664,-0.425019792841041,0.346822834276196,0.314646583992148,0.672147451336675,-0.58034438952156,0.14044880024783,0.339791972046702,-0.0384985468676835,-0.276486267278109,-0.0682785702976616,-0.139764673351714,0.277826264244546,0.0589944628582359,0.312144888108882,-0.447774076866517,0.566222844400408,0.00279699685586821,-0.421870626244189,-0.230709846370439,0.212128092295435,-0.0752927542250942,0.249315663134291,-0.0641192412707053,0.283465276393277,-0.242455352946775,0.22463015220863,-0.0930939098293566,-0.134673058082534,0.0437706541721264,0.10592597532245,-0.552918659090568,0.0199300213047896,0.3278135735314,-0.0439620646212562,-0.206019828999268,0.142990729225909,0.00625352166473978,0.323429426141537,0.0114155658170881,0.00272667060129977,-0.15407265404841,0.230594619873331,0.00966654058183359,0.121084048960365,-0.444076696538249,-0.280434886098894,-0.0203766552434973,0.0398280700571376,-0.128721827854857,-0.359002606138923,0.223364509895525,0.188188303458032,-0.162658869825859,-0.157984177941005,-0.102445307081827,0.244672605814845,-0.190156221844305,-0.396715324373414,0.302413094062433,0.37153310399687)),
                               class = "data.frame",
                               row.names = c("1.1","1.2","1.3","10.1","10.2","10.3","11.1","11.2","11.3","12.1","12.2","12.3","13.1","13.2","13.3","14.1","14.2","14.3","15.1","15.2","15.3","16.1","16.2","16.3","17.1","17.2","17.3","18.1","18.2","18.3","19.1","19.2","19.3","2.1","2.2","2.3","20.1","20.2","20.3","21.1","21.2","21.3","22.1","22.2","22.3","23.1","23.2","23.3","24.1","24.2","24.3","25.1","25.2","25.3","3.1","3.2","3.3","4.1","4.2","4.3","5.1","5.2","5.3","6.1","6.2","6.3","7.1","7.2","7.3","8.1","8.2","8.3","9.1","9.2","9.3")), 
               batch = structure(list('(Intercept)'=c(0.166070218345111,-0.114003724399762,-0.0971460442350244,0.137348274951939,0.134441636696807,-0.221391877550108,0.0416616428741636,-0.0449575131730436,0.0315021142218759,-0.0108872110178352,-0.257997119800447,0.222531525257832,-0.0831047450321626,0.226595766190439,-0.285141327790195,-0.0645218667183461,0.102678869600867,0.0927353617437159,-0.0328857910029793,0.256233210864148,-0.0947788278859101,-0.234215587499144,0.0951757740271408,0.152772411066265,-0.162799707419394,0.300211204012577,-0.098781430844689,-0.292291950027528,0.173805601185707,0.0760633571905189,0.107652700091309,0.227726872706668,-0.253921593446057,0.0531335806324594,-0.0791601969138935,-0.0259398706492175,0.138067233713131,-0.0775099914359328,0.198881688039271,-0.0837277522959336,-0.130505247635422,0.152831064368202,0.196133519820784,-0.366628565839104,0.178811704518833,-0.0397511455850197,-0.156585815880064,0.0827573807232534,0.422970165435738,0.00212943841559413,-0.484327130545092,0.00941322924340147,0.0100887075925886,0.169321331282844,-0.0519996394006836,0.131189215693038,-0.0353195917386138,0.0289044327975289,-0.179414611636477,-0.0229390070911693,-0.0565311747831045,-0.0914610394579688,0.00351855046107388,-0.35655073897538,0.134101624335313,0.0761063343751364,0.0466356702335856,-0.344990521421639,0.225031348741327,-0.0690346634128251,0.110475659487739,-0.0679732016924462,0.19031130503112,-0.0616392332148874,0.0247957314782996)),
                                 class = "data.frame",
                                 row.names = c("1.1","1.2","1.3","10.1","10.2","10.3","11.1","11.2","11.3","12.1","12.2","12.3","13.1","13.2","13.3","14.1","14.2","14.3","15.1","15.2","15.3","16.1","16.2","16.3","17.1","17.2","17.3","18.1","18.2","18.3","19.1","19.2","19.3","2.1","2.2","2.3","20.1","20.2","20.3","21.1","21.2","21.3","22.1","22.2","22.3","23.1","23.2","23.3","24.1","24.2","24.3","25.1","25.2","25.3","3.1","3.2","3.3","4.1","4.2","4.3","5.1","5.2","5.3","6.1","6.2","6.3","7.1","7.2","7.3","8.1","8.2","8.3","9.1","9.2","9.3"))
               )
    test <- ranef(eCRI3.lmm0)
    expect_equal(as.double(test[test$variable=="patient","estimate"]), as.double(GS$patient[test[test$variable=="patient","level"],1]), tol = 1e-6)
    expect_equal(as.double(test[test$variable=="day","estimate"]), as.double(GS$day[test[test$variable=="day","level"],1]), tol = 1e-6)
    expect_equal(as.double(test[test$variable=="batch","estimate"]), as.double(GS$batch[test[test$variable=="batch","level"],1]), tol = 1e-6)

    ## ** random effects (conditional variance)
    ## GS <- as.data.frame(VarCorr(eCRI3.lmer))
    GS <- data.frame("grp" = c("batch", "day", "patient", "Residual"), 
                     "var1" = c("(Intercept)", "(Intercept)", "(Intercept)", "NA"), 
                     "var2" = c("NA", "NA", "NA", "NA"), 
                     "vcov" = c(0.06942091, 0.12540494, 0.64105869, 0.06404538), 
                     "sdcor" = c(0.2634785, 0.3541256, 0.8006614, 0.2530719))
    test <- ranef(eCRI3.lmm0, effects = "variance", format = "wide", simplify = FALSE)
    expect_equal(as.double(test[-1,"absolute.1"]), as.double(GS[,"vcov"]), tol = 1e-3)
})

## * Nested random intercept model (2 levels)
## https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified
## df.school <- read.table("http://bayes.acs.unt.edu:8083/BayesContent/class/Jon/R_SC/Module9/lmm.data.txt",
##                  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
## eNRI2.lmer <- lmer(extro ~ open + agree + social + (1 | school/class), data = dt)
## dt.red <- as.data.table(dt)[,.SD[1:10,.(id,extro=round(extro,2))],by = c("class","school")]
## setkeyv(dt.red, c("school","class"))
## dt.red$id <- 1:NROW(dt.red)

df.red <- data.frame("class" = rep(unlist(lapply(c("a","b","c","d"),rep,10)), 6), 
                     "school" = unlist(lapply(c("I","II","III","IV","V","VI"),rep,40)), 
                     "id" = 1:240, 
                     "extro" = c(41.10, 43.24, 42.15, 43.72, 40.59, 41.38, 35.43, 39.46, 41.06, 38.96, 47.24, 45.77, 44.32, 47.34, 45.06, 47.01, 46.07, 44.41, 45.30, 45.20, 47.65, 47.76, 48.59, 47.46, 48.63, 48.31, 49.28, 47.81, 49.01, 48.00, 50.97, 49.28, 49.66, 50.32, 49.33, 49.87, 50.87, 51.43, 50.95, 50.87, 52.02, 51.67, 52.62, 51.96, 51.88, 52.56, 51.54, 52.08, 52.76, 52.69, 53.69, 54.10, 54.07, 53.70, 53.03, 53.88, 53.80, 53.50, 53.42, 53.17, 55.17, 54.58, 55.06, 54.40, 54.86, 54.45, 55.27, 54.60, 55.33, 54.81, 55.63, 55.65, 55.96, 56.01, 55.58, 56.01, 56.08, 55.47, 56.07, 55.56, 56.68, 57.30, 56.84, 57.24, 56.57, 56.94, 56.81, 57.11, 57.04, 56.35, 57.75, 57.84, 57.66, 58.03, 58.22, 57.50, 57.54, 57.42, 58.11, 58.08, 59.03, 58.41, 59.04, 59.09, 58.94, 58.39, 58.30, 58.80, 58.95, 58.35, 60.15, 59.45, 59.43, 59.98, 59.40, 60.05, 60.00, 59.90, 59.93, 59.80, 60.65, 60.84, 60.66, 60.31, 60.98, 60.75, 60.60, 60.42, 60.74, 60.17, 62.29, 61.84, 61.57, 61.17, 62.12, 61.59, 61.42, 61.36, 61.67, 62.12, 62.97, 63.15, 62.82, 62.84, 63.29, 62.68, 62.45, 63.21, 62.80, 62.75, 63.69, 64.25, 64.18, 64.14, 63.72, 63.42, 64.18, 63.68, 63.48, 64.04, 64.62, 65.13, 65.40, 65.23, 64.56, 64.63, 64.81, 64.31, 64.76, 64.44, 65.92, 66.49, 66.28, 65.68, 66.03, 66.46, 65.62, 66.01, 65.61, 66.10, 66.57, 66.98, 66.55, 66.96, 66.64, 67.24, 67.23, 67.34, 66.94, 67.30, 67.94, 68.01, 68.73, 68.61, 68.89, 68.91, 68.84, 68.06, 68.13, 68.93, 69.48, 70.12, 70.17, 70.67, 70.60, 69.62, 70.16, 70.65, 70.34, 69.87, 71.79, 71.36, 72.45, 71.07, 70.87, 70.94, 71.25, 72.48, 72.04, 71.67, 74.35, 74.70, 73.14, 74.44, 73.98, 75.03, 74.21, 76.51, 75.94, 73.42, 79.74, 78.15, 78.19, 83.34, 80.11, 79.73, 80.64, 77.07, 78.02, 82.25))

test_that("Nested random intercept model (2 levels)",{

    ## ** fit
    ## eNRI2.lmer <- lmer(extro ~ (1|school/class), data = df.red)
    eNRI2.lmm0 <- lmm(extro ~ (1|school/class), data = df.red, df = FALSE)
    ## eNRI2.lmm <- lmm(extro ~ (1|school/class), data = df.red, df = FALSE, control = list(init = "lmer"))

        ## ** iteration
    expect_equal(eNRI2.lmm0$opt$n.iter,9)
    ## expect_equal(eNRI2.lmm$opt$n.iter,3)

    ## ** likelihood
    ## expect_equal(as.double(logLik(eNRI2.lmer)), as.double(logLik(eNRI2.lmm0)), tol = 1e-6)
    ## expect_equal(as.double(logLik(eNRI2.lmer)), as.double(logLik(eNRI2.lmm)), tol = 1e-6)
    expect_equal(-351.01958598, as.double(logLik(eNRI2.lmm0)), tol = 1e-6)
    ## expect_equal(-351.01958598, as.double(logLik(eNRI2.lmm)), tol = 1e-6)

    ## ** random effects (conditional mean)
    ## GS <- ranef(eNRI2.lmer)
    GS <- list(`class:school` = structure(list(`(Intercept)` = c(-5.78223961693662,-2.04010426217362,-1.44752110028082,-1.61025488289359,-1.64692110248124,-3.54928781003974,-0.76308025475413,-0.594729118965627,-0.528548015401662,-0.51680578621429,-0.426580360856386,-2.13761826002042,1.69346268411378,0.611732848993305,0.378528977763854,0.653967906603711,0.520150271354621,0.816577848868358,3.78023543888253,1.55251543534749,1.44818589101586,1.62746474517042,2.03690196484449,5.92396655792706)),
                                          class = "data.frame",
                                          row.names = c("a:I","a:II","a:III","a:IV","a:V","a:VI","b:I","b:II","b:III","b:IV","b:V","b:VI","c:I","c:II","c:III","c:IV","c:V","c:VI","d:I","d:II","d:III","d:IV","d:V","d:VI")),
               school = structure(list(`(Intercept)` = c(-13.7183795009558,-6.02420112625015,-1.91196029908475,1.976194907378,6.19018139383447,13.4881646229664)),
                                  class = "data.frame",
                                  row.names = c("I", "II", "III", "IV", "V", "VI")
                                  )
               )
    test <- ranef(eNRI2.lmm0)
    test$label <- ifelse(!is.na(test$class),paste0(test$class,":",test$school),test$school)
    expect_equal(as.double(test[is.na(test$class),"estimate"]), as.double(GS$school[,1]), tol = 1e-4)
    expect_equal(as.double(test[match(rownames(GS$class),test$label),"estimate"]), as.double(GS$class[,1]), tol = 1e-4)

    ## ** random effects (conditional variance)
    ## GS <- as.data.frame(VarCorr(eNRI2.lmer))
    GS <- data.frame("grp" = c("class:school", "school", "Residual"), 
                     "var1" = c("(Intercept)", "(Intercept)", "NA"), 
                     "var2" = c("NA", "NA", "NA"), 
                     "vcov" = c( 7.2056645, 92.2434064,  0.6293901), 
                     "sdcor" = c(2.6843369, 9.6043431, 0.7933411))
    test <- ranef(eNRI2.lmm0, effects = "variance", format = "wide", simplify = FALSE)
    expect_equal(as.double(test[-1,"absolute.1"]), as.double(GS[c(2,1,3),"vcov"]), tol = 1e-3)
})


## * Nested random intercept model (3 levels)
Sigma.NRI3 <- matrix(0, nrow = 8, ncol = 8)
Sigma.NRI3[5:8,1:4] <- Sigma.NRI3[1:4,5:8] <- 0.25^2
Sigma.NRI3[1:2,3:4] <- Sigma.NRI3[3:4,1:2] <- Sigma.NRI3[5:6,7:8] <- Sigma.NRI3[7:8,5:6] <- 0.5^2
Sigma.NRI3[1:2,1:2] <- Sigma.NRI3[3:4,3:4] <- Sigma.NRI3[5:6,5:6] <- Sigma.NRI3[7:8,7:8] <- 0.8^2
diag(Sigma.NRI3) <- 1

## c(sqrt(0.25^2), sqrt(0.5^2-0.25^2), sqrt(0.8^2-0.5^2-0.25^2))
n <- 1000
set.seed(10)
dfL.NRI3 <- reshape2::melt(data.frame(patient = 1:n,
                                      mvtnorm::rmvnorm(n, mean = 1:8, sigma = Sigma.NRI3)
                                      ),
                           id.vars = "patient")
dfL.NRI3$day <- dfL.NRI3$variable %in% paste0("X",5:8)+1
dfL.NRI3$session <- paste(dfL.NRI3$day,dfL.NRI3$variable %in% c(paste0("X",c(3:4,7:8)))+1,sep=".")
dfL.NRI3 <- dfL.NRI3[order(dfL.NRI3$patient),]

## head(dfL.NRI3,10)

test_that("Nested random intercept model (2 levels)",{

    ## ** fit
    ## eNRI3.lmer <- lmer(value ~ session + (1|patient/day/session), data = dfL.NRI3)
    eNRI3.lmm0 <- lmm(value ~ session + (1|patient/day/session), data = dfL.NRI3, df = FALSE)
    ## eNRI3.lmm <- lmm(value ~ session + (1|patient/day/session), data = dfL.NRI3, df = FALSE, control = list(init = "lmer"))
    
    ## ** iteration
    expect_equal(eNRI3.lmm0$opt$n.iter,4)
    ## expect_true(eNRI3.lmm$opt$n.iter<=2)

    ## ** likelihood
    ## expect_equal(as.double(logLik(eNRI3.lmer)), as.double(logLik(eNRI3.lmm0)), tol = 1e-6)
    ## expect_equal(as.double(logLik(eNRI3.lmer)), as.double(logLik(eNRI3.lmm)), tol = 1e-6)
    expect_equal(-12046.61, as.double(logLik(eNRI3.lmm0)), tol = 1e-6)
    ## expect_equal(-12046.61, as.double(logLik(eNRI3.lmm)), tol = 1e-6)

    ## ** random effects (conditional mean)
    ## slow!!
    ## GS <- ranef(eNRI3.lmer)
    ## test <- ranef(eNRI3.lmm0)
    ## test$variable2 <- sapply(test$variable,paste, collapse=":")
    ## expect_equal(as.double(test[test$variable2=="patient","estimate"]), as.double(GS$patient[,1]), tol = 1e-4)
    ## expect_equal(as.double(test[test$variable2=="patient:day","estimate"]), as.double(GS$day[,1]), tol = 1e-4)
    ## expect_equal(as.double(test[test$variable2=="patient:day:session","estimate"]), as.double(GS$session[,1]), tol = 1e-4)

    ## ** random effects (conditional variance)
    ## GS <- as.data.frame(VarCorr(eNRI3.lmer))
    GS <- data.frame("grp" = c("session:(day:patient)", "day:patient", "patient", "Residual"), 
                     "var1" = c("(Intercept)", "(Intercept)", "(Intercept)", "NA"), 
                     "var2" = c("NA", "NA", "NA", "NA"), 
                     "vcov" = c(0.14345485, 0.16865260, 0.09309796, 0.88206248), 
                     "sdcor" = c(0.3787543, 0.4106733, 0.3051196, 0.9391818))
    test <- ranef(eNRI3.lmm0, effects = "variance", format = "wide", simplify = FALSE)
    expect_equal(as.double(test[-1,"absolute.1"]), as.double(GS[c(3,2,1,4),"vcov"]), tol = 1e-3)
})


##----------------------------------------------------------------------
### test-auto-mixed-model.R ends here
