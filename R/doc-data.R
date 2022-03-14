### doc-data.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 21 2020 (13:42) 
## Version: 
## Last-Updated: mar 14 2022 (09:41) 
##           By: Brice Ozenne
##     Update #: 79
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * blandAltman
## ** blandAltmanW
#' @title Data From The Bland Altman Study (Wide Format)
#'
#' @description  Data From The Bland Altman Study where two methods to measure the peak expiratory flow rate (PEFR) where compared.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item id Patient identifier
#' \item wright1 First measurement made with a Wright peak flow meter.
#' \item wright2 Second measurement made with a Wright peak flow meter.
#' \item mini1 First measurement made with a mini Wright peak flow meter.
#' \item mini2 Second measurement made with a mini Wright peak flow meter.
#' }
#' 
#' @name blandAltmanW
#' @docType data
#' @usage data(blandAltmanW)
#' @references Bland & Altman, Statistical methods for assessing agreement between two methods of clinical measurement, Lancet, 1986; i: 307-310.
#' @keywords data
NULL

## ** calciumL
#' @title Data From The Bland Altman Study (Long Format)
#'
#' @description  Data From The Bland Altman Study where two methods to measure the peak expiratory flow rate (PEFR) where compared.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item id Patient identifier.
#' \item replicate Index of the measurement (first or second).
#' \item method Device used to make the measurement (Wright peak flow meter or mini Wright peak flow meter).
#' \item pefr Measurement (peak expiratory flow rate).
#' }
#' 
#' @name blandAltmanL
#' @docType data
#' @usage data(blandAltmanL)
#' @references Bland & Altman, Statistical methods for assessing agreement between two methods of clinical measurement, Lancet, 1986; i: 307-310.
#' @keywords data
NULL

## * bloodpressure
## ** bloodpressureL
#' @title Data From The Blood Pressure Study (Long Format)
#'
#' @description  Data from a cross-over trial comparing the impact of three formulations of a drug on the blood pressure.
#' The study was conducted on 12 male volunteers randomly divided into tree groups
#' and receiving each of the three formulations with a wash-out period of one week.
#'
#' \itemize{
#' \item id Patient identifier
#' \item sequence sequence of treatment 
#' \item treatment Formulation of the treatment:
#' A (50 mg tablet)
#' B (100 mg tablet)
#' C (sustained-release formulation capsule)
#' \item period time period (in weeks)
#' \item duration duration of the drug (in hours)
#' }
#' 
#' @name bloodpressureL
#' @docType data
#' @usage data(bloodpressureL)
#' @references TO ADD
#' @keywords data
NULL

## * calcium
## ** calciumW
#' @title Data From The Calcium Supplements Study (Wide Format)
#'
#' @description  Data from a randomized study including 112 girls at age 11 investigate the effect of a calcium supplement (n=55) vs. placebo (n=57)
#' on bone mineral density over a 2 year follow-up. The clinical question is: does a calcium supplement help to increase bone gain in adolescent women?
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item girl Patient identifier
#' \item grp Treatment group: calcium supplement (coded \code{C}) or placebo (coded \code{P})
#' \item obstime1 Time after the start of the study at which the first visit took place (in years).
#' \item obstime2 Time after the start of the study at which the second visit took place (in years).
#' \item obstime3 Time after the start of the study at which the third visit took place (in years).
#' \item obstime4 Time after the start of the study at which the fourth visit took place (in years).
#' \item obstime5 Time after the start of the study at which the fifth visit took place (in years).
#' \item bmd1 Bone mineral density measured at the first visit (in mg/cm3).
#' \item bmd2 Bone mineral density measured at the second visit (in mg/cm3).
#' \item bmd3 Bone mineral density measured at the third visit (in mg/cm3).
#' \item bmd4 Bone mineral density measured at the fourth visit (in mg/cm3).
#' \item bmd5 Bone mineral density measured at the fifth visit (in mg/cm3).
#' }
#' 
#' @name calciumW
#' @docType data
#' @usage data(calciumW)
#' @references Vonesh and Chinchilli 1997. Linear and Nonlinear models for the analysis of repeated measurement (Table 5.4.1 on page 228). New York: Marcel Dekker.
#' @keywords data
NULL
## calciumW <- read.table("inst/dataTXT/calcium1.txt", header = TRUE, na.string = ".")
## calciumW$obstime1 <- as.numeric(calciumW$obstime1)
## calciumW$dropout <- NULL
## calciumW$girl <- as.factor(calciumW$girl)
## calciumW$grp <- as.factor(calciumW$grp)
## save(calciumW, file = "data/calciumW.rda")

## ** calciumL
#' @title Data From The Calcium Supplements Study (Long Format)
#'
#' @description  Data from a randomized study including 112 girls at age 11 investigate the effect of a calcium supplement (n=55) vs. placebo (n=57)
#' on bone mineral density over a 2 year follow-up. The clinical question is: does a calcium supplement help to increase bone gain in adolescent women?
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item girl Patient identifier
#' \item grp Treatment group: calcium supplement (coded \code{C}) or placebo (coded \code{P})
#' \item visit Visit index
#' \item bmd Bone mineral density (mg/cm3)
#' \item time.obs Visit time (in years)
#' \item time.num Scheduled visit time (numeric variable, in years)
#' \item time.fac Scheduled visit time (factor variable)
#' }
#' 
#' @name calciumL
#' @docType data
#' @usage data(calciumL)
#' @references TO ADD
#' @keywords data
NULL
## data("calciumW")
## dtW <- data.table::as.data.table(calciumW)
## dtL <- data.table::melt(dtW, id.vars = c("girl","grp"),
##                         measure.vars = patterns("obstime","bmd1"),
##                         value.name = c("time.obs","bmd"), variable.name = "visit")
## calciumL <- as.data.frame(dtL)
## calciumL$time.num <- sapply(as.character(calciumL$visit), switch,
##                             "1" = 0.0,
##                             "2" = 0.5,
##                             "3" = 1.0,
##                             "4" = 1.5,
##                             "5" = 2.0)
## calciumL$time.fac <- factor(calciumL$visit, levels = 1:5,
##                             labels = c("0 years","0.5 years","1 years","1.5 years","2 years")) 
## save(calciumL, file = "data/gastricbypassL.rda")
##
## str(calciumL)

## * ckd
## ** ckdW
#' @title CKD wide
#'
#' @description TODO
#'
#' \itemize{
#' \item id Patient identifier
#' \item allocation
#' \item sex
#' \item age
#' \item pwv0
#' \item pwv12
#' \item pwv24
#' \item aix0
#' \item aix12
#' \item aix24
#' \item dropout
#' }
#' 
#' @name ckdW
#' @docType data
#' @usage data(ckdW)
#' @references TO ADD
#' @keywords data
NULL

## ** ckdL
#' @title CKD long
#'
#' @description TODO
#' 
#' \itemize{
#' \item id Patient identifier
#' \item allocation
#' \item sex
#' \item age
#' \item visit
#' \item time
#' \item pwv
#' \item aix
#' \item dropout
#' }
#' @name ckdL
#' @docType data
#' @usage data(ckdL)
#' @references TO ADD
#' @keywords data
NULL

## * gastricbypass

## ** gastricbypassW
#' @title Data From The Gastric Bypass Study (Wide Format)
#'
#' @description  Data from the gastric bypass study
#' where the bodyweight and serum glucagon (a gut hormone) were measured in 20 obese subjects prior and after gastric bypass surgery.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item id Patient identifier
#' \item weight1 Bodyweight (in kg) 3 months before surgery.
#' \item weight2 Bodyweight (in kg) 1 week before surgery.
#' \item weight3 Bodyweight (in kg) 1 week after surgery.
#' \item weight4 Bodyweight (in kg) 3 months after surgery.
#' \item glucagonAUC1 Glucagon value 3 months before surgery.
#' \item glucagonAUC2 Glucagon value 1 week before surgery.
#' \item glucagonAUC3 Glucagon value 1 week after surgery.
#' \item glucagonAUC4 Glucagon value 3 months after surgery.
#' }
#' 
#' @name gastricbypassW
#' @docType data
#' @usage data(gastricbypassW)
#' @references The effect of Roux-en-Y gastric bypass surgery on the gut mucosal gene expression profile and circulating gut hormones. \url{https://easddistribute.m-anage.com/from.storage?image=4iBH9mRQm1kfeEHULC2CxovdlyCtA1EHeVDdoffnZrAUGG9SHTO-U4ItnLU078eVkF1ZUZgYTy7THlTW3KSgFA2}
#' @keywords data
NULL

## ** gastricbypassL
#' @title Data From The Gastric Bypass Study (Long Format)
#'
#' @description  Data from the gastric bypass study
#' where the bodyweight and serum glucagon (a gut hormone) were measured in 20 obese subjects prior and after gastric bypass surgery.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item id Patient identifier
#' \item visit The visit index.
#' \item time The time at which the visit took place.
#' \item weight Bodyweight (in kg) measured during the visit.
#' \item glucagonAUC Glucagon measured during the visit.
#' }
#' 
#' @name gastricbypassL
#' @docType data
#' @usage data(gastricbypassL)
#' @references The effect of Roux-en-Y gastric bypass surgery on the gut mucosal gene expression profile and circulating gut hormones. \url{https://easddistribute.m-anage.com/from.storage?image=4iBH9mRQm1kfeEHULC2CxovdlyCtA1EHeVDdoffnZrAUGG9SHTO-U4ItnLU078eVkF1ZUZgYTy7THlTW3KSgFA2}
#' @keywords data
NULL
## data("gastricbypassW")
## dtW <- data.table::as.data.table(gastricbypassW)
## dtL <- data.table::melt(dtW, id.vars = "id",
##                         measure.vars = patterns("weight","glucagonAUC"),
##                         value.name = c("weight","glucagonAUC"), variable.name = "time")
## gastricbypassL <- as.data.frame(dtL)
## gastricbypassL$visit <- gastricbypassL$time
## gastricbypassL$time <- factor(gastricbypassL$visit, levels = 1:4,
##                               labels = c("3 months before","1 week before",
##                                          "1 week after","3 months after"))
## gastricbypassL <- gastricbypassL[,c("id","visit","time","weight","glucagonAUC")]
## save(gastricbypassL, file = "data/gastricbypassL.rda")
##
## str(gastricbypassL)

## * ncgs
## ** ncgsW
#' @title Data From National Cooperative Gallstone Study (Wide Format)
#'
#' @description  Data from the National Cooperative Gallstone Study (NCGS),
#' a randomized study where the level of serum cholesterol was measured at baseline and after intake of high-dose chenondiol (750mg/day) or placebo.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item group Treatment group: highdose or placebo.
#' \item id Patient identifier
#' \item cholest1 cholesterol measurement at baseline (before treatment).
#' \item cholest2 cholesterol measurement at 6 months (after treatment).
#' \item cholest3 cholesterol measurement at 12 months (after treatment).
#' \item cholest4 cholesterol measurement at 20 months (after treatment).
#' \item cholest5 cholesterol measurement at 24 months (after treatment).
#' }
#' 
#' @name ncgsW
#' @docType data
#' @usage data(ncgsW)
#' @references Grundy SM, Lan SP, Lachin J. The effects of chenodiol on biliary lipids and their association with gallstone dissolution in the National Cooperative Gallstone Study (NCGS). J Clin Invest. 1984 Apr;73(4):1156-66. doi: 10.1172/JCI111301.  
#' @keywords data
NULL
## ncgsW <- read.table("inst/dataTXT/ncgs.txt", header = TRUE, na.string = ".")
## ncgsW$group <- as.factor(ncgsW$group)
## ncgsW$id <- as.factor(ncgsW$id)
## save(ncgsW, file = "data/ncgsW.rda")
##
## str(ncgsW)

## ** ncgsL
#' @title Data From National Cooperative Gallstone Study (Long Format)
#'
#' @description  Data from the National Cooperative Gallstone Study (NCGS),
#' a randomized study where the level of serum cholesterol was measured at baseline and after intake of high-dose chenondiol (750mg/day) or placebo.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item group Treatment group: highdose or placebo.
#' \item id Patient identifier
#' \item visit visit index.
#' \item cholest cholesterol measurement.
#' \item time time after the start of the study at which the measurement has been done (in month). Treatment is given at 0+.
#' }
#' 
#' @name ncgsL
#' @docType data
#' @usage data(ncgsL)
#' @references Grundy SM, Lan SP, Lachin J. The effects of chenodiol on biliary lipids and their association with gallstone dissolution in the National Cooperative Gallstone Study (NCGS). J Clin Invest. 1984 Apr;73(4):1156-66. doi: 10.1172/JCI111301.  
#' @keywords data
NULL
## data("ncgsW")
## ncgsL <- reshape2::melt(ncgsW, id.vars = c("group","id"),
##                         measure.vars = paste0("cholest",1:5),
##                         value.name = c("cholest"), variable.name = "visit")
## ncgsL$visit <- as.factor(sapply(ncgsL$visit, gsub,
##                       pattern = "cholest", replacement = ""))
## ncgsL$time <- sapply(as.character(ncgsL$visit), switch,
##                             "1" = 0,
##                             "2" = 6,
##                             "3" = 12,
##                             "4" = 18,
##                             "5" = 24)
## save(ncgsL, file = "data/ncgsL.rda")
##
## str(ncgsL)

## * potassium
## ** potassiumSingleW
#' @title Data From The Potassium Intake Study (Wide Format)
#'
#' @description  Data from the potassium intake study,
#' a randomized placebo-controlled crossover study where the effect of potassium supplement (90 mmol/day) on the renin-angiostensin-aldosteron system (RAAS) was assessed.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item id Patient identifier
#' \item sequence Treatment group to which the patient has been randomized.
#' \item treatment1 Treatment during the first time period.
#' \item treatment2 Treatment during the second time period
#' \item auc1 Area under the curve of ?? during the first time period
#' \item auc2 Area under the curve of ?? during the second time period 
#' \item bsauc1 ??
#' \item aldo1 ??
#' \item aldo2 ??
#' }
#' 
#' @name potassiumSingleW
#' @docType data
#' @usage data(potassiumSingleW)
#' @references Dreier et al. Effect of increased potassium intake on the reninangiotensinaldosterone system and subcutaneous resistance arteries: a randomized crossover study,
#' Nephrol Dial Transplant (2020) 110. doi: 10.1093/ndt/gfaa114
#' @keywords data
NULL

## ** potassiumSingleL
#' @title Data From The Potassium Intake Study (Long Format)
#'
#' @description  Data from the potassium intake study,
#' a randomized placebo-controlled crossover study where the effect of potassium supplement (90 mmol/day) on the renin-angiostensin-aldosteron system (RAAS) was assessed.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item id Patient identifier
#' \item sequence Treatment group to which the patient has been randomized.
#' \item period Time period.
#' \item treatment Treatment during the time period
#' \item auc Area under the curve of ?? during the time period
#' \item bsauc ??
#' \item aldo ??
#' }
#' 
#' @name potassiumSingleL
#' @docType data
#' @usage data(potassiumSingleL)
#' @references Dreier et al. Effect of increased potassium intake on the reninangiotensinaldosterone system and subcutaneous resistance arteries: a randomized crossover study,
#' Nephrol Dial Transplant (2020) 110. doi: 10.1093/ndt/gfaa114
#' @keywords data
NULL

## ** potassiumRepeatedL
#' @title Data From The Potassium Intake Study (Long Format with intermediate measurements)
#'
#' @description  Data from the potassium intake study,
#' a randomized placebo-controlled crossover study where the effect of potassium supplement (90 mmol/day) on the renin-angiostensin-aldosteron system (RAAS) was assessed.
#' This dataset is in the long format (i.e. one line per measurement) and contains measurement over 6 timepoints for each time period.
#'
#' \itemize{
#' \item id Patient identifier
#' \item sequence Treatment group to which the patient has been randomized.
#' \item period Time period.
#' \item treatment Treatment during the time period
#' \item time Time within each period
#' \item aldo ??
#' }
#' 
#' @name potassiumRepeatedL
#' @docType data
#' @usage data(potassiumRepeatedL)
#' @references Dreier et al. Effect of increased potassium intake on the reninangiotensinaldosterone system and subcutaneous resistance arteries: a randomized crossover study,
#' Nephrol Dial Transplant (2020) 110. doi: 10.1093/ndt/gfaa114
#' @keywords data
NULL

## * swabs
## ** swabsW
#' @title Data From The SWABS Study (Wide Format)
#'
#' @description Data from the swabs study,
#' where the pneumococcus was studied in 18 families with different space available for the household.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item crowding Space available in the household.
#' \item family Family serial number
#' \item mother number of times the swab measurement was positive for the mother.
#' \item father number of times the swab measurement was positive for the father.
#' \item child1 number of times the swab measurement was positive for the first child.
#' \item child2 number of times the swab measurement was positive for the second child.
#' \item child3 number of times the swab measurement was positive for the third child.
#' }
#' 
#' @name swabsW
#' @docType data
#' @usage data(swabsW)
#' @references Grundy SM, Lan SP, Lachin J. The effects of chenodiol on biliary lipids and their association with gallstone dissolution in the National Cooperative Gallstone Study (SWABS). J Clin Invest. 1984 Apr;73(4):1156-66. doi: 10.1172/JCI111301.  
#' @keywords data
NULL
## library(reshape2)
## data(swabsL)
## swabsW <- dcast(swabsL, formula = crowding+family~name, value.var = "swabs")
## save(swabsW, file = "data/swabsW.rda")
## str(swabsW)

## ** swabsL
#' @title Data From The SWABS Study (Long Format)
#'
#' @description  Data from the swabs study,
#' where the pneumococcus was studied in 18 families with different space available for the household.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item crowding Space available in the household.
#' \item family Family serial number
#' \item name Type of family member.
#' \item swabs number of times the swab measurement was positive.
#' }
#' 
#' @name swabsL
#' @docType data
#' @usage data(swabsL)
#' @references TODO
#' @keywords data
NULL
## swabsL <- read.table("inst/dataTXT/swabs.txt", header = TRUE, na.string = ".")
## swabsL$family <- as.factor(swabsL$family)
## swabsL$crowding <- factor(swabsL$crowding, 
##             levels = c("uncrow","crow","overcrow"))
## swabsL$name <- factor(swabsL$name,
##    levels = c("mother","father", "child1","child2","child3"))
## swabsL <- swabsL[order(swabsL$crowding,swabsL$family,swabsL$name),]
## save(swabsL, file = "data/swabsL.rda")
## str(swabsL)


## * vasscores
## ** vasscoresW
#' @title Data From The VAS Study (Wide Format)
#'
#' @description  Data from the VAS Study,
#' a randomized controlled clinial trial assessing the healing effect of topical zink sulfate on epidermal wound.
#' The study includes 30 heatlhy volunteers with induced wounds on each buttock which where subsequently treated with a different treatment for each wound.
#' Then the VAS-score (pain sensation on a 0-100mm visual analogue scale) was assessed after each treatment application and summarized by area under the curve.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item id Patient identifier.
#' \item group Treatment group to which the patient has been randomized.
#' \item vasA VAS-score when using a zink shower gel.
#' \item vasB VAS-score when using a placebo treatment (shower gel without zink).
#' \item vasC VAS-score when using a control treatment with demineralized water.
#' }
#' 
#' @name vasscoresW
#' @docType data
#' @usage data(vasscoresW)
#' @references TODO
#' @keywords data
NULL
## vasscoresW <- read.table("inst/dataTXT/vasscores.txt", header = TRUE, na.string = ".")
## vasscoresW$id <- factor(vasscoresW$id)
## vasscoresW$group <- factor(vasscoresW$group)
## save(vasscoresW, file = "data/vasscoresW.rda")

## ** vasscoresL
#' @title Data From The VAS Study (Long Format)
#'
#' @description  Data from the VAS Study,
#' a randomized controlled clinial trial assessing the healing effect of topical zink sulfate on epidermal wound.
#' The study includes 30 heatlhy volunteers with induced wounds on each buttock which where subsequently treated with a different treatment for each wound.
#' Then the VAS-score (pain sensation on a 0-100mm visual analogue scale) was assessed after each treatment application and summarized by area under the curve.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item id Patient identifier.
#' \item group Treatment group to which the patient has been randomized.
#' \item treat.num 
#' \item vas VAS-score relative to the wound.
#' \item treatment Treatment used on the wound.
#' A: active treatment (zink shower gel),
#' B: placebo treatment (shower gel without zink),
#' C: control treatment (demineralized water).
#' }
#' 
#' @name vasscoresL
#' @docType data
#' @usage data(vasscoresL)
#' @references TODO
#' @keywords data
NULL
## data("vasscoresW")
## ## transform to long format:
## vasscoresL <- reshape(vasscoresW, 
##                direction="long", 
##                idvar=c("id","group"), 
##                varying=c("vasA","vasB","vasC"),
##                v.names=c("vas"),
##                timevar="treat.num")
##
## ## Make a categorical version of the treatment variable:
## vasscoresL$treatment <- factor(vasscoresL$treat.num, labels=c('A','B','C'))
## ## Fix attributes
## rownames(vasscoresL) <- NULL
## attr(vasscoresL, "reshapeLong") <- NULL
## ## Export
## save(vasscoresL, file="data/vasscoresL.rda")
## str(vasscoresL)

## * vitamin
## ** vitaminW
#' @title Data From The Vitamin Study (Wide Format)
#'
#' @description  Data from the vitamin Study,
#' a randomized study where the growth of guinea pigs was monitored before and after intake of vitamin E/placebo.
#' The weight of each guinea pig was recorded at the end of week 1, 3, 4, 5, 6, and 7. Vitamin E/placebo is given at the beginning of week 5.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item group Treatment group: vitamin or placebo.
#' \item animal Identifier
#' \item weight1 weight (in g) of the pig at the end of week 1 (before treatment).
#' \item weight3 weight (in g) of the pig at the end of week 3 (before treatment).
#' \item weight4 weight (in g) of the pig at the end of week 4 (before treatment).
#' \item weight5 weight (in g) of the pig at the end of week 5 (after treatment).
#' \item weight6 weight (in g) of the pig at the end of week 6 (after treatment).
#' \item weight7 weight (in g) of the pig at the end of week 7 (after treatment).
#' }
#' 
#' @name vitaminW
#' @docType data
#' @usage data(vitaminW)
#' @references TODO
#' @keywords data
NULL
## vitaminW <- read.table("inst/dataTXT/vitamin.txt", header = TRUE, na.string = ".")
## vitaminW$group <- as.factor(vitaminW$group)
## vitaminW$animal <- as.factor(vitaminW$animal)
## save(vitaminW, file = "data/vitaminW.rda")

## ** vitaminL
#' @title Data From The Vitamin Study (Long Format)
#'
#' @description  Data from the vitamin Study,
#' a randomized study where the growth of guinea pigs was monitored before and after intake of vitamin E/placebo.
#' The weight of each guinea pig was recorded at the end of week 1, 3, 4, 5, 6, and 7. Vitamin E/placebo is given at the beginning of week 5.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item group Treatment group: vitamin or placebo.
#' \item animal Identifier
#' \item weight1 weight (in g) of the pig at the end of week 1 (before treatment).
#' \item weight3 weight (in g) of the pig at the end of week 3 (before treatment).
#' \item weight4 weight (in g) of the pig at the end of week 4 (before treatment).
#' \item weight5 weight (in g) of the pig at the end of week 5 (after treatment).
#' \item weight6 weight (in g) of the pig at the end of week 6 (after treatment).
#' \item weight7 weight (in g) of the pig at the end of week 7 (after treatment).
#' }
#' 
#' @name vitaminL
#' @docType data
#' @usage data(vitaminL)
#' @references Crowder and Hand (1990, p. 27) Analysis of Repeated Measures.
#' @keywords data
NULL
## data("vitaminW")
## vitaminL <- reshape2::melt(vitaminW, id.vars = c("group","animal"),
##                         measure.vars = paste0("weight",c(1,3:7)),
##                         value.name = c("weight"), variable.name = "visit")
## vitaminL$visit <- as.factor(as.numeric(as.factor(sapply(vitaminL$visit, gsub,
##                              pattern = "weight", replacement = ""))))
## vitaminL$time <- sapply(as.character(vitaminL$visit), switch,
##                             "1" = 1,
##                             "2" = 3,
##                             "3" = 4,
##                             "4" = 5,
##                             "5" = 6,
##                             "6" = 7)
## save(vitaminL, file = "data/vitaminL.rda")


######################################################################
### doc-data.R ends here
