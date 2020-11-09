### doc-data.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 21 2020 (13:42) 
## Version: 
## Last-Updated: nov  9 2020 (10:49) 
##           By: Brice Ozenne
##     Update #: 39
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

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
#' @examples
#' \dontrun{
#' calciumW <- read.table("inst/dataTXT/calcium1.txt", header = TRUE, na.string = ".")
#' calciumW$obstime1 <- as.numeric(calciumW$obstime1)
#' calciumW$dropout <- NULL
#' calciumW$girl <- as.factor(calciumW$girl)
#' calciumW$grp <- as.factor(calciumW$grp)
#' save(calciumW, file = "data/calciumW.rda")
#' }
NULL

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
#' @examples
#' \dontrun{
#' data("calciumW")
#' dtW <- data.table::as.data.table(calciumW)
#' dtL <- data.table::melt(dtW, id.vars = c("girl","grp"),
#'                         measure.vars = patterns("obstime","bmd1"),
#'                         value.name = c("time.obs","bmd"), variable.name = "visit")
#' calciumL <- as.data.frame(dtL)
#' calciumL$time.num <- sapply(as.character(calciumL$visit), switch,
#'                             "1" = 0.0,
#'                             "2" = 0.5,
#'                             "3" = 1.0,
#'                             "4" = 1.5,
#'                             "5" = 2.0)
#' calciumL$time.fac <- factor(calciumL$visit, levels = 1:5, labels = c("0 years","0.5 years","1 years","1.5 years","2 years")) 
#' save(calciumL, file = "data/gastricbypassL.rda")
#' }
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
#' @examples
#' \dontrun{
#' gastricbypassW <- read.table("inst/dataTXT/gastricbypass.txt", header = TRUE, na.string = ".")
#' gastricbypassW$id <- as.factor(gastricbypassW$id)
#' save(gastricbypassW, file = "data/gastricbypassW.rda")
#' }
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
#' \item glucagon Glucagon measured during the visit.
#' }
#' 
#' @name gastricbypassL
#' @docType data
#' @usage data(gastricbypassL)
#' @references The effect of Roux-en-Y gastric bypass surgery on the gut mucosal gene expression profile and circulating gut hormones. \url{https://easddistribute.m-anage.com/from.storage?image=4iBH9mRQm1kfeEHULC2CxovdlyCtA1EHeVDdoffnZrAUGG9SHTO-U4ItnLU078eVkF1ZUZgYTy7THlTW3KSgFA2}
#' @keywords data
#' @examples
#' \dontrun{
#' data("gastricbypassW")
#' dtW <- data.table::as.data.table(gastricbypassW)
#' dtL <- data.table::melt(dtW, id.vars = "id",
#'                         measure.vars = patterns("weight","glucagonAUC"),
#'                         value.name = c("weight","glucagon"), variable.name = "time")
#' gastricbypassL <- as.data.frame(dtL)
#' gastricbypassL$visit <- gastricbypassL$time
#' gastricbypassL$time <- factor(gastricbypassL$visit, levels = 1:4, labels = c("3 months before surgery","1 week before surgery","1 week after surgery","3 months after surgery"))
#' gastricbypassL <- gastricbypassL[,c("id","visit","time","weight","glucagon")]
#' save(gastricbypassL, file = "data/gastricbypassL.rda")
#' }
NULL

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
#' @examples
#' \dontrun{
#' ncgsW <- read.table("inst/dataTXT/ncgs.txt", header = TRUE, na.string = ".")
#' ncgsW$group <- as.factor(ncgsW$group)
#' ncgsW$id <- as.factor(ncgsW$id)
#' save(ncgsW, file = "data/ncgsW.rda")
#' }
NULL

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
#' @examples
#' \dontrun{
#' data("ncgsW")
#' ncgsL <- reshape2::melt(ncgsW, id.vars = c("group","id"),
#'                         measure.vars = paste0("cholest",1:5),
#'                         value.name = c("cholest"), variable.name = "visit")
#' ncgsL$visit <- as.factor(sapply(ncgsL$visit, gsub,
#'                       pattern = "cholest", replacement = ""))
#' ncgsL$time <- sapply(as.character(ncgsL$visit), switch,
#'                             "1" = 0,
#'                             "2" = 6,
#'                             "3" = 12,
#'                             "4" = 18,
#'                             "5" = 24)
#' save(ncgsL, file = "data/ncgsL.rda")
#' }
NULL

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
#' \item weigth1 weight (in g) of the pig at the end of week 1 (before treatment).
#' \item weigth3 weight (in g) of the pig at the end of week 3 (before treatment).
#' \item weigth4 weight (in g) of the pig at the end of week 4 (before treatment).
#' \item weigth5 weight (in g) of the pig at the end of week 5 (after treatment).
#' \item weigth6 weight (in g) of the pig at the end of week 6 (after treatment).
#' \item weigth7 weight (in g) of the pig at the end of week 7 (after treatment).
#' }
#' 
#' @name vitaminW
#' @docType data
#' @usage data(vitaminW)
#' @references TODO
#' @keywords data
#' @examples
#' \dontrun{
#' vitaminW <- read.table("inst/dataTXT/vitamin.txt", header = TRUE, na.string = ".")
#' vitaminW$group <- as.factor(vitaminW$group)
#' vitaminW$animal <- as.factor(vitaminW$animal)
#' save(vitaminW, file = "data/vitaminW.rda")
#' }
NULL

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
#' \item weigth1 weight (in g) of the pig at the end of week 1 (before treatment).
#' \item weigth3 weight (in g) of the pig at the end of week 3 (before treatment).
#' \item weigth4 weight (in g) of the pig at the end of week 4 (before treatment).
#' \item weigth5 weight (in g) of the pig at the end of week 5 (after treatment).
#' \item weigth6 weight (in g) of the pig at the end of week 6 (after treatment).
#' \item weigth7 weight (in g) of the pig at the end of week 7 (after treatment).
#' }
#' 
#' @name vitaminL
#' @docType data
#' @usage data(vitaminL)
#' @references Crowder and Hand (1990, p. 27) Analysis of Repeated Measures.
#' @keywords data
#' @examples
#' \dontrun{
#' data("vitaminW")
#' vitaminL <- reshape2::melt(vitaminW, id.vars = c("group","animal"),
#'                         measure.vars = paste0("weight",c(1,3:7)),
#'                         value.name = c("weight"), variable.name = "visit")
#' vitaminL$visit <- as.factor(as.numeric(as.factor(sapply(vitaminL$visit, gsub,
#'                              pattern = "weight", replacement = ""))))
#' vitaminL$time <- sapply(as.character(vitaminL$visit), switch,
#'                             "1" = 1,
#'                             "2" = 3,
#'                             "3" = 4,
#'                             "4" = 5,
#'                             "5" = 6,
#'                             "6" = 7)
#' save(vitaminL, file = "data/vitaminL.rda")
#' }
NULL


######################################################################
### doc-data.R ends here
