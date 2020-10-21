### doc-data.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 21 2020 (13:42) 
## Version: 
## Last-Updated: okt 21 2020 (14:44) 
##           By: Brice Ozenne
##     Update #: 10
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * gastricbypassW
#' @title Data From The Gastric Bypass Study (Wide Format)
#'
#' @description  Data from the gastric bypass study
#' where the bodyweight and serum glucagon (a gut hormone) were measured in 20 obese subjects prior and after gastric bypass surgery.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item id. Patient identifier
#' \item weight1. Bodyweight (in kg) 3 months before surgery.
#' \item weight2. Bodyweight (in kg) 1 week before surgery.
#' \item weight3. Bodyweight (in kg) 1 week after surgery.
#' \item weight4. Bodyweight (in kg) 3 months after surgery.
#' \item glucagonAUC1. Glucagon value 3 months before surgery.
#' \item glucagonAUC2. Glucagon value 1 week before surgery.
#' \item glucagonAUC3. Glucagon value 1 week after surgery.
#' \item glucagonAUC4. Glucagon value 3 months after surgery.
#' }
#' 
#' @name gastricbypassW
#' @docType data
#' @usage data(gastricbypassW)
#' @references The effect of Roux-en-Y gastric bypass surgery on the gut mucosal gene expression profile and circulating gut hormones. \url{https://easddistribute.m-anage.com/from.storage?image=4iBH9mRQm1kfeEHULC2CxovdlyCtA1EHeVDdoffnZrAUGG9SHTO-U4ItnLU078eVkF1ZUZgYTy7THlTW3KSgFA2}
#' @keywords data
NULL

#### For reshaping the datasets ####
## library(data.table)
## gastricbypassW$id <- as.factor(gastricbypassW$id)
## save(gastricbypassW, file = "data/gastricbypassW.rda")
##
## gastricbypassL <- melt(as.data.table(gastricbypassW), id.vars = "id", measure.vars = patterns("weight","glucagon"), value.name = c("weight","glucagon"), variable.name = "time")
## gastricbypassL[, visit := time]
## gastricbypassL[, time := factor(time, levels = 1:4, labels = c("3 months before surgery","1 week before surgery","1 week after surgery","3 months after surgery"))]
## gastricbypassL[, id := as.factor(id)]
## setcolorder(gastricbypassL, c("id","visit","time","weight","glucagon"))
## gastricbypassL <- as.data.frame(gastricbypassL)
## save(gastricbypassL, file = "data/gastricbypassL.rda")
#####################################

## * gastricbypassL
#' @title Data From The Gastric Bypass Study (Long Format)
#'
#' @description  Data from the gastric bypass study
#' where the bodyweight and serum glucagon (a gut hormone) were measured in 20 obese subjects prior and after gastric bypass surgery.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item id. Patient identifier
#' \item visit. The visit index.
#' \item time. The time at which the visit took place.
#' \item weight. Bodyweight (in kg) measured during the visit.
#' \item glucagon. Glucagon measured during the visit.
#' }
#' 
#' @name gastricbypassL
#' @docType data
#' @usage data(gastricbypassL)
#' @references The effect of Roux-en-Y gastric bypass surgery on the gut mucosal gene expression profile and circulating gut hormones. \url{https://easddistribute.m-anage.com/from.storage?image=4iBH9mRQm1kfeEHULC2CxovdlyCtA1EHeVDdoffnZrAUGG9SHTO-U4ItnLU078eVkF1ZUZgYTy7THlTW3KSgFA2}
#' @keywords data
NULL

######################################################################
### doc-data.R ends here
