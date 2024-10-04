### doc-data.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 21 2020 (13:42) 
## Version: 
## Last-Updated: sep 30 2024 (13:59) 
##           By: Brice Ozenne
##     Update #: 160
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * abeta
## ** abetaW
#' @title Data From The abeta Study (Wide Format)
#' @name abetaW
#' @rdname data-abetaW
#'
#' @description Extract data from a longitudinal case control study including 87 patients newly diagnosed with bipolar disorder and 44 age and sex matched healthy controls.
#' Contains demographic data and lifestyle factors at baseline, as well as measures of psychosocial functioning at baseline and 1 year follow-up.
#' This dataset is in the wide format (i.e. one line per participant).
#'
#' \itemize{
#' \item \code{id}: study participant.
#' \item \code{sex}: male (M) or female (F).
#' \item \code{age}: age in years.
#' \item \code{group}: bipolar disorder (BD) or healthy control (HC).
#' \item \code{episode}: whether the patient experience an affective episode during follow-up.
#' \item \code{fast0},\code{fast1}: functioning assessment short test at baseline and follow-up.
#' \item \code{qol0},\code{qol1}: WHO quality of life score at baseline and follow-up.
#' \item \code{pss0},\code{pss1}: perceived stress score at baseline and follow-up.
#' \item \code{educationyears}: years of education including basic school.
#' \item \code{alcohol}: daily alcohol consumption.
#' \item \code{missingreason}: reason of drop out or missed visit.
#' }
#' 
#' @docType data
#' @usage data(abetaW)
#' @references Pech, Josefine, et al. "The impact of a new affective episode on psychosocial functioning, quality of life and perceived stress in newly diagnosed patients with bipolar disorder: A prospective one-year case-control study."Journal of Affective Disorders 277 (2020): 486-494.
#' @keywords datasets
NULL

## ** abetaL
#' @title Data From The Bland Altman Study (Long Format)
#' @name abetaL
#' @rdname data-abetaL
#'
#' @description Extract data from a longitudinal case control study including 87 patients newly diagnosed with bipolar disorder and 44 age and sex matched healthy controls.
#' Contains demographic data and lifestyle factors at baseline, as well as measures of psychosocial functioning at baseline and 1 year follow-up.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item \code{id}: study participant.
#' \item \code{sex}: male (M) or female (F).
#' \item \code{age}: age in years.
#' \item \code{group}: bipolar disorder (BD) or healthy control (HC).
#' \item \code{episode}: whether the patient experience an affective episode during follow-up.
#' \item \code{visit}: index of time at which pss, fast, and qol measurements where performed.
#' \item \code{year}: time at which pss, fast, and qol measurements where performed.
#' \item \code{pss}: perceived stress score.
#' \item \code{fast}: functioning assessment short test.
#' \item \code{qol}: WHO quality of life score.
#' \item \code{educationyears}: years of education including basic school.
#' \item \code{alcohol}: daily alcohol consumption.
#' \item \code{missingreason}: reason of drop out or missed visit.
#' }
#' 
#' @docType data
#' @usage data(abetaL)
#' @references Pech, Josefine, et al. The impact of a new affective episode on psychosocial functioning, quality of life and perceived stress in newly diagnosed patients with bipolar disorder: A prospective one-year case-control study.Journal of Affective Disorders 277 (2020): 486-494.
#' @keywords datasets
NULL

## * blandAltman
## ** blandAltmanW
#' @title Data From The Bland Altman Study (Wide Format)
#' @name blandAltmanW
#' @rdname data-blandAltmanW
#'
#' @description  Data From The Bland Altman Study where two methods to measure the peak expiratory flow rate (PEFR) where compared.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{wright1}: first measurement made with a Wright peak flow meter.
#' \item \code{wright2}: second measurement made with a Wright peak flow meter.
#' \item \code{mini1}: first measurement made with a mini Wright peak flow meter.
#' \item \code{mini2}: second measurement made with a mini Wright peak flow meter.
#' }
#' 
#' @docType data
#' @usage data(blandAltmanW)
#' @references Bland & Altman, Statistical methods for assessing agreement between two methods of clinical measurement, Lancet, 1986; i: 307-310.
#' @keywords datasets
NULL

## ** blandAltmanL
#' @title Data From The Bland Altman Study (Long Format)
#' @name blandAltmanL
#' @rdname data-blandAltmanL
#'
#' @description  Data From The Bland Altman Study where two methods to measure the peak expiratory flow rate (PEFR) where compared.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{replicate}: index of the measurement (first or second).
#' \item \code{method}: device used to make the measurement (Wright peak flow meter or mini Wright peak flow meter).
#' \item \code{pefr}: measurement (peak expiratory flow rate).
#' }
#' 
#' @docType data
#' @usage data(blandAltmanL)
#' @references Bland & Altman, Statistical methods for assessing agreement between two methods of clinical measurement, Lancet, 1986; i: 307-310.
#' @keywords datasets
NULL

## * bloodpressure
## ** bloodpressureL
#' @title Data From The Blood Pressure Study (Long Format)
#' @name bloodpressureL
#' @rdname data-bloodpressureL
#'
#' @description  Data from a cross-over trial comparing the impact of three formulations of a drug on the blood pressure.
#' The study was conducted on 12 male volunteers randomly divided into tree groups
#' and receiving each of the three formulations with a wash-out period of one week.
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{sequence}: sequence of treatment .
#' \item \code{treatment}: formulation of the treatment
#' A (50 mg tablet)
#' B (100 mg tablet)
#' C (sustained-release formulation capsule)
#' \item \code{period}: time period (in weeks).
#' \item \code{duration}: duration of the drug (in hours).
#' }
#' 
#' @docType data
#' @usage data(bloodpressureL)
#' @references TO ADD
#' @keywords datasets
NULL

## * calcium
## ** calciumW
#' @title Data From The Calcium Supplements Study (Wide Format)
#' @name calciumW
#' @rdname data-calciumW
#'
#' @description  Data from a randomized study including 112 girls at age 11 investigate the effect of a calcium supplement (n=55) vs. placebo (n=57)
#' on bone mineral density over a 2 year follow-up. The clinical question is: does a calcium supplement help to increase bone gain in adolescent women?
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item \code{girl}: patient identifier
#' \item \code{grp}: treatment group: calcium supplement (coded \code{C}) or placebo (coded \code{P}).
#' \item \code{obstime1}: time after the start of the study at which the first visit took place (in years).
#' \item \code{obstime2}: time after the start of the study at which the second visit took place (in years).
#' \item \code{obstime3}: time after the start of the study at which the third visit took place (in years).
#' \item \code{obstime4}: time after the start of the study at which the fourth visit took place (in years).
#' \item \code{obstime5}: time after the start of the study at which the fifth visit took place (in years).
#' \item \code{bmd1}: bone mineral density measured at the first visit (in mg/cm3).
#' \item \code{bmd2}: bone mineral density measured at the second visit (in mg/cm3).
#' \item \code{bmd3}: bone mineral density measured at the third visit (in mg/cm3).
#' \item \code{bmd4}: bone mineral density measured at the fourth visit (in mg/cm3).
#' \item \code{bmd5}: bone mineral density measured at the fifth visit (in mg/cm3).
#' }
#' 
#' @docType data
#' @usage data(calciumW)
#' @references Vonesh and Chinchilli 1997. Linear and Nonlinear models for the analysis of repeated measurement (Table 5.4.1 on page 228). New York: Marcel Dekker.
#' @keywords datasets
NULL

## ** calciumL
#' @title Data From The Calcium Supplements Study (Long Format)
#' @name calciumL
#' @rdname data-calciumL
#'
#' @description  Data from a randomized study including 112 girls at age 11 investigate the effect of a calcium supplement (n=55) vs. placebo (n=57)
#' on bone mineral density over a 2 year follow-up. The clinical question is: does a calcium supplement help to increase bone gain in adolescent women?
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item \code{girl}: patient identifier.
#' \item \code{grp}: treatment group: calcium supplement (coded \code{C}) or placebo (coded \code{P}).
#' \item \code{visit}: visit index.
#' \item \code{bmd}: bone mineral density (mg/cm3).
#' \item \code{time.obs}: visit time (in years).
#' \item \code{time.num}: scheduled visit time (numeric variable, in years).
#' \item \code{time.fac}: scheduled visit time (factor variable).
#' }
#' 
#' @docType data
#' @usage data(calciumL)
#' @references TO ADD
#' @keywords datasets
NULL

## * ckd
## ** ckdW
#' @title CKD wide
#' @name ckdW
#' @rdname data-ckdW
#'
#' @description TODO
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{allocation}:
#' \item \code{sex}:
#' \item \code{age}:
#' \item \code{pwv0}:
#' \item \code{pwv12}:
#' \item \code{pwv24}:
#' \item \code{aix0}:
#' \item \code{aix12}:
#' \item \code{aix24}:
#' \item \code{dropout}:
#' }
#' 
#' @docType data
#' @usage data(ckdW)
#' @references TO ADD
#' @keywords datasets
NULL

## ** ckdL
#' @title CKD long
#' @name ckdL
#' @rdname data-ckdL
#'
#' @description TODO
#' 
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{allocation}:
#' \item \code{sex}:
#' \item \code{age}:
#' \item \code{visit}:
#' \item \code{time}:
#' \item \code{pwv}:
#' \item \code{aix}:
#' \item \code{dropout}:
#' }
#' @docType data
#' @usage data(ckdL)
#' @references TO ADD
#' @keywords datasets
NULL

## * gastricbypass
## ** gastricbypassW
#' @title Data From The Gastric Bypass Study (Wide Format)
#' @name gastricbypassW
#' @rdname data-gastricbypassW
#'
#' @description  Data from the gastric bypass study
#' where the bodyweight and serum glucagon (a gut hormone) were measured in 20 obese subjects prior and after gastric bypass surgery.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{weight1}: bodyweight (in kg) 3 months before surgery.
#' \item \code{weight2}: bodyweight (in kg) 1 week before surgery.
#' \item \code{weight3}: bodyweight (in kg) 1 week after surgery.
#' \item \code{weight4}: bodyweight (in kg) 3 months after surgery.
#' \item \code{glucagonAUC1}: glucagon  (in pmol/l x hours) 3 months before surgery.
#' \item \code{glucagonAUC2}: glucagon  (in pmol/l x hours) 1 week before surgery.
#' \item \code{glucagonAUC3}: glucagon  (in pmol/l x hours) 1 week after surgery.
#' \item \code{glucagonAUC4}: glucagon  (in pmol/l x hours) 3 months after surgery.
#' }
#' 
#' @docType data
#' @usage data(gastricbypassW)
#' @references The effect of Roux-en-Y gastric bypass surgery on the gut mucosal gene expression profile and circulating gut hormones. \url{https://easddistribute.m-anage.com/from.storage?image=4iBH9mRQm1kfeEHULC2CxovdlyCtA1EHeVDdoffnZrAUGG9SHTO-U4ItnLU078eVkF1ZUZgYTy7THlTW3KSgFA2}
#' @keywords datasets
NULL

## ** gastricbypassL
#' @title Data From The Gastric Bypass Study (Long Format)
#' @name gastricbypassL
#' @rdname data-gastricbypassL
#'
#' @description  Data from the gastric bypass study
#' where the bodyweight and serum glucagon (a gut hormone) were measured in 20 obese subjects prior and after gastric bypass surgery.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{visit}: the visit index (factor).
#' \item \code{time}: the week at which the visit took place (numeric).
#' \item \code{weight}: bodyweight (in kg) measured during the visit.
#' \item \code{glucagonAUC}: glucagon measured during the visit (in pmol/l x hours).
#' }
#' 
#' @docType data
#' @usage data(gastricbypassL)
#' @references The effect of Roux-en-Y gastric bypass surgery on the gut mucosal gene expression profile and circulating gut hormones. \url{https://easddistribute.m-anage.com/from.storage?image=4iBH9mRQm1kfeEHULC2CxovdlyCtA1EHeVDdoffnZrAUGG9SHTO-U4ItnLU078eVkF1ZUZgYTy7THlTW3KSgFA2}
#' @keywords datasets
NULL

## * ncgs
## ** ncgsW
#' @title Data From National Cooperative Gallstone Study (Wide Format)
#' @name ncgsW
#' @rdname data-ncgsW
#'
#' @description  Data from the National Cooperative Gallstone Study (NCGS),
#' a randomized study where the level of serum cholesterol was measured at baseline and after intake of high-dose chenondiol (750mg/day) or placebo.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item \code{group}: treatment group (highdose or placebo).
#' \item \code{id}: patient identifier.
#' \item \code{cholest1}: cholesterol measurement at baseline (before treatment).
#' \item \code{cholest2}: cholesterol measurement at 6 months (after treatment).
#' \item \code{cholest3}: cholesterol measurement at 12 months (after treatment).
#' \item \code{cholest4}: cholesterol measurement at 20 months (after treatment).
#' \item \code{cholest5}: cholesterol measurement at 24 months (after treatment).
#' }
#' 
#' @docType data
#' @usage data(ncgsW)
#' @references Grundy SM, Lan SP, Lachin J. The effects of chenodiol on biliary lipids and their association with gallstone dissolution in the National Cooperative Gallstone Study (NCGS). J Clin Invest. 1984 Apr;73(4):1156-66. doi: 10.1172/JCI111301.  
#' @keywords datasets
NULL

## ** ncgsL
#' @title Data From National Cooperative Gallstone Study (Long Format)
#' @name ncgsL
#' @rdname data-ncgsL
#'
#' @description  Data from the National Cooperative Gallstone Study (NCGS),
#' a randomized study where the level of serum cholesterol was measured at baseline and after intake of high-dose chenondiol (750mg/day) or placebo.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item \code{group}: treatment group (highdose or placebo).
#' \item \code{id}: patient identifier.
#' \item \code{visit}: visit index.
#' \item \code{cholest}: cholesterol measurement.
#' \item \code{time}: time after the start of the study at which the measurement has been done (in month). Treatment is given at 0+.
#' }
#' 
#' @docType data
#' @usage data(ncgsL)
#' @references Grundy SM, Lan SP, Lachin J. The effects of chenodiol on biliary lipids and their association with gallstone dissolution in the National Cooperative Gallstone Study (NCGS). J Clin Invest. 1984 Apr;73(4):1156-66. doi: 10.1172/JCI111301.  
#' @keywords datasets
NULL

## * onycholysis
## ** onycholysisW
#' @title Data From The toenail onycholysis Study (Wide Format)
#' @name onycholysisW
#' @rdname data-onycholysisW
#'
#' @description  Data from the toenail onychomycosis study,
#' a randomized double-blind comparative study comparing Terbinafine (250 mg/day) against Itraconazole (200 mg/day) over 7 weeks.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{group}: treatment arm to which the patient has been randomized.
#' \item \code{response1},...,\code{response7}: degree of onycholysis (separation of the nail plate from the nail-bed). Can be 0 for none or mild, 1 for moderate or severe.
#' \item \code{time1},...,\code{time7}: time elapsed between baseline and the measurement (in months).
#' }
#' 
#' @docType data
#' @usage data(onycholysisW)
#' @references De Backer et al. Twelve weeks of continuous oral therapy for toenail onychomycosis caused by dermatophytes: A double-blind comparative trial of terbinafine 250 mg/day versus itraconazole 200 mg/day,  Journal of the American Academy of Dermatology, (2020) 110. doi: 10.1016/s0190-9622(98)70486-4
#' @keywords datasets
NULL


## ** onycholysisL
#' @title Data From The toenail onycholysis Study (Long Format)
#' @name onycholysisL
#' @rdname data-onycholysisL
#'
#' @description  Data from the toenail onychomycosis study,
#' a randomized double-blind comparative study comparing Terbinafine (250 mg/day) against Itraconazole (200 mg/day) over 7 weeks.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{group}: treatment arm to which the patient has been randomized (\code{"itraconazole"} or \code{"terbinafine"}).
#' \item \code{treatment}: treatment recieved by the patient at a given timepoint (\code{"none"}, \code{"itraconazole"}, or \code{"terbinafine"}).
#' \item \code{visit}: index of time at which the measurement was taken (\code{1} to \code{7}).
#' \item \code{time}: scheduled time, in months, for each measurement (\code{0}, \code{4}, \code{8}, \code{12}, \code{24}, \code{36}, \code{48}).
#' \item \code{obstime}: time, in months, at which the measurement was taken (numeric between \code{0} and \code{18.5}).
#' \item \code{response}: degree of onycholysis (separation of the nail plate from the nail-bed). Can be 0 for none or mild, 1 for moderate or severe.
#' }
#' 
#' @docType data
#' @usage data(onycholysisL)
#' @references De Backer et al. Twelve weeks of continuous oral therapy for toenail onychomycosis caused by dermatophytes: A double-blind comparative trial of terbinafine 250 mg/day versus itraconazole 200 mg/day,  Journal of the American Academy of Dermatology, (2020) 110. doi: 10.1016/s0190-9622(98)70486-4
#' @keywords datasets
NULL

## * potassium
## ** potassiumSingleW
#' @title Data From The Potassium Intake Study (Wide Format)
#' @name potassiumSingleW
#' @rdname data-potassiumSingleW
#'
#' @description  Data from the potassium intake study,
#' a randomized placebo-controlled crossover study where the effect of potassium supplement (90 mmol/day) on the renin-angiostensin-aldosteron system (RAAS) was assessed.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{sequence}: treatment group to which the patient has been randomized.
#' \item \code{treatment1}: treatment during the first time period.
#' \item \code{treatment2}: treatment during the second time period.
#' \item \code{auc1}: area under the curve of ?? during the first time period.
#' \item \code{auc2}: area under the curve of ?? during the second time period.
#' \item \code{bsauc1}: ??
#' \item \code{aldo1}: ??
#' \item \code{aldo2}: ??
#' }
#' 
#' @docType data
#' @usage data(potassiumSingleW)
#' @references Dreier et al. Effect of increased potassium intake on the reninangiotensinaldosterone system and subcutaneous resistance arteries: a randomized crossover study,
#' Nephrol Dial Transplant (1998) 110. doi: 10.1093/ndt/gfaa114
#' @keywords datasets
NULL

## ** potassiumSingleL
#' @title Data From The Potassium Intake Study (Long Format)
#' @name potassiumSingleL
#' @rdname data-potassiumSingleL
#'
#' @description  Data from the potassium intake study,
#' a randomized placebo-controlled crossover study where the effect of potassium supplement (90 mmol/day) on the renin-angiostensin-aldosteron system (RAAS) was assessed.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{sequence}: treatment group to which the patient has been randomized.
#' \item \code{period}: time period.
#' \item \code{treatment}: treatment during the time period.
#' \item \code{auc}: area under the curve of ?? during the time period.
#' \item \code{bsauc}: ??
#' \item \code{aldo}: ??
#' }
#' 
#' @docType data
#' @usage data(potassiumSingleL)
#' @references Dreier et al. Effect of increased potassium intake on the reninangiotensinaldosterone system and subcutaneous resistance arteries: a randomized crossover study,
#' Nephrol Dial Transplant (2020) 110. doi: 10.1093/ndt/gfaa114
#' @keywords datasets
NULL

## ** potassiumRepeatedL
#' @title Data From The Potassium Intake Study (Long Format with intermediate measurements)
#' @name potassiumRepeatedL
#' @rdname data-potassiumRepeatedL
#'
#' @description  Data from the potassium intake study,
#' a randomized placebo-controlled crossover study where the effect of potassium supplement (90 mmol/day) on the renin-angiostensin-aldosteron system (RAAS) was assessed.
#' This dataset is in the long format (i.e. one line per measurement) and contains measurement over 6 timepoints for each time period.
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{sequence}: treatment group to which the patient has been randomized.
#' \item \code{period}: time period.
#' \item \code{treatment}: treatment during the time period.
#' \item \code{time}: time within each period.
#' \item \code{aldo}: ??
#' }
#' 
#' @docType data
#' @usage data(potassiumRepeatedL)
#' @references Dreier et al. Effect of increased potassium intake on the reninangiotensinaldosterone system and subcutaneous resistance arteries: a randomized crossover study,
#' Nephrol Dial Transplant (2020) 110. doi: 10.1093/ndt/gfaa114
#' @keywords datasets
NULL

## * school
## ** schoolL
#' @title Simulated Data with 3-level struture (Long Format)
#' @name schoolL
#' @rdname data-schoolL
#'
#' @description Simulated data a nested structure: Student/Class/School and one outcome.
#'
#' \itemize{
#' \item \code{school}:
#' \item \code{class}: 
#' \item \code{student}:
#' \item \code{outcome}:  
#' }
#' 
#' @docType data
#' @usage data(schoolL)
#' @keywords datasets
NULL


## * sleep
## ** sleepL (embargo)
#' @title Data from a sleep study 
#' @name sleepL
#' @rdname data-sleepL
#' 
#' @description  
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{deprivation}: 
#' \item \code{condition}: 
#' \item \code{vigilance}: 
#' \item \code{signal.34} and \code{signal.98}: outcome.
#' }
#' 
#' @docType data
#' @usage data(sleepL)
#' @keywords datasets
NULL

## ## ** sleepL (embargo)
## #' @title MREG Data from a sleep deprivation experiment 
## #' @name sleepL
## #' @rdname data-sleepL
## #' 
## #' @description  Data from sleep deprivation cross-over study where all subjects have scan at baseline and two follow-up.
## #' They received placebo at one follow-up and carvedilol at the other.
## #' Sequence (Baseline-Placebo-Carvedilol or Baseline-Placebo-Carvedilol) is randomized and the current dataset do not contain the Carvedilol measurement.
## #' At each follow-up they under-go two scans: one awake and one asleep. This was not always possible leading to missing data.
## #' The signal measured relates to spectral characteristics of magnetic resonance encephalography (MREG) and has been modified compared to the real data.
## #'
## #' \itemize{
## #' \item \code{id}: patient identifier.
## #' \item \code{deprivation}: is the patient sleep deprived?
## #' \item \code{condition}: treatment received: none at baseline and then placebo or carvedilol.
## #' \item \code{vigilance}: is the patient awake or sleeping?
## #' \item \code{signal.34} and \code{signal.98}: modified MREG signal.
## #' }
## #' 
## #' @docType data
## #' @usage data(sleepL)
## #' @keywords datasets
## NULL

## * swabs
## ** swabsW
#' @title Data From The SWABS Study (Wide Format)
#' @name swabsW
#' @rdname data-swabsW
#'
#' @description Data from the swabs study,
#' where the pneumococcus was studied in 18 families with different space available for the household.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item \code{crowding}: space available in the household.
#' \item \code{family}: family serial number
#' \item \code{mother}: number of times the swab measurement was positive for the mother.
#' \item \code{father}: number of times the swab measurement was positive for the father.
#' \item \code{child1}: number of times the swab measurement was positive for the first child.
#' \item \code{child2}: number of times the swab measurement was positive for the second child.
#' \item \code{child3}: number of times the swab measurement was positive for the third child.
#' }
#' 
#' @docType data
#' @usage data(swabsW)
#' @references Grundy SM, Lan SP, Lachin J. The effects of chenodiol on biliary lipids and their association with gallstone dissolution in the National Cooperative Gallstone Study (SWABS). J Clin Invest. 1984 Apr;73(4):1156-66. doi: 10.1172/JCI111301.  
#' @keywords datasets
NULL

## ** swabsL
#' @title Data From The SWABS Study (Long Format)
#' @name swabsL
#' @rdname data-swabsL
#'
#' @description  Data from the swabs study,
#' where the pneumococcus was studied in 18 families with different space available for the household.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item \code{crowding}: space available in the household.
#' \item \code{family}: family serial number
#' \item \code{name}: type of family member.
#' \item \code{swabs}: number of times the swab measurement was positive.
#' }
#' 
#' @docType data
#' @usage data(swabsL)
#' @references TODO
#' @keywords datasets
NULL


## * vasscores
## ** vasscoresW
#' @title Data From The VAS Study (Wide Format)
#' @name vasscoresW
#' @rdname data-vasscoresW
#'
#' @description  Data from the VAS Study,
#' a randomized controlled clinial trial assessing the healing effect of topical zink sulfate on epidermal wound.
#' The study includes 30 heatlhy volunteers with induced wounds on each buttock which where subsequently treated with a different treatment for each wound.
#' Then the VAS-score (pain sensation on a 0-100mm visual analogue scale) was assessed after each treatment application and summarized by area under the curve.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{group}: treatment group to which the patient has been randomized.
#' \item \code{vasA}: VAS-score when using a zink shower gel.
#' \item \code{vasB}: VAS-score when using a placebo treatment (shower gel without zink).
#' \item \code{vasC}: VAS-score when using a control treatment with demineralized water.
#' }
#' 
#' @docType data
#' @usage data(vasscoresW)
#' @references TODO
#' @keywords datasets
NULL

## ** vasscoresL
#' @title Data From The VAS Study (Long Format)
#' @name vasscoresL
#' @rdname data-vasscoresL
#'
#' @description  Data from the VAS Study,
#' a randomized controlled clinial trial assessing the healing effect of topical zink sulfate on epidermal wound.
#' The study includes 30 heatlhy volunteers with induced wounds on each buttock which where subsequently treated with a different treatment for each wound.
#' Then the VAS-score (pain sensation on a 0-100mm visual analogue scale) was assessed after each treatment application and summarized by area under the curve.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item \code{id}: patient identifier.
#' \item \code{group}: treatment group to which the patient has been randomized.
#' \item \code{treat.num}:
#' \item \code{vas}: VAS-score relative to the wound.
#' \item \code{treatment}: Treatment used on the wound.
#' A: active treatment (zink shower gel),
#' B: placebo treatment (shower gel without zink),
#' C: control treatment (demineralized water).
#' }
#' 
#' @docType data
#' @usage data(vasscoresL)
#' @references TODO
#' @keywords datasets
NULL

## * vitamin
## ** vitaminW
#' @title Data From The Vitamin Study (Wide Format)
#' @name vitaminW
#' @rdname data-vitaminW
#'
#' @description  Data from the vitamin Study,
#' a randomized study where the growth of guinea pigs was monitored before and after intake of vitamin E/placebo.
#' The weight of each guinea pig was recorded at the end of week 1, 3, 4, 5, 6, and 7. Vitamin E/placebo is given at the beginning of week 5.
#' This dataset is in the wide format (i.e. one line per patient).
#'
#' \itemize{
#' \item \code{group}: treatment group: vitamin or placebo.
#' \item \code{animal}: identifier
#' \item \code{weight1}: weight (in g) of the pig at the end of week 1 (before treatment).
#' \item \code{weight3}: weight (in g) of the pig at the end of week 3 (before treatment).
#' \item \code{weight4}: weight (in g) of the pig at the end of week 4 (before treatment).
#' \item \code{weight5}: weight (in g) of the pig at the end of week 5 (after treatment).
#' \item \code{weight6}: weight (in g) of the pig at the end of week 6 (after treatment).
#' \item \code{weight7}: weight (in g) of the pig at the end of week 7 (after treatment).
#' }
#' 
#' @docType data
#' @usage data(vitaminW)
#' @references TODO
#' @keywords datasets
NULL

## ** vitaminL
#' @title Data From The Vitamin Study (Long Format)
#' @name vitaminL
#' @rdname data-vitaminL
#'
#' @description  Data from the vitamin Study,
#' a randomized study where the growth of guinea pigs was monitored before and after intake of vitamin E/placebo.
#' The weight of each guinea pig was recorded at the end of week 1, 3, 4, 5, 6, and 7. Vitamin E/placebo is given at the beginning of week 5.
#' This dataset is in the long format (i.e. one line per measurement).
#'
#' \itemize{
#' \item \code{group}: treatment group: vitamin or placebo.
#' \item \code{animal}: identifier
#' \item \code{weight1}: weight (in g) of the pig at the end of week 1 (before treatment).
#' \item \code{weight3}: weight (in g) of the pig at the end of week 3 (before treatment).
#' \item \code{weight4}: weight (in g) of the pig at the end of week 4 (before treatment).
#' \item \code{weight5}: weight (in g) of the pig at the end of week 5 (after treatment).
#' \item \code{weight6}: weight (in g) of the pig at the end of week 6 (after treatment).
#' \item \code{weight7}: weight (in g) of the pig at the end of week 7 (after treatment).
#' }
#' 
#' @docType data
#' @usage data(vitaminL)
#' @references Crowder and Hand (1990, p. 27) Analysis of Repeated Measures.
#' @keywords datasets
NULL


######################################################################
### doc-data.R ends here
