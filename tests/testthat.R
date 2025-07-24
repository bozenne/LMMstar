## * load packages
library(testthat)
library(LMMstar)

library(data.table)
library(ggplot2)
library(lava)
library(lme4)
library(lmerTest)
library(mice)
library(multcomp)
library(nlme)
library(numDeriv)
library(qqtest)

## * run tests
test_check("LMMstar")
## test_check("LMMstar", filter = "anova")
