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
library(reshape2)
library(qqtest)

## * run tests
## setwd("~/Documents/GitHub/LMMstar/tests/")
test_check("LMMstar")
