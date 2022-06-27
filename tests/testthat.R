## * load packages
library(testthat)
library(LMMstar)

library(data.table)
library(emmeans)
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
## setwd("c:/Users/hpl802/Documents/Github/LMMstar/tests/")
test_check("LMMstar")
