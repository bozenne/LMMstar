## * load packages
library(testthat)
library(LMMstar)

library(lava)
library(nlme)
library(lme4)
library(lmerTest)
library(emmeans)
library(numDeriv)
library(data.table)
library(multcomp)
library(reshape2)
library(qqtest)
library(ggplot2)

## * run tests
## setwd("~/Documents/GitHub/LMMstar/tests/")
## setwd("c:/Users/hpl802/Documents/Github/LMMstar/tests/")
test_check("LMMstar")
