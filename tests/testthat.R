## * load packages
library(testthat)
library(LMMstar)

library(lava)
library(nlme)
library(emmeans)
library(numDeriv)
library(data.table)
library(multcomp)
library(qqtest)
library(ggplot2)

## * run tests
## setwd("~/Documents/GitHub/LMMstar/tests/")
## setwd("c:/Users/hpl802/Documents/Github/LMMstar/tests/")
test_check("LMMstar")
