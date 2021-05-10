## * load packages
library(testthat)
library(LMMstar)

library(lava)
library(numDeriv)
library(data.table)
library(multcomp)

## * run tests
## setwd("~/Documents/GitHub/LMMstar/tests/")
test_check("LMMstar")
