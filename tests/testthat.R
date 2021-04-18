## * load packages
library(testthat)
library(LMMstar)

library(lava)
library(numDeriv)
library(data.table)

## * run tests
## setwd("~/Documents/GitHub/LMMstar/tests/")
test_check("LMMstar")
