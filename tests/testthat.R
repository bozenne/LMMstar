## * load packages
library(testthat)
library(repeated)

library(lava)
library(data.table)

## * run tests
## setwd("~/Documents/GitHub/repeated/tests/")
test_check("repeated")

beta <- 0

m <- lvm()
latent(m) <- ~eta
distribution(m,"X1") <- lava::binomial.lvm(size = 1, p = 0.5)
regression(m, eta ~ X1) <- beta

set.seed(10)
sim(m,10)


m <- lvm()
m <- lava::`latent<-`(m,value = ~eta)
m <- lava::`distribution<-`(m,"X1", value = lava::binomial.lvm(size = 1, p = 0.5))
m <- lava::`regression<-`(m, eta ~ X1, value = beta)

set.seed(10)
sim(m,10)


