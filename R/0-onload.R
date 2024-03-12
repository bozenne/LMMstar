### 0onload.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 16 2021 (11:59) 
## Version: 
## Last-Updated: mar 12 2024 (09:55) 
##           By: Brice Ozenne
##     Update #: 8
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
LMMstar.env <- new.env() # create a specific environment for the package

.onAttach <- function(lib, pkg="LMMstar") {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(desc$Package, " version ",desc$Version)
    LMMstar.options(reinitialise = TRUE) # generate .LMMstar-options when loading the package   
}

.getXlevels <- get(".getXlevels", envir = asNamespace("stats"), inherits = FALSE)

##----------------------------------------------------------------------
### 0onload.R ends here
