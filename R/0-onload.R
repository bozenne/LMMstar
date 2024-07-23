### 0onload.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 16 2021 (11:59) 
## Version: 
## Last-Updated: Jul 23 2024 (10:59) 
##           By: Brice Ozenne
##     Update #: 9
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
