### 0onload.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 16 2021 (11:59) 
## Version: 
## Last-Updated: Apr 16 2021 (12:15) 
##           By: Brice Ozenne
##     Update #: 5
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


##----------------------------------------------------------------------
### 0onload.R ends here
