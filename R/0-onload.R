### 0onload.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 16 2021 (11:59) 
## Version: 
## Last-Updated: May 27 2021 (17:31) 
##           By: Brice Ozenne
##     Update #: 7
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
    if (requireNamespace("emmeans", quietly = TRUE)){
        emmeans::.emm_register("lmm", "LMMstar")
    }
    LMMstar.options(reinitialise = TRUE) # generate .LMMstar-options when loading the package   
}

.getXlevels <- get(".getXlevels", envir = asNamespace("stats"), inherits = FALSE)

##----------------------------------------------------------------------
### 0onload.R ends here
