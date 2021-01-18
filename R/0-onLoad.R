'.onAttach' <- function(libname, pkgname="LMMstar") {
    desc <- utils::packageDescription(pkgname)
    packageStartupMessage(desc$Package, " version ",desc$Version)
}


