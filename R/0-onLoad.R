'.onAttach' <- function(libname, pkgname="handrem") {
    desc <- utils::packageDescription(pkgname)
    packageStartupMessage(desc$Package, " version ",desc$Version)
}


