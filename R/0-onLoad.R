'.onAttach' <- function(libname, pkgname="repeated") {
    desc <- utils::packageDescription(pkgname)
    packageStartupMessage(desc$Package, " version ",desc$Version)
}


