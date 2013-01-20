.onLoad <- function(libname, pkgname) {
  library.dynam(pkgname, pkgname, NULL)
#  --> error on R-Forge
#  but
#http://heuristically.wordpress.com/2009/12/24/
#                  error-onload-failed-loadnamespace-rweka/
# says
#  Solution: Just close R and re-open it.
#  Cause: Apparently the installation requires some initialization.
# *****
#  LibLoc <- system.file(package=pkgname)
#  library.dynam(pkgname, pkgname, LibLoc)
#   --> error on my Windows 7
}
