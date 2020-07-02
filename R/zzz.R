# zzz.R
# ::gamsel::


.onAttach <- function(libname, pkgname) {
  
  gamsel.version <- packageVersion(pkgname)
  packageStartupMessage(paste0("  ~ ", pkgname, " ", gamsel.version,
                               ": Welcome, ", Sys.getenv("USER")))
}
