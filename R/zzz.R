# zzz.R
# ::gamsel::

gamsel.version <- packageVersion("gamsel")

.onAttach <- function(libname, pkgname) {
  
  packageStartupMessage(paste0("  ~ ", pkgname, " ", gamsel.version,
                               ": Welcome, ", Sys.getenv("USER")))
}
