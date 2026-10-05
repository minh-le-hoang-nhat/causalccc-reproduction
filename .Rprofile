# Read automatically when R starts in this folder.
# Uses CRAN as it was on 2026-10-03 (Posit Package Manager snapshot), so everyone installs identical versions.
# On Linux, pre-compiled binaries for the local distribution are used (e.g. noble = Ubuntu 24.04).
local({
  snapshot <- "2026-10-03"
  base <- "https://packagemanager.posit.co/cran"
  codename <- character(0)
  if (file.exists("/etc/os-release")) {
    line <- grep("^VERSION_CODENAME=", readLines("/etc/os-release"), value = TRUE)
    codename <- sub("^VERSION_CODENAME=", "", line)
  }
  url <- if (length(codename) == 1 && nzchar(codename)) {
    paste(base, "__linux__", codename, snapshot, sep = "/")
  } else {
    paste(base, snapshot, sep = "/")
  }
  options(repos = c(CRAN = url))
})
