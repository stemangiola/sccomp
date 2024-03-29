core <- c("tidyseurat", "tidySingleCellExperiment")

core_unloaded <- function() {
  search <- paste0("package:", core)
  core[!search %in% search()]
}

# Attach the package from the same library it was loaded from before.
# [source: https://github.com/tidy-biology/tidyverse/issues/171]
same_library <- function(pkg) {
  loc <- if (pkg %in% loadedNamespaces()) 
    dirname(getNamespaceInfo(pkg, "path"))
  library(pkg, lib.loc=loc, character.only=TRUE, warn.conflicts=FALSE)
}

tidyverse_attach <- function() {
  to_load <- core_unloaded()
  
  suppressPackageStartupMessages(
    lapply(to_load, same_library))
  
  invisible(to_load)
}