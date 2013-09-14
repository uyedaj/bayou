require(devtools)
require(roxygen2)
require(testthat)
load_all()
document()
test_dir("inst/tests/")

install("~/bayOU_1.0")

phytools.funs <- c("nodeHeights","plotSimmap","phenogram")
cat(paste0("importFrom(phytools, ", paste(phytools.funs, collapse=", "), ")"),
    paste0("export(", paste(phytools.funs, collapse=", "), ")"),
    file = "NAMESPACE",
    sep = "\n")

