# Package environment to store a couple of variables that must be global
meta.env <- new.env(parent=emptyenv())
assign("VERBOSE",NULL,envir=meta.env)
assign("LOGGER",NULL,envir=meta.env)
#.onLoad <- function(...) {
#    packageStartupMessage("\n\nWelcome to metaseqR 1.5.3. Please do not ",
#        "include \"bayseq\" in the statistics\nargument in metaseqr calls. ",
#        "There is an incompatibility with the latest\nversion of baySeq which ",
#        "we are fixing right now...\n\n")
#}
