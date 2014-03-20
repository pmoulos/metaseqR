# Package environment to store a couple of variables that must be global
meta.env <- new.env(parent=emptyenv())
assign("VERBOSE",NULL,envir=meta.env)
assign("LOGGER",NULL,envir=meta.env)
