\name{check.main.args}
\alias{check.main.args}
\title{Main argument validator}
\usage{
    check.main.args(main.args)
}
\arguments{
    \item{main.args}{ a list of parameters 
    with which metaseqr is called (essentially,
    the output of \code{\link{match.call}}.}
}
\description{
    Checks if the arguments passed to 
    \code{\link{metaseqr}} are valid and throws a 
    warning about the invalid ones (which are 
    ignored anyway because of the \code{...} in 
    \code{\link{metaseqr}}. However, for this reason 
    this function is useful as some important 
    parameter faults might go unnoticed in the 
    beginning and cause a failure afterwards.
}
\author{
    Panagiotis Moulos
}

