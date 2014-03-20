\name{check.file.args}
\alias{check.file.args}
\title{File argument validator}
\usage{
    check.file.args(arg.name, arg.value)
}
\arguments{
    \item{arg.name}{argument name to display in a possible
    error.}

    \item{arg.value}{the filename to check.}
}
\description{
    Checks if a file exists for specific arguments requiring
    a file input. Internal use only.
}
\examples{
# OK
check.file.args("file",system.file("metaseqr_report.html",
    package="metaseqR"))
## Error!
#check.file.args("file",system.file("metaseqr_report.htm",
#    package="metaseqR"))
}
\author{
    Panagiotis Moulos
}
