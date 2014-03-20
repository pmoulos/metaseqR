\name{check.parallel}
\alias{check.parallel}
\title{Parallel run validator}
\usage{
    check.parallel(rc)
}
\arguments{
    \item{rc}{fraction of available cores to use.}
}
\description{
    Checks existence of multiple cores and loads multicore
    package.
}
\examples{
multic <- check.parallel(0.8)
}
\author{
    Panagiotis Moulos
}

