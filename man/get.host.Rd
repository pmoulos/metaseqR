\name{get.host}
\alias{get.host}
\title{Annotation downloader helper}
\usage{
    get.host(org)
}
\arguments{
    \item{org}{the organism for which to return the host
    address.}
}
\value{
    A string with the host address.
}
\description{
    Returns the appropriate Ensembl host address to get
    different versions of annotation from. Internal use.
}
\examples{
mm9.hist <- get.host("mm9")
}
\author{
    Panagiotis Moulos
}

