\name{get.valid.chrs}
\alias{get.valid.chrs}
\title{Annotation downloader helper}
\usage{
    get.valid.chrs(org)
}
\arguments{
    \item{org}{the organism for which to return the
    chromosomes.}
}
\value{
    A character vector of chromosomes.
}
\description{
    Returns a vector of chromosomes to maintain after
    annotation download. Internal use.
}
\examples{
hg18.chr <- get.valid.chrs("hg18")
}
\author{
    Panagiotis Moulos
}

