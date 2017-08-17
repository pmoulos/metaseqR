\name{get.transcript.utr.attributes}
\alias{get.transcript.utr.attributes}
\title{Annotation downloader helper}
\usage{
    get.transcript.utr.attributes(org)
}
\arguments{
    \item{org}{one of the supported organisms.}
}
\value{
    A character vector of Ensembl transcript attributes.
}
\description{
    Returns a vector of genomic annotation attributes which
    are used by the biomaRt package in order to fetch the
    transcript annotation for each organism. It has no 
    parameters. Internal use.
}
\examples{
tran.attr <- get.transcript.utr.attributes("mm9")
}
\author{
    Panagiotis Moulos
}
