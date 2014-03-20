\name{get.exon.attributes}
\alias{get.exon.attributes}
\title{Annotation downloader helper}
\usage{
    get.exon.attributes()
}
\value{
    A character vector of Ensembl exon attributes.
}
\description{
    Returns a vector of genomic annotation attributes which
    are used by the biomaRt package in order to fetch the
    exon annotation for each organism. It has no parameters.
    Internal use.
}
\examples{
exon.attr <- get.exon.attributes()
}
\author{
    Panagiotis Moulos
}

