\name{get.dataset}
\alias{get.dataset}
\title{Annotation downloader helper}
\usage{
    get.dataset(org)
}
\arguments{
    \item{org}{the organism for which to return the
    identifier.}
}
\value{
    A string with the dataset identifier.
}
\description{
    Returns a dataset (gene or exon) identifier for each
    organism recognized by the Biomart service for Ensembl.
    Internal use.
}
\examples{
\donttest{
dm3.id <- get.dataset("dm3")
}
}
\author{
    Panagiotis Moulos
}

