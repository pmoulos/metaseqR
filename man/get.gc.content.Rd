\name{get.gc.content}
\alias{get.gc.content}
\title{Return a named vector of GC-content for each genomic
    region}
\usage{
    get.gc.content(ann, org)
}
\arguments{
    \item{ann}{a data frame which can be converted to a 
    GRanges object, that means it has at least the 
    chromosome, start, end fields. Preferably, the output 
    of \code{link{get.ucsc.annotation}}.}

    \item{org}{one of metaseqR supported organisms.}
}
\value{
    A named numeric vector.
}
\description{
    Returns a named numeric vector (names are the 
    genomic region names, e.g. genes) given a data 
    frame which can be converted to a GRanges object 
    (e.g. it has at least chromosome, start, end 
    fields). This function works best when the input
    annotation data frame has been retrieved using one 
    of the SQL queries generated from 
    \code{\link{get.ucsc.query}}, used in 
    \code{\link{get.ucsc.annotation}}.
}
\examples{
\donttest{
ann <- get.ucsc.annotation("mm9","gene","ucsc")
gc <- get.gc.content(ann,"mm9")
}
}
\author{
    Panagiotis Moulos
}

