\name{get.ucsc.organism}
\alias{get.ucsc.organism}
\title{Return a proper formatted organism alias}
\usage{
    get.ucsc.organism(org)
}
\arguments{
    \item{org}{one of the \code{\link{metaseqr}}
    supported organism.}
}
\value{
    A proper organism alias.
}
\description{
    Returns the proper UCSC Genome Browser database 
    organism alias based on what is given to metaseqR.
    Internal use.
}
\examples{
\donttest{
org <- get.ucsc.organism("danrer7")
}
}
\author{
    Panagiotis Moulos
}

