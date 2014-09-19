\name{load.bs.genome}
\alias{load.bs.genome}
\title{Loads (or downloads) the required BSGenome package}
\usage{
    load.bs.genome(org)
}
\arguments{
    \item{org}{one of \code{\link{metaseqr}} supported 
    organisms.}
}
\value{
    A proper BSgenome package name.
}
\description{
    Retrieves the required BSgenome package when the 
    annotation source is \code{"ucsc"} or 
    \code{"refseq"}. These packages are required in order 
    to estimate the GC-content of the retrieved genes from 
    UCSC or RefSeq.
}
\examples{
\donttest{
bs.obj <- load.bs.genome("mm9")
}
}
\author{
    Panagiotis Moulos
}

