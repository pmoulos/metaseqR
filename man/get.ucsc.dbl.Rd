\name{get.ucsc.dbl}
\alias{get.ucsc.dbl}
\title{Download annotation from UCSC servers, according to organism and 
    source}
\usage{
    get.ucsc.dbl(org, type, refdb="ucsc")
}
\arguments{
    \item{org}{one of metaseqR supported organisms.}    

    \item{type}{either \code{"gene"} or \code{"exon"}.}

    \item{refdb}{one of \code{"ucsc"} or \code{"refseq"} 
    to use the UCSC or RefSeq annotation sources 
    respectively.}
}
\value{
    An SQLite database.
}
\description{
    Directly downloads UCSC and RefSeq annotation files from UCSC servers 
    to be used with metaseqR. This functionality is used when the package 
    RMySQL is not available for some reason, e.g. Windows machines. It
    created an SQLite database where the same queries can be used.
}
\examples{
\donttest{
db.file <- get.ucsc.dbl("hg18","gene","ucsc")
}
}
\author{
    Panagiotis Moulos
}

