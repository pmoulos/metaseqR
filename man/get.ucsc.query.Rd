\name{get.ucsc.query}
\alias{get.ucsc.query}
\title{Return queries for the UCSC Genome Browser database,
    according to organism and source}
\usage{
    get.ucsc.query(org, type, refdb="ucsc")
}
\arguments{
    \item{org}{one of metaseqR supported organisms.}    

    \item{type}{either \code{"gene"} or \code{"exon"}.}

    \item{refdb}{one of \code{"ucsc"} or \code{"refseq"} 
    to use the UCSC or RefSeq annotation sources 
    respectively.}
}
\value{
    A valid SQL query.
}
\description{
    Returns an SQL query to be used with a connection to 
    the UCSC Genome Browser database and fetch metaseqR 
    supported organism annotations. This query is constructed 
    based on the data source and data type to be returned.
}
\examples{
\donttest{
db.query <- get.ucsc.query("gene","ucsc")
}
}
\author{
    Panagiotis Moulos
}

