\name{get.ucsc.tabledef}
\alias{get.ucsc.tabledef}
\title{Get SQLite UCSC table defintions, according to organism and source}
\usage{
    get.ucsc.tabledef(org, type, refdb="ucsc", what="queries")
}
\arguments{
    \item{org}{one of metaseqR supported organisms.}    

    \item{type}{either \code{"gene"} or \code{"exon"}.}

    \item{refdb}{one of \code{"ucsc"} or \code{"refseq"} 
    to use the UCSC or RefSeq annotation sources 
    respectively.}
    
    \item{what}{either \code{"queries"} for SQLite table definitions
    or \code{"fields"} for only a vector of table field names.}
}
\value{
    A list with SQLite table definitions.
}
\description{
    Creates a list of UCSC Genome Browser database tables and 
    their SQLite definitions with the purpose of creating a 
    temporary SQLite database to be used used with metaseqR. 
    This functionality is used when the package RMySQL is not 
    available for some reason, e.g. Windows machines.
}
\examples{
\donttest{
db.tabledefs <- get.ucsc.tabledef("hg18","gene","ucsc")
}
}
\author{
    Panagiotis Moulos
}
