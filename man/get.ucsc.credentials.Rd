\name{get.ucsc.credentials}
\alias{get.ucsc.credentials}
\title{Return host, username and password for UCSC 
    Genome Browser database}
\usage{
    get.ucsc.credentials()
}
\value{
    A named character vector.
}
\description{
    Returns a character vector with a hostname, username 
    and password to connect to the UCSC Genome Browser 
    database to retrieve annotation. Internal use mostly.
}
\examples{
\donttest{
db.creds <- get.ucsc.credentials()
}
}
\author{
    Panagiotis Moulos
}

