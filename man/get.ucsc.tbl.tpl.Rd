\name{get.ucsc.tbl.tpl}
\alias{get.ucsc.tbl.tpl}
\title{Create SQLite UCSC table template defintions}
\usage{
    get.ucsc.tbl.tpl(tab, what="queries")
}
\arguments{
    \item{tab}{name of UCSC database table.}
    
    \item{what}{\code{"queries"} for SQLite table 
    definitions or \code{"fields"} for table column names.}
}
\value{
    An SQLite table definition.
}
\description{
    Returns an SQLIte table template defintion, 
    according to  UCSC Genome Browser database 
    table schemas. This functionality is used 
    when the package RMySQL is not available 
    for some reason, e.g. Windows machines. 
    Internal use only.
}
\examples{
\donttest{
db.table.tmpl <- get.ucsc.tbl.tpl("knownCanonical")
}
}
\author{
    Panagiotis Moulos
}

