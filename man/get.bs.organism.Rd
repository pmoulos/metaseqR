\name{get.bs.organism}
\alias{get.bs.organism}
\title{Return a proper formatted BSgenome organism name}
\usage{
    get.bs.organism(org)
}
\arguments{
    \item{org}{one of metaseqR supported organisms.}
}
\value{
    A proper BSgenome package name.
}
\description{
    Returns a properly formatted BSgenome package name 
    according to metaseqR's supported organism. Internal 
    use.
}
\examples{
\donttest{
bs.name <- get.bs.organism("hg18")
}
}
\author{
    Panagiotis Moulos
}

