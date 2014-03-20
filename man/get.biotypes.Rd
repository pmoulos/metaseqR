\name{get.biotypes}
\alias{get.biotypes}
\title{Biotype converter}
\usage{
    get.biotypes(a)
}
\arguments{
    \item{a}{the annotation data frame (output of
    \code{\link{get.annotation}}).}
}
\value{
    A character vector of biotypes.
}
\description{
    Returns biotypes as character vector. Internal use.
}
\examples{
\donttest{
hg18.genes <- get.annotation("hg18","gene")
hg18.bt <- get.biotypes(hg18.genes)
}
}
\author{
    Panagiotis Moulos
}

