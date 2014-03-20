\name{make.grid}
\alias{make.grid}
\title{Optimize rectangular grid plots}
\usage{
    make.grid(n)
}
\arguments{
    \item{n}{An integer, denoting the total number of plots
    to be created.}
}
\value{
    A 2-element vector with the dimensions of the grid.
}
\description{
    Returns a vector for an optimized m x m plot grid to be
    used with e.g. \code{par(mfrow)}. m x m is as close as
    possible to the input n. Of course, there will be empty
    grid positions if n < m x m.
}
\examples{
\donttest{
g1 <- make.grid(16) # Returns c(4,4)
g2 <- make.grid(11) # Returns c(4,3)
}
}
\author{
    Panagiotis Moulos
}

