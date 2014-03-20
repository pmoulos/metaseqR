\name{filter.high}
\alias{filter.high}
\title{Filtering helper}
\usage{
    filter.high(x, f)
}
\arguments{
    \item{x}{a data numeric matrix.}

    \item{f}{a threshold.}
}
\description{
    High score filtering function. Internal use.
}
\examples{
data("mm9.gene.data",package="metaseqR")
counts <- as.matrix(mm9.gene.counts[,9:12])
f <- filter.low(counts,median(counts))
}
\author{
    Panagiotis Moulos
}

