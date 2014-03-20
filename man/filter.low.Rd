\name{filter.low}
\alias{filter.low}
\title{Filtering helper}
\usage{
    filter.low(x, f)
}
\arguments{
    \item{x}{a data numeric matrix.}

    \item{f}{a threshold.}
}
\description{
    Low score filtering function. Internal use.
}
\examples{
data("mm9.gene.data",package="metaseqR")
counts <- as.matrix(mm9.gene.counts[,9:12])
f <- filter.low(counts,median(counts))
}
\author{
    Panagiotis Moulos
}

