\name{make.venn.pairs}
\alias{make.venn.pairs}
\title{Helper for Venn diagrams}
\usage{
    make.venn.pairs(algs)
}
\arguments{
    \item{algs}{a vector with the names of the sets (up to
    length 5, if larger, it will be truncated with a
    warning).}
}
\value{
    A list with as many pairs as the comparisons to be made
    for the construction of the Venn diagram. The pairs are
    encoded with the uppercase letters A through E, each one
    corresponding to order of the input sets.
}
\description{
    This function creates a list of pairwise comparisons to
    be performed in order to create an up to 5-way Venn
    diagram using the R package VennDiagram. Internal use
    mostly.
}
\examples{
\donttest{
sets <- c("apple","pear","banana")
pairs <- make.venn.pairs(sets)
}
}
\author{
    Panagiotis Moulos
}

