\name{make.venn.counts}
\alias{make.venn.counts}
\title{Helper for Venn diagrams}
\usage{
    make.venn.counts(n)
}
\arguments{
    \item{n}{the number of the sets used for the Venn
    diagram.}
}
\value{
    A named list, see descritpion.
}
\description{
    This function creates a list with names the arguments of
    the Venn diagram construction functions of the R package
    VennDiagram and list members are initially \code{NULL}.
    They are filled by the \code{\link{diagplot.venn}}
    function. Internal use mostly.
}
\examples{
\donttest{
sets <- c("apple","pear","banana")
counts <- make.venn.counts(length(sets))
}
}
\author{
    Panagiotis Moulos
}

