\name{make.venn.areas}
\alias{make.venn.areas}
\title{Helper for Venn diagrams}
\usage{
    make.venn.areas(n)
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
    VennDiagram and list members the internal encoding
    (uppercase letters A to E and combinations among then)
    used to encode the pairwise comparisons to create the
    intersections needed for the Venn diagrams. Internal use
    mostly.
}
\examples{
\donttest{
sets <- c("apple","pear","banana")
pairs <- make.venn.pairs(sets)
areas <- make.venn.areas(length(sets))
}
}
\author{
    Panagiotis Moulos
}

