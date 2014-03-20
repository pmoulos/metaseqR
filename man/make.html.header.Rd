\name{make.html.header}
\alias{make.html.header}
\title{HTML report helper}
\usage{
    make.html.header(h)
}
\arguments{
    \item{h}{the colnames of a matrix or data frame, usually
    as output of \code{\link{make.html.cells}} function.}
}
\value{
    A character vector with html formatted header of a
    matrix.
}
\description{
    Returns a character vector with an html formatted table
    head row. Essentially, it collapses the input row to a
    single character and puts a <th></th> tag set around. It
    is meant to be applied to the output of
    \code{\link{make.html.cells}}. Internal use.
}
\examples{
\donttest{
data.matrix <- round(1000*matrix(runif(400),100,4))
rownames(data.matrix) <- paste("gene_",1:100,sep="")
colnames(data.matrix) <- c("C1","C2","T1","T2")
the.cells <- make.html.cells(data.matrix)
the.header <- make.html.header(the.cells[1,])
}
}
\author{
    Panagiotis Moulos
}

