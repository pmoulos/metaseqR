\name{make.html.rows}
\alias{make.html.rows}
\title{HTML report helper}
\usage{
    make.html.rows(mat)
}
\arguments{
    \item{mat}{the data matrix, usually the output of
    \code{\link{make.html.cells}} function.}
}
\value{
    A character vector with html formatted rows of a matrix.
}
\description{
    Returns a character vector with html formatted rows.
    Essentially, it collapses every row of a matrix to a
    single character and puts a <tr></tr> tag set around. It
    is meant to be applied to the output of
    \code{\link{make.html.cells}}. Internal use.
}
\examples{
\donttest{
data.matrix <- round(1000*matrix(runif(400),100,4))
rownames(data.matrix) <- paste("gene_",1:100,sep="")
colnames(data.matrix) <- c("C1","C2","T1","T2")
the.cells <- make.html.cells(data.matrix)
the.rows <- make.html.rows(the.cells)
}
}
\author{
    Panagiotis Moulos
}

