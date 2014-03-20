\name{make.html.cells}
\alias{make.html.cells}
\title{HTML report helper}
\usage{
    make.html.cells(mat, type = "numeric", digits = 3)
}
\arguments{
    \item{mat}{the data matrix (numeric or character)}

    \item{type}{the type of data in the matrix
    (\code{"numeric"} or \code{"character"})}

    \item{digits}{the number of digits on the right of the
    decimal points to pass to \code{\link{formatC}}. It has
    meaning when \code{type="numeric"}.}
}
\value{
    A character matrix with html formatted cells.
}
\description{
    Returns a character matrix with html formatted table
    cells. Essentially, it converts the input data to text
    and places them in a <td></td> tag set. Internal use.
}
\examples{
\donttest{
data.matrix <- round(1000*matrix(runif(400),100,4))
rownames(data.matrix) <- paste("gene_",1:100,sep="")
colnames(data.matrix) <- c("C1","C2","T1","T2")
the.cells <- make.html.cells(data.matrix)
}
}
\author{
    Panagiotis Moulos
}

