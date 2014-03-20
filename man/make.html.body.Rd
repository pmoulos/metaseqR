\name{make.html.body}
\alias{make.html.body}
\title{HTML report helper}
\usage{
    make.html.body(mat)
}
\arguments{
    \item{mat}{the character vector produced by
    \code{\link{make.html.rows}}.}
}
\value{
    A character vector with the body of mat formatted in
    html.
}
\description{
    Returns a character vector with an html formatted table.
    Essentially, it collapses the input rows to a single
    character and puts a <tbody></tbody> tag set around. It
    is meant to be applied to the output of
    \code{\link{make.html.rows}}. Internal use.
}
\examples{
\donttest{
data.matrix <- round(1000*matrix(runif(400),100,4))
rownames(data.matrix) <- paste("gene_",1:100,sep="")
colnames(data.matrix) <- c("C1","C2","T1","T2")
the.cells <- make.html.cells(data.matrix)
the.header <- make.html.header(the.cells[1,])
the.rows <- make.html.rows(the.cells)
the.body <- make.html.body(the.rows)
}
}
\author{
    Panagiotis Moulos
}

