\name{make.html.table}
\alias{make.html.table}
\title{HTML report helper}
\usage{
    make.html.table(b, h = NULL, id = NULL)
}
\arguments{
    \item{b}{the table body as produced by
    \code{\link{make.html.body}}.}

    \item{h}{the table header as produced by
    \code{\link{make.html.header}}.}

    \item{id}{the table id attribute.}
}
\value{
    A fully formatted html table.
}
\description{
    Returns a character vector with a fully html formatted
    table. Essentially, it binds the outputs of
    \code{\link{make.html.cells}},
    \code{\link{make.html.rows}},
    \code{\link{make.html.header}} and
    \code{\link{make.html.body}} to the final table and
    optionally assigns an id attribute. The above functions
    are meant to format a data table so as it can be rendered
    by external tools such as DataTables.js during a report
    creation. It is meant for internal use.
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
the.table <- make.html.table(the.body,the.header,id="my_table")
}
}
\author{
    Panagiotis Moulos
}

