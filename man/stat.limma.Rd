\name{stat.limma}
\alias{stat.limma}
\title{Statistical testing with limma}
\usage{
    stat.limma(object, sample.list, contrast.list = NULL,
        stat.args = NULL)
}
\arguments{
    \item{object}{a matrix or an object specific to each
    normalization algorithm supported by metaseqR, containing
    normalized counts. Apart from matrix (also for NOISeq),
    the object can be a SeqExpressionSet (EDASeq),
    CountDataSet (DESeq) or DGEList (edgeR).}

    \item{sample.list}{the list containing condition names
    and the samples under each condition.}

    \item{contrast.list}{a named structured list of contrasts
    as returned by \code{\link{make.contrast.list}} or just
    the vector of contrasts as defined in the main help page
    of \code{\link{metaseqr}}.}

    \item{stat.args}{a list of edgeR statistical algorithm
    parameters. See the result of
    \code{get.defaults("statistics",} \code{"limma")} for an
    example and how you can modify it.}
}
\value{
    A named list of p-values, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over limma statistical
    testing. It accepts a matrix of normalized gene counts or
    an S4 object specific to each normalization algorithm
    supported by metaseqR.
}
\examples{
\donttest{
require(DESeq)
data.matrix <- counts(makeExampleCountDataSet())
sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
norm.data.matrix <- normalize.edger(data.matrix,sample.list)
p <- stat.limma(norm.data.matrix,sample.list,contrast)
}
}
\author{
    Panagiotis Moulos
}

