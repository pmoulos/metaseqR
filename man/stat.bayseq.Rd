\name{stat.bayseq}
\alias{stat.bayseq}
\title{Statistical testing with baySeq}
\usage{
    stat.bayseq(object, sample.list, contrast.list = NULL,
        stat.args = NULL, libsize.list = NULL)
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
    \code{get.defaults("statistics",} \code{"bayseq")} for an
    example and how you can modify it.}

    \item{libsize.list}{an optional named list where names
    represent samples (MUST be the same as the samples in
    \code{sample.list}) and members are the library sizes
    (the sequencing depth) for each sample. If not provided,
    they will be estimated from baySeq.}
}
\value{
    A named list of the value 1-likelihood that a gene is
    differentially expressed, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over baySeq statistical
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
norm.data.matrix <- normalize.edaseq(data.matrix,sample.list,gene.data)
p <- stat.bayseq(norm.data.matrix,sample.list,contrast)
}
}
\author{
    Panagiotis Moulos
}

