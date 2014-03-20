\name{normalize.deseq}
\alias{normalize.deseq}
\title{Normalization based on the DESeq package}
\usage{
    normalize.deseq(gene.counts, sample.list,
        norm.args = NULL, output = c("matrix", "native"))
}
\arguments{
    \item{gene.counts}{a table where each row represents a
    gene and each column a sample. Each cell contains the
    read counts for each gene and sample. Such a table can be
    produced outside metaseqr and is imported during the
    basic metaseqr workflow.}

    \item{sample.list}{the list containing condition names
    and the samples under each condition.}

    \item{norm.args}{a list of DESeq normalization
    parameters. See the result of
    \code{get.defaults("normalization",} \code{"deseq")} for
    an example and how you can modify it.}

    \item{output}{the class of the output object. It can be
    \code{"matrix"} (default) for versatility with other
    tools or \code{"native"} for the DESeq native S4 object
    (CountDataSet). In the latter case it should be handled
    with suitable DESeq methods.}
}
\value{
    A matrix or a CountDataSet with normalized counts.
}
\description{
    This function is a wrapper over DESeq normalization. It
    accepts a matrix of gene counts (e.g. produced by
    importing an externally generated table of counts to the
    main metaseqr pipeline).
}
\examples{
\donttest{
require(DESeq)
data.matrix <- counts(makeExampleCountDataSet())
sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplot.boxplot(data.matrix,sample.list)

norm.data.matrix <- normalize.deseq(data.matrix,sample.list)
diagplot.boxplot(norm.data.matrix,sample.list)
}
}
\author{
    Panagiotis Moulos
}

