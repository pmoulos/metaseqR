\name{stat.noiseq}
\alias{stat.noiseq}
\title{Statistical testing with NOISeq}
\usage{
    stat.noiseq(object, sample.list, contrast.list = NULL,
        stat.args = NULL, gene.data = NULL, log.offset = 1)
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
    \code{get.defaults("statistics",} \code{"noiseq")} for an
    example and how you can modify it.}

    \item{gene.data}{an optional annotation data frame (such
    the ones produced by \code{get.annotation} which contains
    the GC content for each gene and from which the gene
    lengths can be inferred by chromosome coordinates.}

    \item{log.offset}{a number to be added to each element of
    data matrix in order to avoid Infinity on log type data
    transformations.}
}
\value{
    A named list of NOISeq q-values, whose names are the
    names of the contrasts.
}
\description{
    This function is a wrapper over NOISeq statistical
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
lengths <- round(1000*runif(nrow(data.matrix)))
starts <- round(1000*runif(nrow(data.matrix)))
ends <- starts + lengths
gc=runif(nrow(data.matrix))
gene.data <- data.frame(
    chromosome=c(rep("chr1",nrow(data.matrix)/2),
        rep("chr2",nrow(data.matrix)/2)),
    start=starts,end=ends,gene_id=rownames(data.matrix),gc_content=gc
)
norm.data.matrix <- normalize.noiseq(data.matrix,sample.list,gene.data)
p <- stat.noiseq(norm.data.matrix,sample.list,contrast,
    gene.data=gene.data)
}
}
\author{
    Panagiotis Moulos
}

