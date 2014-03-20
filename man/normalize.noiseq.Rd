\name{normalize.noiseq}
\alias{normalize.noiseq}
\title{Normalization based on the NOISeq package}
\usage{
    normalize.noiseq(gene.counts, sample.list,
        norm.args = NULL, gene.data = NULL, log.offset = 1,
        output = c("matrix", "native"))
}
\arguments{
    \item{gene.counts}{a table where each row represents a
    gene and each column a sample. Each cell contains the
    read counts for each gene and sample. Such a table can be
    produced outside metaseqr and is imported during the
    basic metaseqr workflow.}

    \item{sample.list}{the list containing condition names
    and the samples under each condition.}

    \item{norm.args}{a list of NOISeq normalization
    parameters. See the result of
    \code{get.defaults("normalization",} \code{"noiseq")} for
    an example and how you can modify it.}

    \item{gene.data}{an optional annotation data frame (such
    the ones produced by \code{get.annotation} which contains
    the GC content for each gene and from which the gene
    lengths can be inferred by chromosome coordinates.}

    \item{log.offset}{an offset to use to avoid infinity in
    logarithmic data transformations.}

    \item{output}{the class of the output object. It can be
    \code{"matrix"} (default) for versatility with other
    tools or \code{"native"} for the NOISeq native S4 object
    (SeqExpressionSet). In the latter case it should be
    handled with suitable NOISeq methods.}
}
\value{
    A matrix with normalized counts.
}
\description{
    This function is a wrapper over NOISeq normalization. It
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
diagplot.boxplot(norm.data.matrix,sample.list)
}
}
\author{
    Panagiotis Moulos
}

