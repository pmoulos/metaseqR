\name{reduce.gene.data}
\alias{reduce.gene.data}
\title{Reduce the gene annotation in case of not all chromosomes present in counts}
\usage{
    reduce.gene.data(exon.data, gene.data)
}
\arguments{
    \item{exon.data}{the exon annotation already reduced to
    the size of the input exon counts table.}

    \item{gene.data}{an annotation data frame from the same
    organism as \code{exon.counts} (such the ones produced by
    \code{get.annotation}).}
}
\value{
    The \code{gene.data} annotation, reduced to have the same
    chromosomes as in \code{exon.data}, or the original
    \code{gene.data} if \code{exon.data} do contain the
    standard chromosomes.
}
\description{
    This function reduces the gene annotation in case of exon
    reads and when the data to be analyzed do not contain all
    the standard chromosomes of the genome under
    investigation. This can greatly reduce processing time in
    these cases.
}
\examples{
\donttest{
data("hg19.exon.data",package="metaseqR")
gene.data <- get.annotation("hg19","gene","ensembl")
reduced.gene.data <- reduce.gene.data(hg19.exon.counts,
    gene.data)
}
}
\author{
    Panagiotis Moulos
}

