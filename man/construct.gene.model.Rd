\name{construct.gene.model}
\alias{construct.gene.model}
\title{Assemble a gene model based on exon counts}
\usage{
    construct.gene.model(exon.counts, sample.list, gene.data,
        multic = FALSE)
}
\arguments{
    \item{exon.counts}{the exon counts data frame produced by
    reading the exon read counts file.}

    \item{sample.list}{the list containing condition names
    and the samples under each condition.}

    \item{gene.data}{an annotation data frame from the same
    organism as \code{exon.counts} (such the ones produced by
    \code{get.annotation}).}

    \item{multic}{a logical value indicating the presence of
    multiple cores. Defaults to \code{FALSE}. Do not change
    it if you are not sure whether package multicore has been
    loaded or not.}
}
\value{
    A named list where names represent samples. Each list
    member is a also a named list where names correspond to
    gene ids and members are named vectors. Each vector is
    named according to the exons corresponding to each gene
    and contains the read counts for each exon. This
    structure is used for exon filtering and assembling final
    gene counts in the metaseqr pipeline.
}
\description{
    This function assembles gene models (single genes, not
    isoforms) based on the input exon read counts file and a
    gene annotation data frame, either from an external file
    provided by the user, or with the
    \code{\link{get.annotation}} function. The
    \code{gene.data} argument should have a specific format
    and for this reason it's better to use one of the two
    aforementioned ways to supply it. This function is
    intended mostly for internal use but can be used if the
    requirements are met.
}
\examples{
\donttest{
# Takes some time to run...
data("hg19.exon.data",package="metaseqR")
gene.data <- get.annotation("hg19","gene")
reduced.gene.data <- reduce.gene.data(hg19.exon.counts,
    gene.data)
multic <- check.parallel(0.4)
gene.model <- construct.gene.model(hg19.exon.counts,
    sample.list.hg19,gene.data,multic)
}
}
\author{
    Panagiotis Moulos
}

