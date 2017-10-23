\name{filter.genes}
\alias{filter.genes}
\title{Filter gene expression based on gene counts}
\usage{
    filter.genes(gene.counts, gene.data, gene.filters,
    sample.list)
}
\arguments{
    \item{gene.counts}{a matrix of gene counts, preferably
    after the normalization procedure.}

    \item{gene.data}{an annotation data frame usually
    obtained with \code{\link{get.annotation}} containing the
    unique gene accession identifiers.}

    \item{gene.filters}{a named list with gene filters and
    their parameters. See the main help page of
    \code{\link{metaseqr}} for details.}
    
    \item{sample.list}{the list containing condition names 
    and the samples under each condition.}
}
\value{
    a named list with three members. The first member 
    (\code{result} is a named list whose names are the 
    gene filter names and its members are the filtered 
    rownames of \code{gene.data}. The second member 
    (\code{cutoff} is a named list whose names are the 
    gene filter names and its members are the cutoff 
    values corresponding to each filter. The third member 
    is a matrix of binary flags (0 for non-filtered, 1 
    for filtered) for each gene. The rownames of the 
    flag matrix correspond to gene ids.
}
\description{
    This function performs the gene expression filtering
    based on gene read counts and a set of gene filter rules.
    For more details see the main help pages of
    \code{\link{metaseqr}}.
}
\examples{
\donttest{
data("mm9.gene.data",package="metaseqR")
gene.counts <- mm9.gene.counts
sample.list <- sample.list.mm9
gene.counts <- normalize.edger(as.matrix(gene.counts[,9:12]),
    sample.list)
gene.data <- get.annotation("mm9","gene")
gene.filters <- get.defaults("gene.filter","mm9")
filter.results <- filter.genes(gene.counts,gene.data,
    gene.filters,sample.list)
}
}
\author{
    Panagiotis Moulos
}

