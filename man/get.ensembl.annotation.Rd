\name{get.ensembl.annotation}
\alias{get.ensembl.annotation}
\title{Ensembl annotation downloader}
\usage{
    get.ensembl.annotation(org, type)
}
\arguments{
    \item{org}{the organism for which to download
    annotation.}

    \item{type}{either \code{"gene"} or \code{"exon"}.}
}
\value{
    A data frame with the canonical (not isoforms!) genes or
    exons of the requested organism. When
    \code{type="genes"}, the data frame has the following
    columns: chromosome, start, end, gene_id, gc_content,
    strand, gene_name, biotype. When \code{type="exon"} the
    data frame has the following columns: chromosome, start,
    end, exon_id, gene_id, strand, gene_name, biotype. The
    gene_id and exon_id correspond to Ensembl gene and exon
    accessions respectively. The gene_name corresponds to
    HUGO nomenclature gene names.
}
\description{
    This function connects to the EBI's Biomart service using
    the package biomaRt and downloads annotation elements
    (gene co-ordinates, exon co-ordinates, gene
    identifications, biotypes etc.) for each of the supported
    organisms. See the help page of \code{\link{metaseqr}}
    for a list of supported organisms. The function downloads
    annotation for an organism genes or exons.
}
\note{
    The data frame that is returned contains only "canonical"
    chromosomes for each organism. It does not contain
    haplotypes or random locations and does not contain
    chromosome M.
}
\examples{
\donttest{
hg19.genes <- get.ensembl.annotation("hg19","gene")
mm9.exons <- get.ensembl.annotation("mm9","exon")
}
}
\author{
    Panagiotis Moulos
}

