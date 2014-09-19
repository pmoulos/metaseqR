\name{get.ucsc.annotation}
\alias{get.ucsc.annotation}
\title{UCSC/RefSeq annotation downloader}
\usage{
    get.ucsc.annotation(org, type, refdb="ucsc",
        multic=FALSE)
}
\arguments{
    \item{org}{the organism for which to download
    annotation.}

    \item{type}{either \code{"gene"} or \code{"exon"}.}

    \item{refdb}{either \code{"ucsc"} or 
    \code{"refseq"}.}

    \item{multic}{a logical value indicating the 
    presence of multiple cores. Defaults to 
    \code{FALSE}. Do not change it if you are not 
    sure whether package parallel has been loaded 
    or not. It is used in the case of 
    \code{type="exon"} to process the return value of 
    the query to the UCSC Genome Browser database.}
}
\value{
    A data frame with the canonical (not isoforms!) genes or
    exons of the requested organism. When
    \code{type="genes"}, the data frame has the following
    columns: chromosome, start, end, gene_id, gc_content,
    strand, gene_name, biotype. When \code{type="exon"} the
    data frame has the following columns: chromosome, start,
    end, exon_id, gene_id, strand, gene_name, biotype. The
    gene_id and exon_id correspond to UCSC or RefSeq gene 
    and exon accessions respectively. The gene_name corresponds 
    to HUGO nomenclature gene names.
}
\description{
    This function connects to the UCSC Genome Browser public 
    database and downloads annotation elements (gene 
    co-ordinates, exon co-ordinates, gene identifications
    etc.) for each of the supported organisms, but using UCSC 
    instead of Ensembl. See the help page of 
    \code{\link{metaseqr}} for a list of supported organisms.
    The function downloads annotation for an organism genes or 
    exons.
}
\note{
    The data frame that is returned contains only "canonical"
    chromosomes for each organism. It does not contain
    haplotypes or random locations and does not contain
    chromosome M.  Note also that as the UCSC databases do 
    not contain biotype classifications like Ensembl. These 
    will be returned as \code{NA} and as a result, some 
    quality control plots will not be available.
}
\examples{
\donttest{
hg19.genes <- get.ucsc.annotation("hg19","gene","ucsc")
mm9.exons <- get.ucsc.annotation("mm9","exon","refseq")
}
}
\author{
    Panagiotis Moulos
}

