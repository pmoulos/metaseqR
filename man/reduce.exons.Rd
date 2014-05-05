\name{reduce.exons}
\alias{reduce.exons}
\title{Merges exons to create a unique set of exons for each gene}
\usage{
    reduce.exons(gr, multic = FALSE)
}
\arguments{
    \item{gr}{a GRanges object created from the supplied 
    annotation (see also the \code{\link{read2count}} and 
    \code{\link{get.annotation}} functions.}

    \item{multic}{a logical value indicating the presence 
    of multiple cores. Defaults to \code{FALSE}. Do not 
    change it if you are not sure whether package parallel
    has been loaded or not.}
}
\value{
    A GRanges object with virtual merged exons for each 
    gene/transcript.
}
\description{
    This function uses the \code{"reduce"} function of 
    IRanges to construct virtual unique exons for each 
    gene, so as to avoid inflating the read counts for 
    each gene because of multiple possible transcripts. 
    If the user wants transcripts instead of genes, they 
    should be supplied to the original annotation table.
}
\examples{
\donttest{
    require(GenomicRanges)
    multic <- check.parallel(0.8)
    ann <- get.annotation("mm9","exon")
    gr <- makeGRangesFromDataFrame(
        df=ann,
        keep.extra.columns=TRUE,
        seqnames.field="chromosome"
    )
    re <- reduce.exons(gr,multic=multic)
}
}
\author{
    Panagiotis Moulos
}

