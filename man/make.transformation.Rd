\name{make.transformation}
\alias{make.transformation}
\title{Calculates several transformation of counts}
\usage{
    make.transformation(data.matrix, export.scale, 
        scf = NULL, log.offset = 1)
}
\arguments{
    \item{data.matrix}{the raw or normalized counts matrix.
    Each column represents one input sample.}

    \item{export.scale}{a character vector containing one of
    the supported data transformations (\code{"natural"},
    \code{"log2"}, \code{"log10"},\code{"vst"}). See also the
    main help page of metaseqr.}
    
    \item{scf}{a scaling factor for the reads of each gene, 
    for example the sum of exon lengths or the gene length. 
    Divided by each read count when \code{export.scale="rpgm"}. 
    It provides an RPKM-like measure but not the actual RPKM 
    as this normalization is not supported.}

    \item{log.offset}{a number to be added to each element of
    data.matrix in order to avoid Infinity on log type data
    transformations.}
}
\value{
    A named list whose names are the elements in
    export.scale. Each list member is the respective
    transformed data matrix.
}
\description{
    Returns a list of transformed (normalized) counts, based
    on the input count matrix data.matrix. The data
    transformations are passed from the \code{export.scale}
    parameter and the output list is named accordingly. This
    function is intended mostly for internal use but can also
    be used independently.
}
\examples{
\donttest{
data.matrix <- round(1000*matrix(runif(400),100,4))
rownames(data.matrix) <- paste("gene_",1:100,sep="")
colnames(data.matrix) <- c("C1","C2","T1","T2")
tr <- make.transformation(data.matrix,c("log2","vst"))
head(tr$vst)
}
}
\author{
    Panagiotis Moulos
}

