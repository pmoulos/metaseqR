\name{make.avg.expression}
\alias{make.avg.expression}
\title{Calculates fold changes}
\usage{
    make.avg.expression(contrast, sample.list, data.matrix,
        log.offset = 1)
}
\arguments{
    \item{contrast}{the vector of requested statistical
    comparison contrasts.}

    \item{sample.list}{the list containing condition names
    and the samples under each condition.}

    \item{data.matrix}{a matrix of gene expression data whose
    column names are the same as the sample names included in
    the sample list.}

    \item{log.offset}{a number to be added to each element of
    data matrix in order to avoid Infinity on log type data
    transformations.}
}
\value{
    A matrix of fold change ratios, treatment to control, as
    these are parsed from contrast.
}
\description{
    Returns a matrix of average expressions (A in MA plot) based 
    on the requested contrast, the list of all samples and the 
    data matrix which is produced by the metaseqr workflow. For 
    details on the \code{contrast}, \code{sample.list} and 
    \code{log.offset} parameters, see the main usage page of metaseqr. 
 This function is intended mostly for internal use but can also be used 
 independently.
}
\examples{
\donttest{
data.matrix <- round(1000*matrix(runif(400),100,4))
rownames(data.matrix) <- paste("gene_",1:100,sep="")
colnames(data.matrix) <- c("C1","C2","T1","T2")
a <- make.avg.expression("Control_vs_Treatment",list(Control=c("C1","C2"),
    Treatment=c("T1","T2")),data.matrix)
}
}
\author{
    Panagiotis Moulos
}

