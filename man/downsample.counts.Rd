\name{downsample.counts}
\alias{downsample.counts}
\title{Downsample read counts}
\usage{
    downsample.counts(counts,seed=42)
}
\arguments{
    \item{counts}{the read counts table 
    which is subjected to downsampling.}

    \item{seed}{random seed for reproducible 
    downsampling.}
}
\value{
    The downsampled counts matrix.
}
\description{
    This function downsamples the library sizes 
    of a read counts table to the lowest library 
    size, according to the methdology used in 
    (Soneson and Delorenzi, BMC Bioinformatics, 
    2013).
}
\examples{
\donttest{
# Dowload locally the file "bottomly_count_table.txt" from
# the ReCount database
download.file(paste("http://bowtie-bio.sourceforge.net/",
    "recount/countTables/bottomly_count_table.txt",sep=""),
    destfile="~/bottomly_count_table.txt")
M <- as.matrix(read.delim("~/bottomly_count_table.txt",row.names=1))
D <- downsample.counts(M)
}
}
\author{
    Panagiotis Moulos
}

