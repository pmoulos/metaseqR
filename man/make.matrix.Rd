\name{make.matrix}
\alias{make.matrix}
\title{Results output build helper}
\usage{
    make.matrix(samples, data.list, export.scale = "natural")
}
\arguments{
    \item{samples}{a set of samples from the dataset under
    processing. They should match sample names from
    \code{sample.list}. See also the main help page of
    \code{\link{metaseqr}}.}

    \item{data.list}{a list containing natural or transformed
    data, typically an output 
    from \code{\link{make.transformation}}.}

    \item{export.scale}{the output transformations used as
    input also to \code{\link{make.transformation}}.}
}
\value{
    A named list whose names are the elements in
    \code{export.scale}. Each list member is the respective
    sample subest data matrix.
}
\description{
    Returns a list of matrices based on the export scales
    that have been chosen from the main function and a subset
    of samples based on the sample names provided in the
    \code{sample.list} argument of the main
    \code{\link{metaseqr}} function. Internal use.
}
\examples{
\donttest{
data.matrix <- round(1000*matrix(runif(400),100,4))
rownames(data.matrix) <- paste("gene_",1:100,sep="")
colnames(data.matrix) <- c("C1","C2","T1","T2")
tr <- make.transformation(data.matrix,c("log2","vst"))
mm <- make.matrix(c("C1","T1"),tr,"log2")
head(tr$vst)
}
}
\author{
    Panagiotis Moulos
}

