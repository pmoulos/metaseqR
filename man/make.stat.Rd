\name{make.stat}
\alias{make.stat}
\title{Calculates several statistices on read counts}
\usage{
    make.stat(samples, data.list, stat, export.scale)
}
\arguments{
    \item{samples}{a set of samples from the dataset under
    processing. They should match sample names from
    \code{sample.list}. See also the main help page of
    \code{\link{metaseqr}}.}

    \item{data.list}{a list containing natural or transformed
    data, typically an output 
    from \code{\link{make.transformation}}.}

    \item{stat}{the statistics to calculate. Can be one or
    more of \code{"mean"}, \code{"median"}, \code{"sd"},
    \code{"mad"}, \code{"cv"}, \code{"rcv"}. See also the
    main help page of \code{\link{metaseqr}}.}

    \item{export.scale}{the output transformations used as
    input also to \code{\link{make.transformation}}.}
}
\value{
    A matrix of statistics calculated based on the input
    sample names. The different data transformnations are
    appended columnwise.
}
\description{
    Returns a matrix of statistics calculated for a set of
    given samples. Internal use.
}
\examples{
\donttest{
data.matrix <- round(1000*matrix(runif(400),100,4))
rownames(data.matrix) <- paste("gene_",1:100,sep="")
colnames(data.matrix) <- c("C1","C2","T1","T2")
tr <- make.transformation(data.matrix,c("log2","vst"))
st <- make.stat(c("C1","C2"),tr,c("mean","sd"),c("log2","vst"))
}
}
\author{
    Panagiotis Moulos
}

