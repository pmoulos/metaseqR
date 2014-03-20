\name{diagplot.edaseq}
\alias{diagplot.edaseq}
\title{Diagnostic plots based on the EDASeq package}
\usage{
    diagplot.edaseq(x, sample.list, covar = NULL,
        is.norm = FALSE,
        which.plot = c("meanvar", "meandiff", "gcbias", "lengthbias"),
        output = "x11", path = NULL, ...)
}
\arguments{
    \item{x}{the count data matrix.}

    \item{sample.list}{the list containing condition names
    and the samples under each condition.}

    \item{covar}{The covariate to plot counts against.
    Usually \code{"gc"} or \code{"length"}.}

    \item{is.norm}{a logical indicating whether object
    contains raw or normalized data. It is not essential and
    it serves only plot annotation purposes.}

    \item{which.plot}{the EDASeq package plot to generate. It
    can be one or more of \code{"meanvar"},
    \code{"meandiff"}, \code{"gcbias"} or
    \code{"lengthbias"}. Please refer to the documentation of
    the EDASeq package for details on the use of these 
    plots. The \code{which.plot="lengthbias"} case is 
    not covered by EDASeq documentation, however it is 
    similar to the GC-bias plot when the covariate is the 
    gene length instead of the GC content.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filenames of the plot produced in a named list with
    names the which.plot argument. If \code{output="x11"}, no
    output filenames are produced.
}
\description{
    A wrapper around the plotting functions availale in the
    EDASeq normalization Bioconductor package. For analytical
    explanation of each plot please see the vignette of the
    EDASeq package. It is best to use this function through
    the main plotting function
    \code{\link{diagplot.metaseqr}}.
}
\examples{
\donttest{
require(DESeq)
data.matrix <- counts(makeExampleCountDataSet())
sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
diagplot.edaseq(data.matrix,sample.list,which.plot="meandiff")
}
}
\author{
    Panagiotis Moulos
}

