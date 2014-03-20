\name{diagplot.noiseq.saturation}
\alias{diagplot.noiseq.saturation}
\title{Simpler implementation of saturation plots inspired from
    NOISeq package}
\usage{
    diagplot.noiseq.saturation(x, o, tb, path = NULL)
}
\arguments{
    \item{x}{the count data matrix.}

    \item{o}{one or more R plotting device to direct the plot
    result to. Supported mechanisms: \code{"x11"} (default),
    \code{"png"}, \code{"jpg"}, \code{"bmp"}, \code{"pdf"} or
    \code{"ps"}.}

    \item{tb}{the vector of biotypes, one for each row of x.}

    \item{path}{the path to create output files.}
}
\value{
    The filenames of the plots produced in a named list with
    names the \code{which.plot} argument. If
    \code{output="x11"}, no output filenames are produced.
}
\description{
    Helper function for \code{\link{diagplot.noiseq}} to plot
    feature detection saturation as presented in the NOISeq
    package vignette. It has two main outputs: a set of
    figures, one for each input sample depicting the
    saturation for each biotype and one single multiplot
    which depicts the saturation of all samples for each
    biotype. It expands the saturation plots of NOISeq by
    allowing more samples to be examined in a simpler way.
    Don't use this function directly. Use either
    \code{\link{diagplot.metaseqr}} or
    \code{\link{diagplot.noiseq}}.
}
\author{
    Panagiotis Moulos
}

