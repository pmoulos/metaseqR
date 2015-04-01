\name{build.export}
\alias{build.export}
\title{Results export builder}
\usage{
    build.export(gene.data, raw.gene.counts,
        norm.gene.counts, flags, sample.list, cnt, 
        statistics, raw.list, norm.list,
        p.mat = matrix(NA, nrow(gene.data), length(statistics)),
        adj.p.mat = matrix(NA, nrow(gene.data), length(statistics)),
        sum.p = rep(NA, nrow(gene.data)),
        adj.sum.p = rep(NA, nrow(gene.data)),
        export.what = c("annotation", "p.value", "adj.p.value", "meta.p.value", 
                "adj.meta.p.value", "fold.change", "stats", "counts","flags"),
        export.scale = c("natural", "log2", "log10", "rpgm", "vst"),
        export.values = c("raw", "normalized"),
        export.stats = c("mean", "median", "sd", "mad", "cv", "rcv"),
        log.offset = 1, report = TRUE)
}
\arguments{
    \item{gene.data}{an annotation data frame (such the ones
    produced by \code{\link{get.annotation}}).}

    \item{raw.gene.counts}{a matrix of filtering flags (0,1), created 
    by the filtering functions.}

    \item{norm.gene.counts}{a matrix of normalized gene
    counts.}

    \item{flags}{a matrix of normalized gene
    counts.}

    \item{sample.list}{see the documentation of
    \code{\link{metaseqr}}.}

    \item{cnt}{the statistical contrast for which the export
    builder is currently running.}

    \item{statistics}{the statistical tests used (see the
    documentation of \code{\link{metaseqr}}).}

    \item{raw.list}{a list of transformed un-normalized
    counts, see the documentation 
    of \code{\link{make.transformation}}.}

    \item{norm.list}{a list of transformed normalized counts,
    see the documentation 
    of \code{\link{make.transformation}}.}

    \item{p.mat}{a matrix of p-values, see the documentation
    of \code{\link{metaseqr}}.}

    \item{adj.p.mat}{a matrix of adjusted p-values, see the
    documentation of \code{\link{metaseqr}}.}

    \item{sum.p}{a vector of combined p-values, see the
    documentation of \code{\link{metaseqr}}.}

    \item{adj.sum.p}{a vector of adjusted combined p-values,
    see the documentation of \code{\link{metaseqr}}.}

    \item{export.what}{see the documentation of
    \code{\link{metaseqr}}.}

    \item{export.scale}{see the documentation of
    \code{\link{metaseqr}}.}

    \item{export.values}{see the documentation of
    \code{\link{metaseqr}}.}

    \item{export.stats}{see the documentation of
    \code{\link{metaseqr}}.}

    \item{log.offset}{see the documentation of
    \code{\link{metaseqr}}.}

    \item{report}{see the documentation of
    \code{\link{metaseqr}}.}
}
\value{
    A list with three members: a data frame to be exported in
    a text file, a long string with the result in a html
    formatted table (if \code{report=TRUE}) and the column
    names of the output data frame.
}
\description{
    This function help build the output files of the metaseqr
    pipeline based on several elements produced during the
    pipeline execution. It is intended for internal use and
    not available to the users.
}
\examples{
\dontrun{
# Not yet available
}
}
\author{
    Panagiotis Moulos
}

