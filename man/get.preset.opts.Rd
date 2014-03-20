\name{get.preset.opts}
\alias{get.preset.opts}
\title{Return several analysis options given an analysis preset}
\usage{
    get.preset.opts(preset, org)
}
\arguments{
    \item{preset}{preset can be one of \code{"all.basic"},
    \code{"all.normal"}, \code{"all.full"},
    \code{"medium.basic"}, \code{"medium.normal"},}

    \item{org}{one of the supported organisms. See
    \code{\link{metaseqr}} main help page.
    \code{"medium.full"}, \code{"strict.basic"},
    \code{"strict.normal"} or \code{"strict.full"}, each of
    which control the strictness of the analysis and the
    amount of data to be exported. For an explanation of the
    presets, see the main \code{\link{metaseqr}} help page.}
}
\value{
    A named list with names \code{exon.filters},
    \code{gene.filters}, \code{pcut}, \code{export.what},
    \code{export.scale}, \code{export.values} and
    \code{export.stats}, each of which correspond to an
    element of the metaseqr pipeline.
}
\description{
    This is a helper function which returns a set of metaseqr
    pipeline options, grouped together according to a preset
    keyword. It is intended mostly for internal use.
}
\examples{
strict.preset <- get.preset.opts("strict.basic","mm9")
}
\author{
    Panagiotis Moulos
}

