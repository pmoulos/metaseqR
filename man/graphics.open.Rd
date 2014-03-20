\name{graphics.open}
\alias{graphics.open}
\title{Open plotting device}
\usage{
    graphics.open(o, f, ...)
}
\arguments{
    \item{o}{the plotting device, see main metaseqr function}

    \item{f}{a filename, if the plotting device requires it
    (e.g. \code{"pdf"})}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\description{
    Wrapper function to open a plotting device. Internal use
    only.
}
\examples{
\donttest{
graphics.open("pdf","test.pdf",width=12,height=12)
}
}
\author{
    Panagiotis Moulos
}

