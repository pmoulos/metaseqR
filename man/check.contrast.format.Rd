\name{check.contrast.format}
\alias{check.contrast.format}
\title{Contrast validator}
\usage{
    check.contrast.format(cnt, sample.list)
}
\arguments{
    \item{cnt}{contrasts vector.}

    \item{sample.list}{the input sample list.}
}
\description{
    Checks if the contrast vector follows the specified
    format. Internal use only.
}
\examples{
sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
cnt <- c("A_vs_B") # Will work
#cnt <- c("A_vs_C") ## Will throw error!
check.contrast.format(cnt,sample.list)
}
\author{
    Panagiotis Moulos
}

