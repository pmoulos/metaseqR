\name{wp.adjust}
\alias{wp.adjust}
\title{Multiple testing correction helper}
\usage{
    wp.adjust(p, m)
}
\arguments{
    \item{p}{a vector of p-values.}

    \item{m}{the adjustment method. See the help of
    \code{\link{p.adjust}}.}
}
\description{
    A wrapper around the \code{\link{p.adjust}} function to
    include also the qvalue adjustment procedure from the
    qvalue package. Internal use.
}
\author{
    Panagiotis Moulos
}

