\name{mlfo}
\alias{mlfo}
\title{MLE dispersion estimate}
\usage{
    mlfo(phi, y)
}
\arguments{
    \item{phi}{the parameter to be optimized.}

    \item{y}{count samples used to perform 
    the optimization.}
}
\value{
    The objective function value.
}
\description{
    MLE function used to estimate negative 
    binomial dispersions from real RNA-Seq 
    data, as in (Soneson and Delorenzi, BMC 
    Bioinformatics, 2013) and (Robles et al., 
    BMC Genomics, 2012). Internal use.
}
\examples{
\donttest{
# Not yet available
}
}
\author{
    Panagiotis Moulos
}

