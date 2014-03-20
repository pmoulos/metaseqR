\name{check.packages}
\alias{check.packages}
\title{Required packages validator}
\usage{
    check.packages(m, p)
}
\arguments{
    \item{m}{meta-analysis method.}
    
    \item{p}{QC plot types.}
}
\description{
    Checks if all the any required packages, not attached 
    during installation or loading, are present according 
    to metaseqR input options. Internal use only.
}
\examples{
check.packages(c("simes","whitlock"),
    c("gcbias","correl"))
}
\author{
    Panagiotis Moulos
}
