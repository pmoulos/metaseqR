\name{get.strict.biofilter}
\alias{get.strict.biofilter}
\title{Group together a more strict biotype filter}
\usage{
    get.strict.biofilter(org)
}
\arguments{
    \item{org}{one of the supported organisms.}
}
\value{
    A list of booleans, one for each biotype.
}
\description{
    Returns a list with TRUE/FALSE according to the biotypes
    that are going to be filtered in a more strict way than
    the defaults. This is a helper function for the analysis
    presets of metaseqR. Internal use only.
}
\examples{
\donttest{
sf <- get.strict.biofilter("hg18")
}
}
\author{
    Panagiotis Moulos
}

