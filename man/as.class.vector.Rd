\name{as.class.vector}
\alias{as.class.vector}
\title{Create a class vector}
\usage{
    as.class.vector(sample.list)
}
\arguments{
    \item{sample.list}{the list containing condition names
    and the samples under each condition.}
}
\value{
    A vector of condition names.
}
\description{
    Creates a class vector from a sample list. Internal to
    the \code{stat.*} functions. Mostly internal use.
}
\examples{
sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
clv <- as.class.vector(sample.list)
}
\author{
    Panagiotis Moulos
}

