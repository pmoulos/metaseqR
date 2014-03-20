\name{check.libsize}
\alias{check.libsize}
\title{Library size validator}
\usage{
    check.libsize(libsize.list, sample.list)
}
\arguments{
    \item{libsize.list}{the samples-names library size list.}

    \item{sample.list}{the input sample list.}
}
\description{
    Checks the names of the supplied library sizes. Internal
    use only.
}
\examples{
sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
libsize.list.1 <- list(A1=1e+6,A2=1.1e+6,B1=1.2e+6,
B2=1.3e+6,B3=1.5e+6)
libsize.list.2 <- list(A1=1e+6,A2=1.1e+6,B1=1.2e+6,
B2=1.3e+6)
check.libsize(libsize.list.1,sample.list) # Will work
#check.libsize(libsize.list.2,sample.list) # Will throw error!
}
\author{
    Panagiotis Moulos
}

