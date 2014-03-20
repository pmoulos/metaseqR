\name{make.sim.data.tcc}
\alias{make.sim.data.tcc}
\title{Create simulated counts using TCC package}
\usage{
    make.sim.data.tcc(...)
}
\arguments{
    \item{...}{parameters to the \code{simulateReadCounts}
    function.}
}
\value{
    A list with the following members: \code{simdata} holding
    the simulated dataset complying with metaseqr
    requirements, and \code{simparam} holding the simulation
    parameters (see TCC documentation). Note that the produced 
    data are based in an Arabidopsis dataset.
}
\description{
    This function creates simulated RNA-Seq gene expression
    datasets using the \code{simulateReadCounts} function
    from the Bioconductor package TCC and it adds simulated
    annoation elements. For further information please
    consult the TCC package documentation.
}
\examples{
\donttest{
dd <- make.sim.data.tcc(Ngene=10000,PDEG=0.2,
        DEG.assign=c(0.9,0.1),
        DEG.foldchange=c(5,5),replicates=c(3,3))
head(dd$simdata)
}
}
\author{
    Panagiotis Moulos
}

