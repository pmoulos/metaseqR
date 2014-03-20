\name{estimate.aufc.weights}
\alias{estimate.aufc.weights}
\title{Estimate AUFC weights}
\usage{
    estimate.aufc.weights(counts, normalization,
        statistics, nsim = 10, N = 10000, 
        samples = c(3, 3), ndeg = c(500, 500),
        top = 500, model.org = "mm9", fc.basis=1.5,
        seed = NULL, draw.fpc = FALSE, multic = FALSE,
        ...)
}
\arguments{
    \item{counts}{the real raw counts table from
    which the simulation parameters will be 
    estimated. It must not be normalized and must
    contain only integer counts, without any other
    annotation elements and unique gene identifiers
    as the rownames attribute.}

    \item{normalization}{same as \code{normalization} 
    in \code{\link{metaseqr}}.}

    \item{statistics}{same as \code{statistics} in 
    \code{\link{metaseqr}}.}

    \item{nsim}{the number of simulations to perform
    to estimate the weights. It default to 10.}

    \item{N}{the number of genes to produce. 
    See \code{\link{make.sim.data.sd}}.}

    \item{samples}{a vector with 2 integers, which
    are the number of samples for each condition 
    (two conditions currently supported).}

    \item{ndeg}{a vector with 2 integers, which
    are the number of differentially expressed 
    genes to be produced. The first element is 
    the number of up-regulated genes while the 
    second is the number of down-regulated genes.}

    \item{fc.basis}{the minimum fold-change for 
    deregulation.}

    \item{top}{the top \code{top} best ranked 
    (according to p-value) to use, to calculate 
    area under the false discovery curve.}

    \item{model.org}{the organism from which the 
    data are derived. It must be one of 
    \code{\link{metaseqr}} supported organisms.}

    \item{seed}{a list of seed for reproducible 
    simulations. Defaults to \code{NULL}.}

    \item{draw.fpc}{draw the averaged false 
    discovery curves? Default to \code{FALSE}.}

    \item{multic}{whether to run in parallel 
    (if package \code{parallel} is present or
    not.}

    \item{...}{Further arguments to be passed to 
    \code{\link{estimate.sim.params}}.}
}
\value{
    A vector of weights to be used in 
    \code{\link{metaseqr}} with the 
    \code{weights} option.
}
\description{
    This function automatically estimates weights 
    for the \code{"weight"} and \code{"dperm.weight"} 
    options of metaseqR for combining p-values from 
    multiple statistical tests. It creates simulated 
    dataset based on real data and then performs 
    statistical analysis with metaseqR several times 
    in order to derive False Discovery Curves. Then, 
    the average areas under the false discovery curves
    are used to construct weights for each algorithm,
    according to its performance when using simulated 
    data.
}
\examples{
\donttest{
data("mm9.gene.data",package="metaseqR")
multic <- check.parallel(0.8)
weights <- estimate.aufc.weights(
    counts=as.matrix(mm9.gene.counts[,9:12]),
    normalization="edaseq",
    statistics=c("deseq","edger"),
    nsim=3,N=100,ndeg=c(10,10),top=10,model.org="mm9",
    seed=10,multic=multic,libsize.gt=1e+5
)
}
}
\author{
    Panagiotis Moulos
}

