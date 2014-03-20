\name{make.sim.data.sd}
\alias{make.sim.data.sd}
\title{Create simulated counts using the 
    Soneson-Delorenzi method}
\usage{
    make.sim.data.sd(N, param, samples = c(5, 5),
        ndeg = rep(round(0.1*N), 2), fc.basis = 1.5,
    libsize.range = c(0.7, 1.4), libsize.mag = 1e+7,
        model.org = NULL, sim.length.bias = FALSE, 
        seed = NULL)
}
\arguments{
    \item{N}{the number of genes to produce.}

    \item{param}{a named list with negative binomial 
    parameter sets to sample from. The first member is
    the mean parameter to sample from (\code{mu.hat}) 
    and the second the dispersion (\code{phi.hat}). 
    This list can be created with the 
    \code{\link{estimate.sim.params}} function.}

    \item{samples}{a vector with 2 integers, 
    which are the number of samples for each 
    condition (two conditions currently supported).}

    \item{ndeg}{a vector with 2 integers, which are 
    the number of differentially expressed genes to 
    be produced. The first element is the number of 
    up-regulated genes while the second is the 
    number of down-regulated genes.}

    \item{fc.basis}{the minimum fold-change for 
    deregulation.}

    \item{libsize.range}{a vector with 2 numbers 
    (generally small, see the default), as they 
    are multiplied with \code{libsize.mag}. These 
    numbers control the library sized of the 
    synthetic data to be produced.}

    \item{libsize.mag}{a (big) number to multiply 
    the \code{libsize.range} to produce library 
    sizes.}

    \item{model.org}{the organism from which the 
    real data are derived from. It must be one 
    of the supported organisms (see the main 
    \code{\link{metaseqr}} help page). It is used 
    to sample real values for GC content.}

    \item{sim.length.bias}{a boolean to instruct 
    the simulator to create genes whose read counts is
    proportional to their length. This is achieved by 
    sorting in increasing order the mean parameter of 
    the negative binomial distribution (and the 
    dispersion according to the mean) which will cause 
    an increasing gene count length with the sampling. 
    The sampled lengths are also sorted so that in the 
    final gene list, shorter genes have less counts as 
    compared to the longer ones. The default is FALSE.}

    \item{seed}{a seed to use with random number 
    generation for reproducibility.}
}
\value{
    A named list with two members. The first 
    member (\code{simdata}) contains the 
    synthetic dataset 
}
\description{
    This function creates simulated RNA-Seq gene 
    expression datasets using the method presented 
    in (Soneson and Delorenzi, BMC Bioinformatics, 
    2013). For the time being, it creates only 
    simulated datasets with two conditions.
}
\examples{
\donttest{
# File "bottomly_read_counts.txt" from the ReCount database
download.file(paste("http://bowtie-bio.sourceforge.net/recount/",
    "countTables/bottomly_count_table.txt",sep=""),
    destfile="~/bottomly_count_table.txt")
N <- 10000
par.list <- estimate.sim.params("~/bottomly_read_counts.txt")
sim <- make.sim.data.sd(N,par.list)
synth.data <- sim$simdata
true.deg <- which(sim$truedeg!=0)
}
}
\author{
    Panagiotis Moulos
}

