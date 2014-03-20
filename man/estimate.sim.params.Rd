\name{estimate.sim.params}
\alias{estimate.sim.params}
\title{Estimate negative binomial parameters from 
    real data}
\usage{
    estimate.sim.params(real.counts, libsize.gt = 3e+6,
        rowmeans.gt = 5,eps = 1e-11, 
        restrict.cores = 0.8, seed = 42, draw = FALSE)
}
\arguments{
    \item{real.counts}{a text tab-delimited file 
    with real RNA-Seq data. The file should 
    strictly contain a unique gene name (e.g. Ensembl
    accession) in the first column and all other 
    columns should contain read counts for each 
    gene. Each column must be named with a unique 
    sample identifier. See examples in the ReCount 
    database 
    \url{http://bowtie-bio.sourceforge.net/recount/}.}

    \item{libsize.gt}{a library size below which 
    samples are excluded from parameter estimation 
    (default: 3000000).}

    \item{rowmeans.gt}{a row means (mean counts 
    over samples for each gene) below which 
    genes are excluded from parameter estimation
    (default: 5).}

    \item{eps}{the tolerance for the convergence 
    of \code{\link{optimize}} function. Defaults 
    to 1e-11.}

    \item{restrict.cores}{in case of parallel 
    optimization, the fraction of the available 
    cores to use.}

    \item{seed}{a seed to use with random number 
    generation for reproducibility.}

    \item{draw}{boolean to determine whether to 
    plot the estimated simulation parameters 
    (mean and dispersion) or not. Defaults to 
    \code{FALSE} (do not draw a mean-dispersion 
    scatterplot).}
}
\value{
    A named list with two members: \code{mu.hat}
    which contains negative binomial mean 
    estimates and \code{phi.hat} which contains 
    dispersion estimates.
}
\description{
    This function reads a read counts table 
    containing real RNA-Seq data (preferebly 
    with more than 20 samples so as to get as 
    much accurate as possible estimations) and 
    calculates a population of count means and 
    dispersion parameters which can be used to 
    simulate an RNA-Seq dataset with synthetic 
    genes by drawing from a negative binomial 
    distribution. This function works in the 
    same way as described in (Soneson and 
    Delorenzi, BMC Bioinformatics, 2013) and 
    (Robles et al., BMC Genomics, 2012).
}
\examples{
\donttest{
# Dowload locally the file "bottomly_read_counts.txt" from
# the ReCount database
download.file(paste("http://bowtie-bio.sourceforge.net/",
    "recount/countTables/bottomly_count_table.txt",sep=""),
    destfile="~/bottomly_count_table.txt")
# Estimate simulation parameters
par.list <- estimate.sim.params("~/bottomly_count_table.txt")
}
}
\author{
    Panagiotis Moulos
}

