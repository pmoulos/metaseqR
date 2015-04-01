\name{meta.test}
\alias{meta.test}
\title{Meta-analysis using several RNA-Seq statistics}
\usage{
    meta.test(cp.list,
        meta.p = c("simes", "bonferroni", "fisher", "dperm.min", 
        "dperm.max", "dperm.weight", "fperm", "whitlock", 
        "minp", "maxp", "weight", "pandora", "none"), counts, 
        sample.list, statistics, stat.args, libsize.list, 
        nperm = 10000, weight = rep(1/length(statistics), 
        length(statistics)), reprod=TRUE, multic = FALSE)
}
\arguments{
    \item{cp.list}{a named list whose names are the contrasts
    requested from metaseqr. Each member is a p-value matrix
    whose colnames are the names of the statistical tests
    applied to the data. See the main \code{\link{metaseqr}}
    help page.}

    \item{meta.p}{the p-value combination method to use. See
    the main \code{\link{metaseqr}} help page.}

    \item{counts}{the normalized and possibly filtered read
    counts matrix. See the main \code{\link{metaseqr}} help
    page.}

    \item{sample.list}{the list containing condition names
    and the samples under each condition. See the main
    \code{\link{metaseqr}} help page.}

    \item{statistics}{the statistical algorithms used in
    metaseqr. See the main \code{\link{metaseqr}} help page.}

    \item{stat.args}{the parameters for each statistical
    argument. See the main \code{\link{metaseqr}} help page.}

    \item{libsize.list}{a list with library sizes. See the
    main \code{\link{metaseqr}} and the \code{stat.*} help
    pages.}

    \item{nperm}{the number of permutations (Monte Carlo
    simulations) to perform.}

    \item{weight}{a numeric vector of weights for each
    statistical algorithm.}

    \item{reprod}{create reproducible permutations when 
    \code{meta.p="dperm.min"}, \code{meta.p="dperm.max"} 
    or \code{meta.p="dperm.weight"}. Ideally one would 
    want to create the same set of indices for a given 
    dataset so as to create reproducible p-values. If 
    \code{reprod=TRUE}, a fixed seed is used by 
    \code{meta.perm} for all the datasets analyzed with
    \code{metaseqr}. If \code{reprod=FALSE}, then the 
    p-values will not be reproducible, although 
    statistical significance is not expected to change 
    for a large number of resambling. Finally, 
    \code{reprod} can be a numeric vector of seeds with 
    the same length as \code{nperm} so that the user can 
    supply his/her own seeds.}

    \item{multic}{use multiple cores to execute the
    premutations. This is an external parameter and implies
    the existence of multicore package in the execution
    environment. See the main \code{\link{metaseqr}} help
    page.}
}
\value{
    A named list with combined p-values. The names are the
    contrasts and the list members are combined p-value
    vectors, one for each contrast.
}
\description{
    This function calculates the combined p-values when
    multiple statistical algorithms are applied to the input
    dataset. It is a helper and it requires very specific
    arguments so it should not be used individually
}
\examples{
\donttest{
# Not yet available
}
}
\author{
    Panagiotis Moulos
}

