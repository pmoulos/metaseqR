\name{meta.perm}
\alias{meta.perm}
\title{Permutation tests for meta-analysis}
\usage{
    meta.perm(contrast, counts, sample.list, statistics,
        stat.args, libsize.list, nperm = 10000,
        weight = rep(1/ncol(counts), ncol(counts)),
        select = c("min", "max", "weight"), replace = "auto",
        reprod=TRUE, multic = FALSE)
}
\arguments{
    \item{contrast}{the contrasts to be tested by each
    statistical algorithm. See the main
    \code{\link{metaseqr}} help page.}

    \item{counts}{a normalized read counts table, one row for
    each gene, one column for each sample.}

    \item{sample.list}{the list containing condition names
    and the samples under each condition. See the main
    \code{\link{metaseqr}} help page.}

    \item{statistics}{the statistical algorithms used in
    metaseqr. See the main \code{\link{metaseqr}} help page.}

    \item{stat.args}{the parameters for each statistical
    algorithm. See the main \code{\link{metaseqr}} help page.}

    \item{libsize.list}{a list with library sizes. See the
    main \code{\link{metaseqr}} and the \code{stat.*} help
    pages.}

    \item{nperm}{the number of permutations (Monte Carlo
    simulations) to perform.}

    \item{weight}{a numeric vector of weights for each
    statistical algorithm.}

    \item{select}{how to select the initial vector of
    p-values. It can be \code{"min"} to select the minimum
    p-value for each gene (more conservative), \code{"max"}
    to select the maximum p-value for each gene (less
    conservative), \code{"weight"} to apply the weights to
    the p-value vector for each gene and derive a weighted
    p-value.}

    \item{replace}{same as the \code{replace} argument in the
    \code{\link{sample}} function. Implies bootstraping or
    simple resampling without replacement. It can also be
    \code{"auto"}, to determine bootstraping or not with the
    following rule: if \code{ncol(counts)<=6}
    \code{replace=FALSE else} \code{replace=TRUE}. This
    protects from the case of having zero variability across
    resampled conditions. In such cases, most statistical
    tests would crash.}

    \item{reprod}{create reproducible permutations. Ideally 
    one would want to create the same set of indices for a 
    given dataset so as to create reproducible p-values. If 
    \code{reprod=TRUE}, a fixed seed is used by 
    \code{meta.perm} for all the datasets analyzed
    with \code{metaseqr}. If \code{reprod=FALSE}, then the 
    p-values will not be reproducible, although statistical 
    significance is not expected to change for a large 
    number of resambling. Finally, \code{reprod} can be a 
    numeric vector of seeds with the same length as 
    \code{nperm} so that the user can supply his/her own 
    seeds.}

    \item{multic}{use multiple cores to execute the
    premutations. This is an external parameter and implies
    the existence of multicore package in the execution
    environment. See the main \code{\link{metaseqr}} help
    page.}
}
\value{
    A vector of meta p-values
}
\description{
    This function performs permutation tests in order to
    derive a meta p-value by combining several of the
    statistical algorithms of metaseqr. This is probably the
    most accurate way of combining multiple statistical
    algorithms for RNA-Seq data, as this issue is different
    from the classic interpretation of the term
    "meta-analysis" which implies the application of the same
    statistical test on different datasets treating the same
    subject/experiment. For other methods, see also the main
    \code{\link{metaseqr}} help page. You should keep in mind
    that the permutation procedure can take a long time, even
    when executed in parallel.
}
\examples{
\donttest{
# Not yet available
}
}
\author{
    Panagiotis Moulos
}

